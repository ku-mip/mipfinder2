#include <cassert>
#include <concepts>
#include <cmath>
#include <filesystem>
#include <functional>
#include <numeric>
#include <set>
#include <unordered_set>


#include "aliases.h"
#include "ancestor.h"
#include "configuration.h"
#include "easylogging++.h"
#include "fasta.h"
#include "file.h"
#include "go.h"
#include "hmmer.h"
#include "helpers.h"
#include "interpro.h"
#include "mipfinder.h"
#include "protein.h"
#include "proteome.h"
#include "printer.h"

#include <ranges>

using ProteinList = std::vector<Protein>;

INITIALIZE_EASYLOGGINGPP

/* Helper functions */
namespace
{
    /* Ensures that all files required by mipfinder are found. If some
     * dependencies are missing, throws a std::runtime_error */
    void checkFileDependencies(mipfinder::Mipfinder::FileParamaters file_parameters)
    {
        mipfinder::file::exists(file_parameters.input_proteome);
        mipfinder::file::exists(file_parameters.known_microprotein_list);
        mipfinder::file::exists(file_parameters.gene_ontology_database);
        mipfinder::file::exists(file_parameters.interpro_database);
        mipfinder::file::exists(file_parameters.uniprot_to_intepro_id_conversion_file);
        mipfinder::file::exists(file_parameters.uniprot_to_go_id_conversion_file);
    }

    void configureLogger()
    {
        el::Configurations default_logging;
        default_logging.setToDefault();
        default_logging.set(el::Level::Debug,
            el::ConfigurationType::Enabled,
            "true");
        default_logging.set(el::Level::Debug,
            el::ConfigurationType::Filename,
            "debug.log");
        default_logging.set(el::Level::Info,
            el::ConfigurationType::Filename,
            "mipfinder.log");
        default_logging.set(el::Level::Info,
            el::ConfigurationType::Format,
            "%datetime %level %msg");

        el::Loggers::reconfigureLogger("default", default_logging);
    }

    double scoring_algorithm(const mipfinder::HmmerScoringData& data)
    {
        std::size_t cmip_length = data.query->length();
        std::size_t ancestor_length = data.target->length();
        double cmip_coverage = static_cast<double>(cmip_length / ancestor_length);

        double growth_rate = 0.1;
        double midpoint = 50.0;
        double weight = 1 / (1 + std::exp(growth_rate * (cmip_coverage * 100 - midpoint)));

        double assigned_score = (1.0 / std::pow(data.bitscore, cmip_coverage))
            * data.bitscore
            * weight;
        return assigned_score;
    }

     //mipfinder::homology::Results
     //	filterSingleDomainAncestors(const mipfinder::Proteome& proteome,
     //								const mipfinder::homology::Results& results)
     //{
     //	mipfinder::homology::Results filtered;

     //	for (const auto& result : results) {
     //		const auto& ancestor_id = result.target;
     //		const auto ancestor = proteome.find(ancestor_id);
     //		assert(ancestor != nullptr);

     //		unsigned domain_count{0};
     //		for (const auto& interpro_entry : ancestor->interproEntries()) {
     //			if (interpro_entry.type == mipfinder::Interpro::Type::DOMAIN_TYPE ||
     //				interpro_entry.type == mipfinder::Interpro::Type::REPEAT) {
     //				++domain_count;
     //			}
     //		}

     //		/* If it is 0, it is very likely that the annotation of that gene has
     //		* failed since virtually all proteins are made up of small domains.
     //		* If it is exactly 1, we assume that the annotation is correct and it
     //		* really only does have one domain */
     //		if (domain_count != 1) {
     //			filtered.push_back(result);
     //		}
     //	}
     //	return filtered;
     //}
}


namespace detail
{
    //Take all hmmprofile files (ending in ".hmmprofile") in "hmmprofile_directory" and create one file with all the profiles as "output_file"
    void mergeHmmprofileFiles(const std::filesystem::path& hmmprofile_directory, const std::filesystem::path& output_file)
    {
        typedef std::filesystem::directory_iterator DirectoryIter;

        std::ofstream of;
        of.open(output_file, std::ios::trunc);
        for (const auto& directory_entry : DirectoryIter{hmmprofile_directory}) {
            //If the output file is in the same directory as the hmmprofile files, 
            //do not also process the output file! Otherwise it'll be an infinite loop.
            if (directory_entry.path() == output_file) { 
                continue;
            }
            
            LOG(DEBUG) << "Processing " << directory_entry.path();
            if (directory_entry.path().extension() != std::filesystem::path{".hmmprofile"}) {
                continue;
            }
            std::ifstream f;
            f.open(directory_entry.path());
            std::string line;
            while (getline(f, line)) {
                of << line << "\n";
            }
        }
    }

    //Return a collection of proteins from a FASTA file
    template <typename T = mipfinder::ProteinList>
    mipfinder::ProteinList loadProteome(const std::filesystem::path& fasta_file)
    {
        LOG(DEBUG) << "Loading proteome...";
        auto file = mipfinder::file::open(fasta_file);
        const mipfinder::FastaRecords proteome_fasta_records = mipfinder::fasta::extractRecords(file);

        T proteome{};
        const char separator = '_';
        for (const auto& [header, sequence] : proteome_fasta_records) {
            const auto& [protein_id, sequence_version, description, existence_level] = mipfinder::fasta::extractUniprotHeader(header);
            const std::string identifier = protein_id + separator + sequence_version;

            if constexpr (EmplaceableContainer<mipfinder::ProteinList>) {
                proteome.emplace_back(mipfinder::Protein{identifier, sequence, description, std::stoi(existence_level)});
            }
            else if constexpr (InsertableContainer<mipfinder::ProteinList>) {
                proteome.insert(mipfinder::Protein{identifier, sequence, description, std::stoi(existence_level)});
            }            
        }
        LOG(DEBUG) << "Found " << proteome.size() << " proteins";
        return proteome;
    }

    //Compare microproteins to each other to establish which ones have homologues (homologous microproteins)
    //and which ones do not (single-copy microproteins).
    template <typename T>
    void compareMicroproteinsToMicroproteins(const T& potential_microproteins,
                                             const mipfinder::Mipfinder::HmmerParameters& parameters,
                                             const std::filesystem::path& results_output)
    {
        LOG(DEBUG) << "Finding homologous relationship between potential microproteins";
        const auto extra_param = "--mx " + parameters.scoring_matrix;

        const std::filesystem::path results_path = results_output.parent_path();
        const std::filesystem::path query_file_location = results_path / "microprotein_query.fasta";
        mipfinder::proteinToFasta(potential_microproteins, query_file_location);

        mipfinder::homology::phmmer(query_file_location, query_file_location, results_output, extra_param);
        LOG(DEBUG) << "Finished finding microProtein homologues";
    }

    //Filter out proteins whose properties suggest that they are not really present within living cells.
    //
    //@Params
    //proteome - A container of proteins.
    //maximum_allowed_existence_level - Cutoff value for existence level as defined by UniProt
    //
    //@Return - A container whose members are the elements that meet the criteria for real proteins. 
    //          If @proteome is empty, return an empty container.
    template <typename T>
    requires std::ranges::range<T>&& requires (std::ranges::range_value_t<T> v)
    {
        { v.existenceLevel() } -> std::convertible_to<std::size_t>;
    }
    mipfinder::ProteinList removeSpuriousProteins(const T& proteome, const std::size_t maximum_allowed_existence_level)
    {
        LOG(DEBUG) << "Removing proteins with existence level equal to or less than " << maximum_allowed_existence_level;
        auto protein_existence_filter = [=](const auto& protein)
        {
            return (protein.existenceLevel() <= maximum_allowed_existence_level);
        };
        return detail::toContainer<std::unordered_set>(proteome | std::views::filter(protein_existence_filter));
    }

    using HomologyTable = std::unordered_map<std::string, std::unordered_set<std::string>>;

    template <typename T, typename U>
    struct ClassifiedMicroproteins
    {
        T single_copy;
        U homologous;
        HomologyTable homology_table;
    };

    //Create a homology relationship table based on homology search results.
    //
    //@Params
    //homology_search_results - A container of homology search results.
    //
    //@Return - A associative table where each key denotes the identifier of a protein and each
    //          value is a container of homologous protein identifiers of that key. If
    //          @homology_search_results is empty, return an empty table.
    template <typename T>
    requires std::ranges::range<T>&& requires (std::ranges::range_value_t<T> v)
    {
        { v.query } -> std::convertible_to<std::string>;
        { v.target } -> std::convertible_to<std::string>;
    }
    HomologyTable createHomologyTable(T&& homology_search_results)
    {
        LOG(DEBUG) << "Creating a homology relatonship table from homology search results";
        HomologyTable homology_table; //Keys are protein identifiers, values are homologous protein identifiers

        //Associate potential microproteins with their homologues
        for (const auto& result : homology_search_results) {
            //Ignore self-hits, proteins aren't homologous to themselves
            if (result.query == result.target) {
                continue;
            }

            if (!homology_table.contains(result.query)) {
                homology_table[result.query] = std::unordered_set{ result.target };
            }
            else {
                homology_table[result.query].insert(result.target);
            }
        }
        LOG(DEBUG) << "Finished creating a homology table";
        return homology_table;
    }





    //Classify microProteins into single-copy and homologous microproteins based on the homology search results.
    //Single-copy microproteins are proteins that do not have any homologues among other microproteins, while
    //homologous microproteins are proteins that have at least one other homologue among microproteins.
    //
    //@Params
    //microproteins - A container of microproteins to be classified.
    //homology_search_results - A container of homology search results corresponding to the given @microproteins.
    //
    //@Return - Two lists corresponding to single-copy and homologous microproteins. If the @homology_search_results
    //          were not obtained from comparing the @microproteins, the result is undefined. 
    template <typename T, typename U>
    requires std::ranges::range<T>&& std::ranges::range<U>
        && requires (std::ranges::range_value_t<T> t, std::ranges::range_value_t<U> u)
    {
        { t.identifier() } -> std::convertible_to<std::string>;
        { u.query } -> std::convertible_to<std::string>;
        { u.target } -> std::convertible_to<std::string>;
    }
    ClassifiedMicroproteins<mipfinder::ProteinList, mipfinder::ProteinList> classifyMicroproteins(const T& microproteins, const U& homology_search_results)
    {
        LOG(DEBUG) << "Classifying microProteins into single-copy and homologous";
        //Find which potential microproteins have homologues (homologous microProteins) and
        //which do not (single-copy)
        auto homology_table = detail::createHomologyTable(homology_search_results);

        //Ignore potential microproteins for which there is not a corresponding homology
        //query. This can happen when the homology search results were conducted on a 
        //different set of proteins than the one supplied as potential microproteins.
        auto protein_list_contains_homology_result = [&homology_table](const auto& elem)
        {
            return homology_table.contains(elem.identifier());
        };
        auto compared_microproteins = microproteins | std::views::filter(protein_list_contains_homology_result);

        auto is_single_copy_microprotein = [&homology_table](const auto& elem)
        {
            return homology_table.at(elem.identifier()).size() == 0;
        };

        auto unique_microproteins = detail::toContainer<std::unordered_set>(compared_microproteins | std::views::filter(is_single_copy_microprotein));
        auto homologous_microproteins = detail::toContainer<std::unordered_set>(compared_microproteins | std::views::filter(std::not_fn(is_single_copy_microprotein)));
        return detail::ClassifiedMicroproteins{ .single_copy = unique_microproteins, .homologous = homologous_microproteins, .homology_table = homology_table };
    }


    void keepTopHomologues(const mipfinder::homology::Results homology_search_results, const std::size_t maximum_homologues_allowed)
    {
        auto homology_table = detail::createHomologyTable(homology_search_results);
        detail::HomologyTable filtered_homology_table;
        for (const auto& [protein, homologues] : homology_table) {
            if (homologues.size() > maximum_homologues_allowed) {
                continue;
            }
            filtered_homology_table[protein] = homologues;
        }

        //mipfinder::homology::Results filtered_results;
        //for (const auto& homology_result : homology_search_results) {
        //    if (homology_result.query && homology_result.target)
        //}

        ////std::ranges::copy_if(homology_table, std::begin(filtered_table), homologue_count_filter);


    }


    //Find all potential microproteins from a proteome based on predetermined criteria
    //
    //@Params
    //proteome - A container of proteins to find microproteins from.
    //run_params - Parameters specifying microprotein characteristics.
    //hmmer_params - Parameters for running homology search.
    //homology_search_output - Path to a folder where the homology search result output
    //                         file will be created.
    //
    //@Return - A struct containing single-copy and homologous microproteins, as well as a homology relationship table
    //for the homologous microproteins.
    template <typename T>
    requires std::ranges::range<T>&& requires (typename std::ranges::range_value_t<T> t)
    {
        { t.length() } -> std::convertible_to<std::size_t>;
        { t.identifier() } -> std::convertible_to<std::string>;
    }
    ClassifiedMicroproteins<mipfinder::ProteinList, mipfinder::ProteinList>
        findMicroproteins(const T& proteome,
                          const mipfinder::Mipfinder::RunParameters& run_params,
                          const mipfinder::Mipfinder::HmmerParameters hmmer_params,
                          const std::filesystem::path& homology_search_output)
    {
        LOG(DEBUG) << "Finding proteins that meet the size criteria for microProteins";
        //Filter out all proteins in the proteome that are too long to be microproteins
        const std::size_t minimum_allowed_microprotein_length = run_params.minimum_microprotein_length;
        const std::size_t maximum_allowed_microprotein_length = run_params.maximum_microprotein_length;
        auto microprotein_filter = [&](const auto& protein) { return protein.length() <= run_params.maximum_microprotein_length; };
        auto potential_microproteins = detail::toContainer<std::vector>(proteome | std::views::filter(microprotein_filter));

        if (std::ranges::size(potential_microproteins) == 0) {
            throw std::runtime_error("No microProteins found in the proteome, aborting processing");
        }

        //Find homologous microproteins 
        detail::compareMicroproteinsToMicroproteins(potential_microproteins, hmmer_params, homology_search_output);
        //Filter out all microprotein homology results below bitscore_cutoff as these do not denote real
        //homologous relationships
        auto microprotein_homology_results = mipfinder::homology::parseResults(homology_search_output);
        const double lowest_allowed_microprotein_homology_bitscore = run_params.microprotein_homologue_bitscore_cutoff;
        auto strong_homologous_matches = mipfinder::homology::filterByBitscore(microprotein_homology_results, lowest_allowed_microprotein_homology_bitscore);


        //Filter out microproteins with more than @maximum_homologues_allowed homologues. This
        //ensures that we do not pick up large protein families that may contribute to false-positives
        //due to these microproteins containing domains that are very common, e.g. zinc fingers
        const auto maximum_allowed_homologues_per_microprotein = run_params.maximum_homologues_per_microprotein;
        detail::keepTopHomologues(strong_homologous_matches, maximum_allowed_homologues_per_microprotein);

        return detail::classifyMicroproteins(potential_microproteins, strong_homologous_matches);
    }

    //Find all potential ancestors from a proteome based on predetermined criteria.
    //
    //@Params
    //proteome - The proteome to search ancestors from.
    //run_params - Parameters extracted from the configuration file.
    //
    //@Return - 
    template <typename T>
    requires std::ranges::range<T>&& requires (std::ranges::range_value_t<T> v)
    {
        { v.length() } -> std::convertible_to<std::size_t>;
    }
    std::vector<std::ranges::range_value_t<T>>
        findAncestors(const T& proteome,
                      const mipfinder::Mipfinder::RunParameters& run_params)
    {
        const std::size_t minimum_allowed_ancestor_length = run_params.minimum_ancestor_length;
        const std::size_t maximum_allowed_ancestor_length = run_params.maximum_ancestor_length;
        auto ancestor_filter = [&](const auto& protein) { return protein.length() >= minimum_allowed_ancestor_length && protein.length() <= maximum_allowed_ancestor_length; };

        auto real_ancestors = proteome | std::views::filter(ancestor_filter);
        if (std::ranges::distance(real_ancestors) == 0) {
            throw std::runtime_error("Could not find any ancestors, aborting processing");
        }
        return detail::toContainer<std::vector>(real_ancestors);
    }

    //Create the results folder for the miPFinder run. The folder name consists of the current
    //date and time, followed by the organism identifier. The folder will be created in the
    //current working directory.
    //
    //@Params
    //organism_identifier - Name of the organism being analysed.
    //
    //@Return - Path to the created folder. 
    //
    //@Throws
    //std::runtime_error - If the results folder cannot be created.
    std::filesystem::path createMipfinderRunResultsFolder(const std::string organism_identifier)
    {
        LOG(DEBUG) << "Creating miPFinder run results folder";
        auto current_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        const std::string date_format = "%Y_%m_%d_%H_%M_%S"; //Get current date in the YYYY_MM_DD_Hour_Min_Sec format
        std::stringstream formatted_date;
        formatted_date << std::put_time(std::localtime(&current_time), date_format.c_str());
        std::filesystem::path results_folder{formatted_date.str() + "_results_" + organism_identifier};
        LOG(DEBUG) << "Creating " << std::filesystem::current_path() / results_folder;
        std::error_code e;
        if (!std::filesystem::create_directory(results_folder, e)) {
            throw(std::runtime_error("Failed to create the results folder : \"" + e.message() + "\""));
        }
        return results_folder;
    }

    template <typename T, typename U>
    requires std::ranges::range<T>&& std::ranges::range<U>
        void findAncestorsOfSingleCopyMicroproteins(const T& single_copy_microproteins,
                                                    const U& potential_ancestors,
                                                    const mipfinder::Mipfinder::HmmerParameters& parameters,
                                                    const std::filesystem::path& homology_search_output)
    {
        LOG(DEBUG) << "Findind ancestors of single-copy microProteins";
        LOG(DEBUG) << "Comparing " << single_copy_microproteins.size() << " cMIPs";

        const std::filesystem::path results_path = homology_search_output.parent_path();
        const std::filesystem::path query_file_location = results_path / "single_copy_microproteins.fasta";
        const std::filesystem::path database_location = results_path / "ancestors.fasta";
        mipfinder::proteinToFasta(single_copy_microproteins, query_file_location);
        mipfinder::proteinToFasta(single_copy_microproteins, database_location);

        //Compare potential microproteins to potential ancestors
        // if (std::ranges::size(single_copy_microproteins) == 0) {
        //     LOG(DEBUG) << "No single-copy microProteins found, aborting";
        //     return;
        // }

        // if (std::ranges::size(potential_ancestors) == 0) {
        //     LOG(DEBUG) << "No ancestors found, aborting";
        //     return;
        // }

        LOG(INFO) << "Finding ancestors of single-copy microProteins";
        const auto extra_param = "--mx " + parameters.scoring_matrix;
        const std::string extra_phmmer_parameters = "--popen " + std::to_string(parameters.gap_open_probability)
            + " --pextend " + std::to_string(parameters.gap_extension_probability)
            + extra_param;
        mipfinder::homology::phmmer(query_file_location, database_location, homology_search_output, extra_param);
    }

    template <typename Comp>
    mipfinder::homology::Results keepTopHits(const mipfinder::homology::Results& homology_search_results, std::size_t hits_to_keep, Comp comparator)
    {
        struct Homologue
        {
            std::string identifier;
            double bitscore;
        };

        //Link every protein up with its homologue and the score
        std::unordered_map<std::string, std::vector<Homologue>> homology_table;

        for (const auto& elem : homology_search_results)
        {
            if (homology_table.contains(elem.query)) {
                homology_table.at(elem.query).push_back(Homologue{elem.target, elem.bitscore});
            }
            else {
                homology_table[elem.query];
            }
        }

        for (auto& [protein, homologues] : homology_table)
        {
            std::sort(std::begin(homologues), std::end(homologues), comparator);
        }

        mipfinder::homology::Results filtered;
        for (const auto& [protein, homologues] : homology_table) {
            std::size_t count = 0;
            for (auto it = std::begin(homologues); it != std::end(homologues); ++it) {
                if (count != hits_to_keep) {
                    auto homologue = *it;
                    filtered.push_back(mipfinder::homology::Result{.query = protein, .target = homologue.identifier, .bitscore = homologue.bitscore});
                }
            }
        }

        return filtered;
    }

    template <typename T>
    mipfinder::homology::Results
	filterByLengthDifference(const mipfinder::homology::Results& homology_results,
							 const T& proteome,
							 const std::size_t min_length_difference)
    {
        mipfinder::homology::Results filtered_results;

        for (const auto& [protein, homologues] : homology_results)
        {
            for (const auto& homologue : homologues) {
                
            }
        }

	    //mipfinder::homology::Results filtered;

	    //std::unordered_map<std::string, mipfinder::Protein*> lookup_table;
	    //for (const auto& protein : proteins) {
		   // lookup_table[protein->identifier()] = protein;
	    //}

	    //for (const auto& result : results) {
		   // const auto query = lookup_table.at(result.query);
		   // const auto target = lookup_table.at(result.target);

		   // if (target->length() <= query->length()) {
			  //  continue;
		   // }
		   // assert(target->length() > query->length());

		   // const auto length_difference = target->length() - query->length();
		   // if (length_difference < min_length_difference) {
			  //  continue;
		   // }

		   // filtered.push_back(result);
	    //}
	    //return filtered;
    }


    template <typename T, typename U>
    requires std::ranges::range<T>&& std::ranges::range<U>
        auto filterAncestorHomologySearchResults(const T& proteome, const U& microprotein_ancestor_homology_results, const mipfinder::Mipfinder::RunParameters parameters)
    {
        //Remove homologous results with bitscores outside the acceptable range
        auto homologues = mipfinder::homology::parseResults(microprotein_ancestor_homology_results);
        const double minimum_allowed_homology_bitscore = parameters.ancestor_bitscore_cutoff;
        const double maximum_allowed_homology_bitscore = 120;
        auto high_confidence_ancestors = mipfinder::homology::filterByBitscore(homologues, minimum_allowed_homology_bitscore, maximum_allowed_homology_bitscore);


        //Keep the top x ancestors for each MIP to ensure we don't pick up EVERY
        //protein with a similar domain. This would be a problem for very common
        //domains such as kinases, zinc fingers etc.
        const auto maximum_homologous_ancestors_per_microprotein_to_keep = parameters.maximum_homologues_per_microprotein;
        auto comparator = [](const auto& x, const auto& y)
        {
            return x.bitscore > y.bitscore;
        };
        auto top_homologues_only = detail::keepTopHits(high_confidence_ancestors, maximum_homologous_ancestors_per_microprotein_to_keep, comparator);

        //Filter out ancestors that are within x a.a of the cMIP. If a protein is chosen as an ancestor
        //but is only slightly larger than the microprotein, then the ancestor is highly unlikely
        //to contain another domain and cannot function as an ancestor of a microProtein because it
        //does not have another effector domain.
        const auto minimum_length_difference = parameters.minimum_length_difference;
        U filtered_homology_results;
        for (const auto& homology_hit : top_homologues_only) {
            auto contains_protein = [&homology_hit](const auto& protein)
            {
                return (protein.identifier() == homology_hit.query && protein.identifier() == homology_hit.target);
            };



            //if (!std::ranges::find_if(proteome, contains_protein) != std::ranges::end(proteome)) {
            //    continue;
            //}

            //if (proteome)
            //filtered_homology_results
            

        }
        // 
        // 
        // 
        //const int min_length_diff = std::stoi(config["MIP"]["min_length_difference"]);
        //results = mipfinder::homology::filterByLengthDifference(results,
        //														proteome.data(),
        //														min_length_diff);


        return high_confidence_ancestors;
    }



    //Takes a file containing unaligned sequences and creates a Multiple Sequence Alignment (MSA) out of them in FASTA format
    void createMultipleSequenceAlignment(const std::filesystem::path& sequences_to_align, const std::filesystem::path& aligned_msa_output_file)
    {
        LOG(INFO) << "Aligning homologous microproteins";

        const std::string clustalo_command = "clustalo -i "
            + sequences_to_align.string()
            + " -o "
            + aligned_msa_output_file.string();

        int sys_call_result = std::system(clustalo_command.c_str());
        if (sys_call_result != 0) {
            throw std::runtime_error("Could not find clustalo, please ensure that it is installed");
        }
    }




    //Create HMMER profiles for each protein in "homologous_microproteins" that has a homologue.
    template <typename T, typename U>
    requires std::ranges::range<T>&& std::ranges::range<U>&& requires (std::ranges::range_value_t<T> t, std::ranges::range_value_t<U> u)
    {
        { t.identifier() } -> std::convertible_to<std::string>;
        std::convertible_to<std::ranges::range_value_t<U>, std::string>;
        std::convertible_to<typename std::ranges::range_value_t<U>::second_type::value_type, std::string>
            || std::convertible_to<typename std::ranges::range_value_t<U>::second_type, std::string>;
    }
    void createHmmprofiles(const T& homologous_microproteins, const U& homology_relationship_table, const std::filesystem::path& hmmprofile_output_folder)
    {
        if (!std::filesystem::exists(hmmprofile_output_folder))
        {
            std::filesystem::create_directories(hmmprofile_output_folder);
        }

        LOG(INFO) << "Creating HMMER profiles from homologous microProtein sequences";
        std::size_t profiles_built_counter = 0;
        //Group all unique homologues together into one sequence based on homology
        for (const auto& [protein_identifier, homologue_identifiers] : homology_relationship_table) {
            std::unordered_set<std::ranges::range_value_t<std::remove_const_t<T>>> grouped_homologous_proteins;
            if (auto found_protein = std::ranges::find_if(homologous_microproteins, [&protein_identifier](const auto& protein)
            {
                return protein.identifier() == protein_identifier;
            }); found_protein != std::ranges::end(homologous_microproteins)) {
                grouped_homologous_proteins.insert(*found_protein);
            }
            else {
                continue;
            };

            for (const auto& homologue_identifier : homologue_identifiers) {
                if (auto found_protein = std::ranges::find_if(homologous_microproteins, [&homologue_identifier](const auto& protein)
                {
                    return protein.identifier() == homologue_identifier;
                }); found_protein != std::ranges::end(homologous_microproteins)) {
                    grouped_homologous_proteins.insert(*found_protein);
                }
            }

            //Create files of unaligned sequences as input to an alignment program
            const std::filesystem::path unaligned_sequence_dir = hmmprofile_output_folder / "unaligned";
            std::filesystem::create_directory(unaligned_sequence_dir);
            const std::filesystem::path unaligned_grouped_sequences_fasta = unaligned_sequence_dir / (protein_identifier + ".fasta");
            mipfinder::proteinToFasta(grouped_homologous_proteins, unaligned_grouped_sequences_fasta);

            //Align the grouped sequences
            const std::filesystem::path multiple_sequence_alignment_dir = hmmprofile_output_folder / "msa";
            std::filesystem::create_directory(multiple_sequence_alignment_dir);
            const std::filesystem::path output_msa_file = multiple_sequence_alignment_dir / (protein_identifier + ".msa");
            detail::createMultipleSequenceAlignment(unaligned_grouped_sequences_fasta, output_msa_file);

            //Create a hmmprofile from the aligned sequences
            const auto hmmprofile_output_file = hmmprofile_output_folder / (protein_identifier + ".hmmprofile");

            const auto hmmprofile_name_command = "-n " + protein_identifier; //Names the MSA profile as the protein identifier
            mipfinder::homology::buildHmmerProfile(output_msa_file, hmmprofile_output_file, hmmprofile_name_command);
            ++profiles_built_counter;
        }
        LOG(INFO) << "Done creating HMMER profiles";
        LOG(INFO) << "Created " << profiles_built_counter << " rofiles";
    }

    template <typename T, typename U, typename V>
    requires std::ranges::range<T>&& std::ranges::range<U>
        void findAncestorsOfHomologousMicroproteins(const T& homologous_microproteins,
                                                    const U& potential_ancestors,
                                                    const V& homology_relationship_table,
                                                    const mipfinder::Mipfinder::HmmerParameters& parameters, //TODO: Remove, is unused
                                                    const std::filesystem::path& homology_search_output_file)
    {
        const std::filesystem::path parent_folder = homology_search_output_file.parent_path();
        const std::filesystem::path hmmer_output_folder = parent_folder / "hmmprofile";

        detail::createHmmprofiles(homologous_microproteins, homology_relationship_table, hmmer_output_folder);

        LOG(INFO) << "Merging HMMprofiles";
        const std::filesystem::path merged_profile_file = hmmer_output_folder / "mpf_merged.hmmprofile";
        detail::mergeHmmprofileFiles(hmmer_output_folder, merged_profile_file);

        //WIP: Find a better location for this
        const std::filesystem::path ancestor_fasta_file = hmmer_output_folder / "all_ancestors";
        mipfinder::proteinToFasta(potential_ancestors, ancestor_fasta_file);
        LOG(INFO) << "Performing hmmmsearch";
        mipfinder::homology::hmmsearch(merged_profile_file, ancestor_fasta_file, homology_search_output_file);
    }
}



mipfinder::Mipfinder::Mipfinder(const std::filesystem::path& configuration_file)
{
    configureLogger();

    //Set up configuration parameters from file
    LOG(DEBUG) << "Setting up miPFinder run configuration";
    //extractConfiguration(configuration_file);
    Configuration config{ configuration_file };

    LOG(DEBUG) << "Setting HMMER parameters";
    m_hmmer_parameters = HmmerParameters{

        .gap_open_probability = std::stod(config["HMMER"]["gap_open_probability"]),
        .gap_extension_probability = std::stod(config["HMMER"]["gap_extend_probability"]),
        .scoring_matrix = config["HMMER"]["matrix"]
    };
    LOG(DEBUG) << "Setting file parameters";
    m_file_parameters = FileParamaters{
        .input_proteome = config["TARGET"]["input_proteome"],
        .known_microprotein_list = config["MIP"]["known_microprotein_list"],
        .interpro_database = config["INTERPRO"]["interpro_database"],
        .gene_ontology_database = config["GO"]["go_database"],
        .uniprot_to_intepro_id_conversion_file = config["INTERPRO"]["uniprot_to_interpro_id_conversion"],
        .uniprot_to_go_id_conversion_file = config["GO"]["uniprot_to_go_id_conversion"],
    };
    LOG(DEBUG) << "Setting run parameters";
    m_run_parameters = RunParameters{
        .minimum_microprotein_length = std::stoul(config["MIP"]["minimum_microprotein_length"]),
        .maximum_microprotein_length = std::stoul(config["MIP"]["maximum_microprotein_length"]),
        .minimum_ancestor_length = std::stoul(config["MIP"]["minimum_ancestor_length"]),
        .maximum_ancestor_length = std::stoul(config["MIP"]["maximum_ancestor_length"]),
        .maximum_homologues_per_microprotein = std::stoul(config["MIP"]["maximum_homologues_per_microprotein"]),
        .minimum_length_difference = std::stoul(config["MIP"]["minimum_length_difference"]),
        .maximum_ancestor_count = std::stoul(config["MIP"]["maximum_ancestor_count"]),
        .maximum_protein_existence_level = std::stoul(config["TARGET"]["maximum_protein_existence_level"]),
        .microprotein_homologue_bitscore_cutoff = std::stod(config["HMMER"]["microprotein_homology_cutoff"]),
        .ancestor_bitscore_cutoff = std::stod(config["HMMER"]["ancestor_bitscore_cutoff"]),
        .output_format = config["REPORT"]["format"],
        .organism_identifier = config["TARGET"]["organism_identifier"],
    };

    //LOG(DEBUG) << "Checking file dependencies";
    //checkFileDependencies(m_file_parameters);
    m_results_folder = detail::createMipfinderRunResultsFolder(m_run_parameters.organism_identifier);

    ////WIP: DO this at the end
    ///* Provide a copy of the run parameters with the results */
    //std::filesystem::path config_file_copy = results_folder_ / configuration_file;
    //std::filesystem::copy(configuration_file, config_file_copy, std::filesystem::copy_options::overwrite_existing);
}

namespace mipfinder
{

    //////////////////////////////
    //                          //
    //                          //
    //      MIPFINDER V2.0      //
    //                          //
    //                          //
    //////////////////////////////

    void Mipfinder::run()
    {
        LOG(INFO) << "Starting mipfinder v2.0";
        const auto proteome = detail::loadProteome(m_file_parameters.input_proteome);

        //Remove proteins whose existence level implies their transcripts do not
        const std::size_t maximum_allowed_existence_level = m_run_parameters.maximum_protein_existence_level;
        auto real_proteins = detail::removeSpuriousProteins(proteome, maximum_allowed_existence_level);
        LOG(INFO) << "Removed " << std::ranges::distance(proteome) - std::ranges::distance(real_proteins) << " proteins";

        LOG(INFO) << "Searching for all microProteins in the proteome";
        const std::filesystem::path classified_microproteins = m_results_folder / "all_microproteins_vs_microproteins.txt";
        auto potential_microproteins = detail::findMicroproteins(proteome, m_run_parameters, m_hmmer_parameters, classified_microproteins);
        LOG(INFO) << "Found " << potential_microproteins.single_copy.size() << " single-copy microProteins";
        LOG(INFO) << "Found " << potential_microproteins.homologous.size() << " homologous microProteins";

        LOG(INFO) << "Finding ancestors";
        auto all_potential_ancestors = detail::findAncestors(proteome, m_run_parameters);
        LOG(INFO) << "Found " << all_potential_ancestors.size() << " potential ancestors";

        if (std::ranges::size(potential_microproteins.single_copy) != 0) {
            //Deal with single copy cMIPS
            //Take all single-copy cMIPS and compare them against large proteins to find their ancestors
            auto unique_microproteins_vs_ancestors = m_results_folder / "unique_vs_ancestor.txt";
            detail::findAncestorsOfSingleCopyMicroproteins(potential_microproteins.single_copy, all_potential_ancestors, m_hmmer_parameters, unique_microproteins_vs_ancestors);
            detail::filterAncestorHomologySearchResults(proteome, unique_microproteins_vs_ancestors, m_run_parameters);
        }

        if (std::ranges::size(potential_microproteins.homologous) != 0) {
            //Deal with homologous cMIPS
            auto homologous_microproteins_vs_ancestors = m_results_folder / "homologous_vs_ancestor.txt";
            detail::findAncestorsOfHomologousMicroproteins(potential_microproteins.homologous, all_potential_ancestors, potential_microproteins.homology_table, m_hmmer_parameters, homologous_microproteins_vs_ancestors);
            detail::filterAncestorHomologySearchResults(proteome, homologous_microproteins_vs_ancestors, m_run_parameters);
        }

        LOG(INFO) << "Main miPFinder processing steps are finished";
        LOG(INFO) << "Detecting whether extra information needs to be processed...";

        //---------------------------------
        //Optional processing steps. These only get executed if the required files have been provided in the configuration
        //file.
        //---------------------------------

        //WIP: Interpro processing steps

        //WIP: GO processing steps


        //---------------------------------
        //createReport()


        //Scoring part
        //---------------------------------
        //This needs to be redone into a separate class or potentially a separate program!
        ///* Assign KNOWN_MIP type to all known mips */
        //Proteome known_mips(config_["MIP"]["known_mips_fasta"]);
        //for (const auto& mip : known_mips.data()) {
        //	auto mip_in_proteome = proteome_.find(mip->identifier());
        //	if (mip_in_proteome) {
        //		mip_in_proteome->setType(mipfinder::Protein::Type::KNOWN_MIP);
        //	}
        //}


        //assignHmmerScores(proteins_to_score, &scoring_algorithm);
        //assignHmmerScores(proteins_to_score, &scoring_algorithm);

        ///* Filter out proteins with a 0 score */
        //mipfinder::ProteinSet non_zero;
        //for (const auto& protein : proteome_.data()) {
        //	if (protein->score() != 0) {
        //		non_zero.insert(protein);
        //	}
        //}

        ///* Filter out proteins that are not CMIPS */
        //mipfinder::ProteinSet only_cmip_type;
        //for (const auto& protein : non_zero) {
        //	if ((protein->type() == mipfinder::Protein::Type::SINGLE_COPY) ||
        //		(protein->type() == mipfinder::Protein::Type::HOMOLOGOUS) ||
        //		(protein->type() == mipfinder::Protein::Type::KNOWN_MIP)) {
        //		only_cmip_type.insert(protein);
        //	}
        //}

        ///* Sort proteome based on score */
        //auto comp = [](const mipfinder::Protein* lhs, const mipfinder::Protein* rhs)
        //{
        //	return lhs->score() > rhs->score();
        //};
        //std::vector<Protein*> sorted_proteome{only_cmip_type.begin(), only_cmip_type.end()};
        //std::sort(sorted_proteome.begin(), sorted_proteome.end(), comp);

        //std::filesystem::path results_fasta_file{"cmip_results.fasta"};
        //std::filesystem::path results_fasta_location = results_folder_ / results_fasta_file;
        //mipfinder::proteinToFasta(sorted_proteome, results_fasta_location);

        //LOG(INFO) << "Writing final report";
        //const std::string report_format = config_["REPORT"]["format"];
        //std::filesystem::path final_results = results_folder_ / "final_results.txt";
        //mipfinder::printer::createReport(report_format,
        //								 '\t',
        //								 sorted_proteome,
        //								 final_results);

        //LOG(INFO) << "mpf v2.0 has finished.";
    }
}

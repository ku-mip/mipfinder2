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

    /* For each cMIP in @results (HMMER queries), add the identified ancestor
     * (HMMER targets) to the cMIP. Does not add self as an ancestor */
     //void associateAncestorsWithCmips(const mipfinder::Proteome& proteome,
     //								 const mipfinder::homology::Results& results)
     //{
     //	for (const auto& result : results) {
     //		if (result.query == result.target) {
     //			continue;
     //		}

     //		const auto cmip = proteome.find(result.query);
     //		const auto ancestor = proteome.find(result.target);

     //		/* If ancestors do not have any domains or family annotations, ignore them */
     //		const auto all_ancestor_domains = ancestor->interproEntries();
     //		unsigned int domain_family_count{0};
     //		for (const auto& domain : all_ancestor_domains)
     //			if (domain.type == mipfinder::Interpro::Type::DOMAIN_TYPE ||
     //				domain.type == mipfinder::Interpro::Type::FAMILY) {
     //				++domain_family_count;
     //			}

     //		if (domain_family_count == 0) {
     //			continue;
     //		}

     //		assert(cmip != nullptr);
     //		assert(ancestor != nullptr);
     //		cmip->addAncestor(mipfinder::Ancestor{ancestor, result.bitscore});
     //	}
     //}

     ///* For each cMIP in @results (HMMER queries), add the identified homologues
     // * (HMMER targets, proteins) to the cMIP. Does not add self as a homologue */
     //void associateHomologuesWithCmips(const mipfinder::ProteinSet& cmips,
     //								  const mipfinder::homology::Results& results)
     //{
     //	/* Create a lookup table */
     //	std::unordered_map<std::string, mipfinder::Protein*> protein_lookup_table;
     //	for (const auto& cmip : cmips) {
     //		protein_lookup_table[cmip->identifier()] = cmip;
     //	}

     //	for (const auto& result : results) {
     //		if (result.query == result.target) {
     //			continue;
     //		}

     //		if (protein_lookup_table.count(result.query) == 0) {
     //			continue;
     //		}

     //		if (protein_lookup_table.count(result.target) == 0) {
     //			continue;
     //		}

     //		const auto cmip = protein_lookup_table[result.query];
     //		const auto homologue = protein_lookup_table[result.target];
     //		assert(cmip != nullptr);
     //		assert(homologue != nullptr);
     //		cmip->addHomologue(mipfinder::Result{homologue, result.bitscore});
     //	}
     //}

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

     ///* Divides a proteome into potential cMIPs and potential ancestor protein
     // * based on lengths specified in the configuration file. */
     //std::pair<mipfinder::ProteinSet, mipfinder::ProteinSet>
     //	divideProteome(const mipfinder::Proteome& proteome,
     //				   const mipfinder::Configuration& config)
     //{
     //	std::size_t max_mip_length = std::stoi(config["MIP"]["max_mip_length"]);
     //	std::size_t min_ancestor_length = std::stoi(config["MIP"]["min_ancestor_length"]);

     //	const auto cmips = mipfinder::filterByLength(proteome.data(),
     //												 1,
     //												 max_mip_length);
     //	const auto ancestors = mipfinder::filterByLength(proteome.data(),
     //													 min_ancestor_length);

     //	/* Mark respective parts of the proteome as a cMIP or an ancestor */
     //	for (const auto& cmip : cmips) {
     //		assert(cmip->type() == mipfinder::Protein::Type::UNKNOWN);
     //		cmip->setType(mipfinder::Protein::Type::CMIP);
     //	}
     //	for (const auto& ancestor : ancestors) {
     //		assert(ancestor->type() == mipfinder::Protein::Type::UNKNOWN);
     //		ancestor->setType(mipfinder::Protein::Type::ANCESTOR);
     //	}

     //	return std::make_pair(cmips, ancestors);
     //}

     /* Applies a set of filters to results file that represent cMIPS being
      * searched against potential ancestors. Returns a list of filtered HMMER
      * results */
      //mipfinder::homology::Results filterAncestors(mipfinder::homology::Results results,
      //										  const mipfinder::Proteome& proteome,
      //										  const mipfinder::Configuration& config)
      //{
      //	/* Apply filters to the unique cMIP vs ancestor results to eliminate
      //	* low-confidence results */
      //	const double ancestor_bitscore_cutoff = std::stod(config["HMMER"]["ancestor_bitscore_cutoff"]);
      //	results = mipfinder::homology::filterByBitscore(results,
      //												 ancestor_bitscore_cutoff);

      //	results = mipfinder::homology::filterByBitscore(results,
      //												 120,
      //												 std::less_equal<>());

      //	/* Keep the top x ancestors for each MIP to ensure we don't pick up EVERY
      //	* protein with a similar domain. This would be a problem for very common
      //	* domains such as kinases, zinc fingers etc. */
      //	const std::size_t hits_to_keep = std::stoi(config["MIP"]["max_ancestor_count"]);
      //	results = mipfinder::homology::keepTopHits(results, hits_to_keep);

      //	//Filter out ancestors that are within 40 a.a of the cMIP
      //	const int min_length_diff = std::stoi(config["MIP"]["min_length_difference"]);
      //	results = mipfinder::homology::filterByLengthDifference(results,
      //														 proteome.data(),
      //														 min_length_diff);

      //	results = filterSingleDomainAncestors(proteome, results);

      //	return results;
      //}

      //void setCmipTypes(const mipfinder::Proteome& proteome)
      //{
      //	for (const auto& protein : proteome.data()) {
      //		if (protein->type() != mipfinder::Protein::Type::CMIP) {
      //			continue;
      //		}
      //		if (protein->homologues().empty()) {
      //			protein->setType(mipfinder::Protein::Type::SINGLE_COPY);
      //		}
      //		else {
      //			protein->setType(mipfinder::Protein::Type::HOMOLOGOUS);
      //		}
      //	}
      //}

}


namespace detail
{
    template <typename T, typename U>
    struct ClassifiedMicroproteins
    {
        T single_copy;
        U homologous;
    };

    mipfinder::ProteinList loadProteome(const std::filesystem::path& fasta_file)
    {
        LOG(INFO) << "Loading proteome...";
        auto file = mipfinder::file::open(fasta_file);
        const mipfinder::FastaRecords proteome_fasta_records = mipfinder::fasta::extractRecords(file);

        mipfinder::ProteinList proteome;
        const char separator = '_';
        for (const auto& [header, sequence] : proteome_fasta_records) {

            const auto& [protein_id, sequence_version, description, existence_level] = mipfinder::fasta::extractUniprotHeader(header);
            const std::string identifier = protein_id + separator + sequence_version;

            proteome.emplace_back(mipfinder::Protein{ identifier, sequence, description, std::stoi(existence_level) });
        }
        LOG(INFO) << "Found " << proteome.size() << " proteins";
        return proteome;
    }

    //Compare microproteins to each other to establish which ones have homologues (homologous microproteins)
    //and which ones do not (single-copy microproteins).
    template <typename Cont>
    void compareMicroproteinsToMicroproteins(const Cont& potential_microproteins,
        const mipfinder::Mipfinder::HmmerParameters& parameters,
        const std::filesystem::path& results_output)
    {
        LOG(INFO) << "Finding homologous relationship between potential microproteins";
        const auto extra_param = "--mx " + parameters.scoring_matrix;
        mipfinder::homology::phmmer(potential_microproteins, potential_microproteins, results_output, extra_param);
        LOG(INFO) << "Finished";
    }

    //Filter out proteins whose existence level hints suggests that they are not translated transcripts
    template <typename T>
    requires std::ranges::range<T>&& requires (std::ranges::range_value_t<T> v)
    {
        { v.existenceLevel() } -> std::convertible_to<std::size_t>;
    }
    auto removeSpuriousProteins(T&& proteome, const std::size_t maximum_allowed_existence_level)
    {
        LOG(INFO) << "Cleaning up proteome";
        LOG(DEBUG) << "Removing proteins with existence level equal to or less than " << maximum_allowed_existence_level;
        auto protein_existence_filter = [=](const auto& protein)
        {
            return (protein.existenceLevel() <= maximum_allowed_existence_level);
        };
        return (proteome | std::views::filter(protein_existence_filter));
    }

    //WIP REDO THIS DOCUMENTATION
    //Return a table where each key is a protein query and each value is a container containing protein targets of homology search
    template <typename T>
    requires std::ranges::range<T>&& requires (std::ranges::range_value_t<T> v)
    {
        { v.query } -> std::convertible_to<std::string>;
        { v.target } -> std::convertible_to<std::string>;
    }
    auto createHomologyTable(T&& homology_search_results)
    {
        std::unordered_map<std::string, std::unordered_set<std::string>> homology_table; //Keys are protein identifiers, values are homologous protein identifiers

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
        return homology_table;
    }

    //Classify microProteins into single copy and homologous microproteins based on the homology search results
    template <typename T, typename U>
    requires std::ranges::range<T>&& std::ranges::range<U>
        && requires (std::ranges::range_value_t<T> t, std::ranges::range_value_t<U> u)
    {
        { t.identifier() } -> std::convertible_to<std::string>;
        { u.query } -> std::convertible_to<std::string>;
        { u.target } -> std::convertible_to<std::string>;
    }
    auto classifyMicroproteins(T&& potential_microproteins, U& homology_search_results)
    {
        //Find which potential microproteins have homologues (homologous microProteins) and
        //which do not (single-copy)
        auto homology_table = detail::createHomologyTable(homology_search_results);

        auto find_unique_microproteins = [&homology_table](const auto& elem)
        {
            //Check in case the homology results files are from a different protein set that the one
            //used for potential microproteins
            if (homology_table.contains(elem.identifier())) {
                if (homology_table.at(elem.identifier()).size() == 0) {
                    return true;
                }
            }
        };

        auto unique_potential_microproteins = potential_microproteins | std::views::filter(find_unique_microproteins);
        auto homologous_potential_microproteins = potential_microproteins | std::views::filter(find_unique_microproteins);
        return detail::ClassifiedMicroproteins{ .single_copy = unique_potential_microproteins, .homologous = homologous_potential_microproteins };
    }

    //Find all potential microproteins from a proteome based on predetermined criteria
    template <typename T>
    requires std::ranges::range<T>&& requires (typename std::ranges::range_value_t<T> t)
    {
        { t.length() } -> std::convertible_to<std::size_t>;
    }
    auto findMicroproteins(T&& proteome,
        const mipfinder::Mipfinder::RunParameters& run_params,
        const mipfinder::Mipfinder::HmmerParameters hmmer_params,
        const std::filesystem::path& homology_search_results)
    {
        //Filter out all proteins in the proteome that are too long to be microproteins
        const std::size_t minimum_allowed_microprotein_length = run_params.minimum_microprotein_length;
        const std::size_t maximum_allowed_microprotein_length = run_params.maximum_microprotein_length;
        auto microprotein_filter = [&](const auto& protein) { return protein.length() <= run_params.maximum_microprotein_length; };
        auto potential_microproteins = proteome | std::views::filter(microprotein_filter);

        //Find homologous microproteins 
        detail::compareMicroproteinsToMicroproteins(potential_microproteins, hmmer_params, homology_search_results);
        //Filter out all microprotein homology results below bitscore_cutoff as these do not denote real
        //homologous relationships
        auto microprotein_homology_results = mipfinder::homology::parseResults(homology_search_results);
        const double lowest_allowed_microprotein_homology_bitscore = run_params.microprotein_homologue_bitscore_cutoff;
        auto strong_homologous_matches = mipfinder::homology::filterByBitscore(microprotein_homology_results, lowest_allowed_microprotein_homology_bitscore);

        ////Find corresponding proteins in the proteome based on the homology search results
        //auto homology_filter = [&](const auto& protein)
        //{
        //    for (const auto& homology_result : strong_homologous_matches) {
        //        if (homology_result.query == protein.identifier()) { return true; }
        //    }
        //};
        //auto found_microproteins = proteome | std::views::filter(homology_filter);

        return detail::classifyMicroproteins(potential_microproteins, strong_homologous_matches);
    }

    //Find all potential ancestors from a proteome based on predetermined criteria
    template <typename T>
    requires std::ranges::range<T>&& requires (std::ranges::range_value_t<T> v)
    {
        { v.length() } -> std::convertible_to<std::size_t>;
    }
    auto findAncestors(T&& proteome,
        const std::filesystem::path& homology_search_results,
        const mipfinder::Mipfinder::RunParameters& run_params)
    {
        const std::size_t minimum_allowed_ancestor_length = run_params.minimum_ancestor_length;
        const std::size_t maximum_allowed_ancestor_length = run_params.maximum_ancestor_length;
        auto ancestor_filter = [&](const auto& protein) { return protein.length() >= minimum_allowed_ancestor_length && protein.length() <= maximum_allowed_ancestor_length; };

        auto potential_ancestors = proteome | std::views::filter(ancestor_filter);
        return potential_ancestors;
    }

    template <typename T>
    requires std::ranges::range<T>
        auto keepTopHits(T& container, std::size_t hits_to_keep)
    {

        //static_assert(std::same_as < std::ranges::range_value_t<T>, mipfinder::Protein>);
        //std::unordered_map<mipfinder::Result<mipfinder::Protein, mipfinder::Protein>, std::string> count_table;
        //std::unordered_map<mipfinder::Result<int, int>, std::string> m;
        //std::unordered_map<std::ranges::range_value_t<T>, std::size_t> count_table;
        T filtered;
        //for (const auto& elem : t) {
        //	count_table[elem] += 1;
        //	if (count_table[elem] <= hits_to_keep) {
        //		filtered.push_back(elem);
        //	}
        //}
        return filtered;
    }




    //void createFolders()
    //{
    //	//Creates the main results folder for the run
    //	const std::string organism_id = m_run_parameters.organism_identifier;
    //	const std::string folder_name = "results_" + organism_id;

    //	const std::filesystem::path results_folder{folder_name};
    //	std::filesystem::create_directory(results_folder);
    //	results_folder_ = results_folder;

    //	//Create the subfolders for individual software package results
    //	//createResultsFolder() has to be run before this is called, as it sets `results_folder_`
    //	msa_folder_ = results_folder_ / std::filesystem::path{"msa"};
    //	hmmprofile_folder_ = results_folder_ / std::filesystem::path{"hmmprofile"};
    //	homologue_folder_ = results_folder_ / std::filesystem::path{"homologues"};
    //	std::filesystem::create_directory(msa_folder_);
    //	std::filesystem::create_directory(hmmprofile_folder_);
    //	std::filesystem::create_directory(homologue_folder_);
    //}

    std::filesystem::path createMipfinderRunResultsFolder(const std::string organism_identifier)
    {
        LOG(DEBUG) << "Creating miPFinder run results folder";
        auto current_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        const std::string date_format = "%Y_%m_%d_%H_%M_%S"; //Get current date in the YYYY_MM_DD_Hour_Min_Sec format
        std::stringstream formatted_date;
        formatted_date << std::put_time(std::localtime(&current_time), date_format.c_str());

        //Creates the miPFinder run results folder
        std::filesystem::path results_folder{ formatted_date.str() + "_results_" + organism_identifier };
        LOG(DEBUG) << "Creating " << std::filesystem::current_path() / results_folder;
        std::error_code e;
        if (!std::filesystem::create_directory(results_folder, e)) {
            throw(std::runtime_error("Failed to create the results folder : \"" + e.message() + "\""));
        }
        return results_folder;
    }

    template <typename T, typename U>
    requires std::ranges::range<T>&& std::ranges::range<U>
        void findAncestorsOfSingleCopyMicroproteins(T& single_copy_microproteins,
            U& potential_ancestors,
            const mipfinder::Mipfinder::HmmerParameters& parameters,
            const std::filesystem::path& homology_search_output)
    {
        //Compare potential microproteins to potential ancestors
        const auto extra_param = "--mx " + parameters.scoring_matrix;
        const std::string extra_phmmer_parameters = "--popen " + std::to_string(parameters.gap_open_probability)
            + " --pextend " + std::to_string(parameters.gap_extension_probability)
            + extra_param;
        mipfinder::homology::phmmer(single_copy_microproteins, potential_ancestors, homology_search_output, extra_param);
        return;
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
        auto filtered_ancestor_homologues = detail::keepTopHits(high_confidence_ancestors, maximum_homologous_ancestors_per_microprotein_to_keep);

        ////Filter out ancestors that are within 40 a.a of the cMIP
        //const int min_length_diff = std::stoi(config["MIP"]["min_length_difference"]);
        //results = mipfinder::homology::filterByLengthDifference(results,
        //														proteome.data(),
        //														min_length_diff);


        return high_confidence_ancestors;
    }




    template <typename T, typename U>
    requires std::ranges::range<T>&& std::ranges::range<U>
        void findAncestorsOfHomologousMicroproteins(T& homologous_microproteins,
            U& potential_ancestors,
            const mipfinder::Mipfinder::HmmerParameters& parameters,
            const std::filesystem::path& homology_search_output)
    {


        ///* Build HMM profiles of all homologous cMIPS and hmmsearch these profiles
        // * against the ancestors */
        //writeHomologuesToFasta(homologous_cmips);
        //alignHomologousCmips(homologue_folder_);

        //if (config_["DEBUG"]["create_hmm_profiles"] == "true") {
        //	hmmer::createHmmProfiles(msa_folder_, hmmprofile_folder_);
        //}

        //std::filesystem::path filename{"all_hmmprofiles.txt"};
        //std::filesystem::path unified_hmmprofile_file = results_folder_ / filename;
        //const auto all_hmmer_profiles = hmmer::concatenateHmmProfiles(hmmprofile_folder_,
        //															  unified_hmmprofile_file);

        ///* Merges all individual hmm profile files into one large one */
        //const auto homologous_vs_ancestor_search =
        //	hmmsearchHomologousCmips(all_hmmer_profiles);
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
      .ancestor_bitscore_cutoff = std::stod(config["HMMER"]["ancestor_bitscore_cutoff"]),
      .microprotein_homologue_bitscore_cutoff = std::stod(config["HMMER"]["microprotein_homology_cutoff"]),
      .output_format = config["REPORT"]["format"],
      .organism_identifier = config["TARGET"]["organism_identifier"],
    };

    //LOG(DEBUG) << "Checking file dependencies";
    //checkFileDependencies(m_file_parameters);

    //LOG(DEBUG) << "Setting up GO database";
    //mipfinder::Go go_database{config_["GO"]["go_database"]};
    //LOG(DEBUG) << "Adding GO ids to proteins";
    //addGoIdentifiers(proteome_,
    //				 go_database,
    //				 config_["GO"]["uniprot_to_go"]);

    //LOG(DEBUG) << "Setting up InterPro database";
    //mipfinder::Interpro interpro_database{config_["INTERPRO"]["interpro_database"]};
    //LOG(DEBUG) << "Adding InterPro ids to proteins";
    //addInterproIdentifiers(proteome_,
    //					   interpro_database,
    //					   config_["INTERPRO"]["uniprot_to_interpro"]);


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
        auto proteome = detail::loadProteome(m_file_parameters.input_proteome);

        //Remove proteins whose existence level implies their transcripts do not
        const std::size_t maximum_allowed_existence_level = m_run_parameters.maximum_protein_existence_level;
        auto real_proteins = detail::removeSpuriousProteins(proteome, maximum_allowed_existence_level);
        LOG(INFO) << "Removed " << std::ranges::distance(proteome) - std::ranges::distance(real_proteins) << " proteins";

        const std::filesystem::path microprotein_homology_search_results = m_results_folder / "all_microproteins_vs_microproteins.txt";
        auto potential_microproteins = detail::findMicroproteins(proteome, m_run_parameters, m_hmmer_parameters, microprotein_homology_search_results);



        //const std::filesystem::path classified_ancestors = m_results_folder / "all_microproteins_vs_ancestors.txt";
        //auto all_potential_ancestors = detail::findAncestors(proteome, classified_ancestors, m_run_parameters);

        ////Deal with single copy cMIPS
        ////Take all single-copy cMIPS and compare them against large proteins to find their ancestors
        //auto unique_microproteins_vs_ancestors = m_results_folder / "unique_vs_ancestor.txt";
        //detail::findAncestorsOfSingleCopyMicroproteins(unique_potential_microproteins, all_potential_ancestors, m_hmmer_parameters, unique_microproteins_vs_ancestors);
        //detail::filterAncestorHomologySearchResults(proteome, unique_microproteins_vs_ancestors, m_run_parameters);

        ////Deal with homologous cMIPS
        //auto homologous_microproteins_vs_ancestors = m_results_folder / "homologous_vs_ancestor.txt";
        //detail::findAncestorsOfHomologousMicroproteins(homologous_unique_microproteins, all_potential_ancestors, m_hmmer_parameters, homologous_microproteins_vs_ancestors);
        //detail::filterAncestorHomologySearchResults(proteome, homologous_microproteins_vs_ancestors);




        //const auto unique_vs_ancestor_file =
        //	phmmerCmipsVsAncestors(unique_cmips, "unique_vs_proteome.txt");

        //auto unique_vs_ancestor_results = hmmer::parseResults(unique_vs_ancestor_file);

        //unique_vs_ancestor_results = filterAncestors(unique_vs_ancestor_results,
        //											 proteome_,
        //											 config_);

        //associateAncestorsWithCmips(proteome_,
        //							unique_vs_ancestor_results);



        //auto homologous_vs_ancestor_results = hmmer::parseResults(homologous_vs_ancestor_search);

        //homologous_vs_ancestor_results = filterAncestors(homologous_vs_ancestor_results,
        //												 proteome_,
        //												 config_);

        //associateAncestorsWithCmips(proteome_,
        //							homologous_vs_ancestor_results);

        //setCmipTypes(proteome_);

        ///* Assign KNOWN_MIP type to all known mips */
        //Proteome known_mips(config_["MIP"]["known_mips_fasta"]);
        //for (const auto& mip : known_mips.data()) {
        //	auto mip_in_proteome = proteome_.find(mip->identifier());
        //	if (mip_in_proteome) {
        //		mip_in_proteome->setType(mipfinder::Protein::Type::KNOWN_MIP);
        //	}
        //}

        ///* Filter out proteins with more than @maximum_homologues_allowed homologues */
        //mipfinder::ProteinSet proteins_to_score;
        //for (const auto& protein : proteome_.data()) {
        //	const double max_homologues_allowed = std::stod(config_["MIP"]["maximum_homologues"]);
        //	if (protein->homologues().size() > max_homologues_allowed) {
        //		continue;
        //	}
        //	proteins_to_score.insert(protein);
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

    //void Mipfinder::createFolders()
    //{
    //	//Creates the main results folder for the run
    //	const std::string organism_id = m_run_parameters.organism_identifier;
    //	const std::string folder_name = "results_" + organism_id;

    //	const std::filesystem::path results_folder{folder_name};
    //	std::filesystem::create_directory(results_folder);
    //	results_folder_ = results_folder;

    //	//Create the subfolders for individual software package results
    //	//createResultsFolder() has to be run before this is called, as it sets `results_folder_`
    //	msa_folder_ = results_folder_ / std::filesystem::path{"msa"};
    //	hmmprofile_folder_ = results_folder_ / std::filesystem::path{"hmmprofile"};
    //	homologue_folder_ = results_folder_ / std::filesystem::path{"homologues"};
    //	std::filesystem::create_directory(msa_folder_);
    //	std::filesystem::create_directory(hmmprofile_folder_);
    //	std::filesystem::create_directory(homologue_folder_);
    //}

    //std::filesystem::path Mipfinder::phmmerAgainstSelf(const mipfinder::ProteinSet& cmips)
    //{
    //	std::filesystem::path results_file{"all_cmips_against_cmips.txt"};
    //	const auto hmmer_cmips_vs_cmips = results_folder_ / results_file;

    //	if (config_["DEBUG"]["phmmer_cmips"] == "false") {
    //		return std::filesystem::path{hmmer_cmips_vs_cmips};
    //	}

    //	LOG(INFO) << "Grouping cMIPs based on homology";

    //	const auto scoring_matrix = config_["HMMER"]["matrix"];
    //	const auto extra_param = " --mx " + scoring_matrix;

    //	mipfinder::homology::phmmer(cmips, cmips, hmmer_cmips_vs_cmips.string(), extra_param);
    //	LOG(INFO) << "Finished grouping cMIPS";
    //	return hmmer_cmips_vs_cmips;
    //}

    //void Mipfinder::writeHomologuesToFasta(const mipfinder::ProteinSet& homologous_cmips)
    //{
    //	//if (config_["DEBUG"]["write_homologous_groups_fasta"] == "false") { return; }

    //	LOG(INFO) << "Creating FASTA files of homologous cMIPs";
    //	int grouped_fasta_files{0};
    //	for (const auto& protein : homologous_cmips) {
    //		mipfinder::ProteinSet grouped_homologues;
    //		grouped_homologues.insert(protein);

    //		for (const auto& homologue : protein->homologues()) {
    //			grouped_homologues.insert(homologue.protein);
    //		}

    //		const auto id = protein->identifier();

    //		const std::filesystem::path unaligned_filename{id + "_homologues.fasta"};
    //		const std::filesystem::path unaligned_file_location{homologue_folder_
    //			/ unaligned_filename};

    //		proteinToFasta(grouped_homologues, unaligned_file_location.string());
    //		++grouped_fasta_files;
    //	}
    //	LOG(DEBUG) << "Created " << grouped_fasta_files << " homologous cMIP sequence "
    //		<< "FASTA files";
    //}

    //void Mipfinder::alignHomologousCmips(const std::filesystem::path& unaligned_seq_folder)
    //{
    //	//if (config_["DEBUG"]["align_homologous_cmips"] == "false") { return; }

    //	LOG(INFO) << "Aligning groups of homologous cMIPs";

    //	typedef std::filesystem::directory_iterator DirectoryIter;

    //	int msas_created{0};
    //	for (const auto& file : DirectoryIter{unaligned_seq_folder}) {

    //		//Unaligned FASTA files are named as `PROTEIN_ID`_homologues.fasta. We
    //		//need to extract the `PROTEIN_ID` part for the aligned filename.
    //		const auto unaligned_filename = file.path().filename().string();
    //		const auto tokens = mipfinder::tokenise(unaligned_filename, '_');
    //		const auto protein_id = tokens[0];

    //		const std::filesystem::path unaligned_file_path = file.path();

    //		const std::filesystem::path msa_filename{protein_id + "_aligned.fasta"};
    //		const std::filesystem::path msa_file_full_path =
    //			msa_folder_ / msa_filename;

    //		const std::string clustalo_command = "clustalo -i "
    //			+ unaligned_file_path.string()
    //			+ " -o "
    //			+ msa_file_full_path.string();

    //		int sys_call_result = std::system(clustalo_command.c_str());
    //		if (sys_call_result != 0) {
    //			continue;
    //		}
    //		++msas_created;
    //	}
    //	LOG(INFO) << "Created " << msas_created << " MSA from grouped cMIP homologues";
    //}

    //std::filesystem::path
    //	Mipfinder::hmmsearchHomologousCmips(const std::filesystem::path& hmmprofiles_file)
    //{
    //	unsigned int min_ancestor_length =
    //		std::stoi(config_["MIP"]["min_ancestor_length"]);

    //	const auto ancestors =
    //		mipfinder::filterByLength(proteome_.data(), min_ancestor_length);

    //	const std::filesystem::path results_filename{"homologous_vs_proteome.txt"};
    //	const std::filesystem::path results_file_location =
    //		results_folder_ / results_filename;

    //	if (config_["DEBUG"]["find_homologous_cmip_ancestors"] == "true") {
    //		mipfinder::homology::hmmsearch(hmmprofiles_file, ancestors, results_file_location);
    //	}
    //	return results_file_location;
    //}

    //void Mipfinder::assignHmmerScores(const mipfinder::ProteinSet& proteins,
    //								  HmmerScoringAlgorithm algorithm)
    //{
    //	int score_counter{0};
    //	for (const auto& protein : proteins) {
    //		/* Only cMIPs have ancestors */
    //		if (protein->ancestors().empty()) {
    //			continue;
    //		}

    //		for (const auto& ancestor : protein->ancestors()) {
    //			const double bitscore = ancestor.bitscore;

    //			double score = algorithm(HmmerScoringData{protein, ancestor.protein, bitscore});
    //			protein->changeScore(score);
    //			++score_counter;
    //		}
    //	}
    //}

}

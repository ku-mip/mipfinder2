#include <algorithm>
#include <concepts>
#include <cmath>
#include <filesystem>
#include <functional>
#include <numeric>
#include <set>
#include <unordered_set>
#include <ranges>
#include <regex>

#include "configuration.h"
#include "easylogging++.h"
#include "fasta.h"
#include "hmmer.h"
#include "helpers.h"
#include "interpro.h"
#include "mipfinder.h"
#include "protein.h"




/* Helper functions */
namespace
{
    /* Ensures that all files required by mipfinder are found. If any
     * dependencies are missing, throws a std::runtime_error.  */
    void checkFileDependencies(mipfinder::Mipfinder::FileParamaters file_parameters)
    {
        if (!std::filesystem::exists(file_parameters.input_proteome)
            || !std::filesystem::exists(file_parameters.known_microprotein_list)
            || !std::filesystem::exists(file_parameters.gene_ontology_database)
            || !std::filesystem::exists(file_parameters.interpro_database)
            || !std::filesystem::exists(file_parameters.uniprot_to_intepro_id_conversion_file)
            || !std::filesystem::exists(file_parameters.uniprot_to_go_id_conversion_file)) {
            throw std::runtime_error("Could not find the required dependencies");
        }
    }
}


namespace detail
{
    struct UniProtHeader
    {
        unsigned int sequence_version;
        unsigned int existence_level;
        std::string accession;
        std::string description;
    };

    /* Takes a UniProt database FASTA header as input and returns an array with the
     * following elements:
     * 1: UniProt accession name
     * 2: Entry sequence version
     * 3: Entry description
     * 4: Entry protein existence level
     *
     * If the header cannot be successfully parsed, returns an empty array.
     */
    UniProtHeader extractUniprotHeader(const std::string& header)
    {
        /* This regex currently only grabs the accession name, description and
         * sequence version because the other data is not interesting */
        static const std::regex uniprot_regex(R"((?:>)(\w+)(?:\|)(\w+)(?:\|)(\w+)[ ](.+?)[ ](?:OS=).+(?:PE=)([0-9]).+(?:SV=)([0-9]))");

        std::string accession_name;
        std::string sequence_version;
        std::string description;
        std::string existence_level;

        std::smatch matches;
        if (std::regex_search(header, matches, uniprot_regex)) {

            if (matches.size() != 7) {
                return UniProtHeader{};
            }

            /* For a standard UniProt header the matches will as following:
             * Match 0 - whole match
             * Match 1 - database type
             * Match 2 - UniProt accession name
             * Match 3 - Entry name
             * Match 4 - Description
             * Match 5 - Protein existence level
             * Match 6 - Sequence version */
            accession_name = matches[2];
            description = matches[4];
            existence_level = matches[5];
            sequence_version = matches[6];

            return UniProtHeader{ .sequence_version = std::stoul(sequence_version),
                     .existence_level = std::stoul(existence_level),
                     .accession = accession_name,
                     .description = description };
        }
        return UniProtHeader{};
    }


    using ProteinDomains = std::unordered_map<mipfinder::protein::Identifier, std::unordered_set<std::string>>;

    /**
     * @brief  Associate each UniProt identifier with their InterPro domains
     * @param  uniprot_entry_interpro_domains  A tsv-file that specifies which 
     *                                         UniProt entry has which InterPro
     *                                         domains.
     *                                         Column 1: UniProt accession 
     *                                         without sequence version.
     *                                         Column 2: Sequence version of the
     *                                         corresponding UniProt accession.
     *                                         Column 3: Commma (;) separated list
     *                                         of InterPro entry identifiers.
     * @return  Associative array where each key is a protein identifier and
     *          the values are a list of unique InterPro identifiers.
     *
     * If the input file is not in a correct format, the behaviour is unspecified.
     */
    ProteinDomains parseProteinDomainList(const std::filesystem::path& uniprot_entry_interpro_domains)
    {
        std::ifstream file{ uniprot_entry_interpro_domains };
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open " +
                uniprot_entry_interpro_domains.string() +
                ", aborting...");
        }

        ProteinDomains parsed_list{};
        std::string line;
        while (std::getline(file, line)) {
            auto tokens = mipfinder::tokenise(line, '\t');
            std::string uniprot_accession = tokens[0];
            unsigned int sequence_version = std::stoul(tokens[1]);
            mipfinder::protein::Identifier id{ uniprot_accession, sequence_version };
            auto interpro_identifiers = mipfinder::tokenise(tokens[2], ';');

            for (const auto& interpro_id : interpro_identifiers) {
                parsed_list[id].insert(interpro_id);
            }
        }
        return parsed_list;
    }


    //Take all hmmprofile files (ending in ".hmmprofile") in "hmmprofile_directory" and create one file with all the profiles as "output_file"
    void mergeHmmprofileFiles(const std::filesystem::path& hmmprofile_directory, const std::filesystem::path& output_file)
    {
        typedef std::filesystem::directory_iterator DirectoryIter;

        std::ofstream of;
        of.open(output_file, std::ios::trunc);
        for (const auto& directory_entry : DirectoryIter{ hmmprofile_directory }) {
            //If the output file is in the same directory as the hmmprofile files, 
            //do not also process the output file! Otherwise it'll be an infinite loop.
            if (directory_entry.path() == output_file) {
                continue;
            }

            LOG(DEBUG) << "Processing " << directory_entry.path();
            if (directory_entry.path().extension() != std::filesystem::path{ ".hmmprofile" }) {
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

    /**
     * @brief  Extract protein data from a FASTA file.
     * @param  fasta_file  Location of the FASTA file.
     *
     * @throws  std::runtime_error  If @a fasta_file cannot be opened.
     * @return  A collection of Proteins sorted by their identifier.
     *
     * The FASTA headers in the given @a fasta_file need to be in UniProt
     * format. If they are not, the behaviour is unspecified.
     */
    mipfinder::protein::ProteinList loadProteins(const std::filesystem::path& fasta_file)
    {
        LOG(DEBUG) << "Entering loadProteins() with argument: " << fasta_file;

        const mipfinder::fasta::Entries fasta_entries = mipfinder::fasta::parse(fasta_file);
        mipfinder::protein::ProteinList proteins{};
        for (const auto& entry : fasta_entries) {
            const auto& header = detail::extractUniprotHeader(entry.header);
            const mipfinder::protein::Identifier identifier{
                header.accession,
                header.sequence_version };

            proteins.emplace_back(mipfinder::protein::Protein{
                identifier,
                entry.sequence,
                header.description,
                header.existence_level });
        }
        LOG(DEBUG) << "Found " << proteins.size() << " proteins";

        LOG(DEBUG) << "Sorting proteins based on their identifier";
        std::sort(std::begin(proteins), std::end(proteins),
            [](const mipfinder::protein::Protein& first, const mipfinder::protein::Protein& second)
            { return first.identifier() < second.identifier(); });
        LOG(DEBUG) << "Finished sorting";

        LOG(DEBUG) << "Exiting loadProteins()";
        return proteins;
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
        mipfinder::protein::proteinToFasta(potential_microproteins, query_file_location);

        mipfinder::homology::phmmer(query_file_location, query_file_location, results_output, extra_param);
        LOG(DEBUG) << "Finished finding microProtein homologues";
    }

    /**
     * @brief  Filter out proteins whose existence level is lower than allowed.
     * @param  proteins  List of proteins to filter.
     * @param  max_allowed_existence_level  Cutoff value for existence level as
                                            defined by UniProt.
     * @return  A list of proteins which meet the criteria for real proteins.
     */
    template <typename T>
    requires std::ranges::range<T>&& requires (std::ranges::range_value_t<T> t)
    {
        { t.existenceLevel() } -> std::convertible_to<std::size_t>;
    }
    mipfinder::protein::ProteinList filterByExistenceLevel(const T& proteins, const std::size_t max_allowed_existence_level)
    {
        LOG(DEBUG) << "Entering filterByExistenceLevel()";
        LOG(DEBUG) << "Removing proteins with existence level equal to or less than " << max_allowed_existence_level;

        T filtered_proteins{};
        for (const auto& protein : proteins) {
            if (protein.existenceLevel() <= max_allowed_existence_level) {
                filtered_proteins.push_back(protein);
            }
        }

        LOG(DEBUG) << "Leaving filterByExistenceLevel()";
        return filtered_proteins;
    }

    struct ClassifiedMicroproteins
    {
        mipfinder::protein::ProteinList single_copy;
        mipfinder::protein::ProteinList homologous;
        mipfinder::homology::Results homology_table;
    };

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
    ClassifiedMicroproteins classifyMicroproteins(const mipfinder::protein::ProteinList& microproteins,
        const mipfinder::homology::Results& homology_search_results)
    {
        LOG(DEBUG) << "Classifying microProteins into single-copy and homologous";
        std::unordered_map<mipfinder::protein::Identifier, std::size_t> count_table;
        //for (const auto& result : homology_search_results) {
        //    ++count_table[mipfinder::protein::Identifier{result.query}];
        //}

        mipfinder::protein::ProteinList single_copy_microproteins;
        mipfinder::protein::ProteinList homologous_microproteins;
        for (const auto& protein : microproteins) {
            if (!count_table.contains(protein.identifier())) {
                continue;
            }

            if (count_table.at(protein.identifier()) == 1) {
                single_copy_microproteins.push_back(protein);
            }
            else {
                homologous_microproteins.push_back(protein);
            }
        }
        return detail::ClassifiedMicroproteins{ .single_copy = single_copy_microproteins, .homologous = homologous_microproteins, .homology_table = homology_search_results };
    }

    /**
     * @brief  Remove proteins whose domain count is outside allowed range.
     * @return  A list of proteins 
     * 
     * Only considers InterPro entries that are annotated as 'Domain' type,
     * and ignores the rest. I.e. if a protein has 5 InterPro annotations but
     * none of them are domains (but rather 'active site', 'family' etc), it
     * will get removed.
     */
    template <typename T, typename U, typename V>
    //TODO: Sensible constraints
    requires std::ranges::range<T> &&
        requires (typename std::ranges::range_value_t<T> t)
    {
        {t.identifier() } -> std::same_as < mipfinder::protein::Identifier>;
    }
    mipfinder::protein::ProteinList filterByDomainCount(const T& microproteins,
        std::size_t minimum_allowed_domains,
        std::size_t maximum_allowed_domains,
        const V& interpro_database,
        const U& uniprot_to_interpro_table)
    {
        mipfinder::protein::ProteinList filtered;
        for (const auto& microprotein : microproteins) {
            if (!uniprot_to_interpro_table.contains(microprotein.identifier())) {
                continue;
            }

            if (uniprot_to_interpro_table.at(microprotein.identifier()).size() < minimum_allowed_domains
                || uniprot_to_interpro_table.at(microprotein.identifier()).size() > minimum_allowed_domains) {
                continue;
            }

            for (const auto& interpro_entry_accession : uniprot_to_interpro_table.at(microprotein.identifier())) {
                if (auto interpro_entry = mipfinder::find(interpro_database, interpro_entry_accession); interpro_entry != interpro_database.cend()) {
                    if (interpro_entry->type != mipfinder::Interpro::Entry::Type::domain_type) {
                        continue;
                    }
                    filtered.push_back(microprotein);
                }
            }
        }
        return filtered;
    }

    /**
     * @brief  Finds all proteins that satisfy the length criteria for a
     *         microprotein.
     * @param  proteins  A collection of proteins to filter.
     * @param  minimum_protein_length  Minimum length of a protein to keep,
     *                                 inclusive.
     * @param  maximum_protein_length  Maximum length of a protein to keep,
     *                                 inclusive.
     * @return  A collection of proteins of [@a minimum_protein_length,
     *          @a maximum_protein_length].
     */
    template <typename T>
    requires std::ranges::range<T>&&
        requires (typename std::ranges::range_value_t<T> t)
    {
        { t.sequence().length() } -> std::convertible_to<std::size_t>;
    }
    T filterProteinsByLength(const T& proteins,
        std::size_t minimum_protein_length,
        std::size_t maximum_protein_length)
    {
        LOG(DEBUG) << "Entering filterProteinsByLength()";
        T potential_microproteins{};
        for (const auto& elem : proteins) {
            if (elem.sequence().length() >= minimum_protein_length &&
                elem.sequence().length() <= maximum_protein_length) {
                potential_microproteins.push_back(elem);
            }
        }
        return potential_microproteins;
        LOG(DEBUG) << "Removed " << proteins.size() - potential_microproteins.size()
                   << " proteins that did not match length criteria.";
        LOG(DEBUG) << "Exiting filterProteinsByLength()";
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
        { t.sequence().length() } -> std::convertible_to<std::size_t>;
        { t.identifier().to_string() } -> std::convertible_to<std::string>;
    }
    ClassifiedMicroproteins findMicroproteins(
        const T& proteins,
        const mipfinder::Mipfinder::RunParameters& run_params,
        const mipfinder::Mipfinder::HmmerParameters& hmmer_params,
        const mipfinder::Mipfinder::FileParamaters& file_params,
        const std::filesystem::path& homology_search_output)
    {
        LOG(DEBUG) << "Entering filterProteinsByLength()";

        LOG(DEBUG) << "Finding microProteins that meet the length criteria as defined in configuration";
        const std::size_t minimum_allowed_microprotein_length = run_params.minimum_microprotein_length;
        const std::size_t maximum_allowed_microprotein_length = run_params.maximum_microprotein_length;
        auto potential_microproteins = detail::filterProteinsByLength(
            proteins,
            minimum_allowed_microprotein_length,
            maximum_allowed_microprotein_length);

        if (std::ranges::size(potential_microproteins) == 0) {
            throw std::runtime_error("After initial filtering, no potential microProteins were found in the proteome. Stopping miPFinder");
        }

        //Find homologous microproteins 
        detail::compareMicroproteinsToMicroproteins(potential_microproteins, hmmer_params, homology_search_output);
        //Filter out all microprotein homology results below bitscore_cutoff as these do not denote real
        //homologous relationships
        auto microprotein_homology_results = mipfinder::homology::parseHmmerTabularResults(homology_search_output);
        const double lowest_allowed_microprotein_homology_bitscore = run_params.microprotein_homologue_bitscore_cutoff;
        auto strong_homologous_matches = mipfinder::homology::filterByBitscore(microprotein_homology_results, lowest_allowed_microprotein_homology_bitscore);

        //Filter out microproteins with more than @maximum_homologues_allowed homologues. This
        //ensures that we do not pick up large protein families that may contribute to false-positives
        //due to these microproteins containing domains that are very common, e.g. zinc fingers
        const auto maximum_allowed_homologues_per_microprotein = run_params.maximum_homologues_per_microprotein;
        auto no_large_protein_families = mipfinder::homology::keepTopHomologues(strong_homologous_matches, maximum_allowed_homologues_per_microprotein);

        LOG(DEBUG) << "Exiting filterProteinsByLength()";
        return detail::classifyMicroproteins(potential_microproteins, no_large_protein_families);
    }

    /**
     * @brief  Find all proteins that meet the criteria for a microprotein
     *         ancestor.
     * @param  proteins  A collection of proteins to find ancestors from.
     * @return  A list of proteins that qualify as ancestors.
     */
    template <typename T>
    requires std::ranges::range<T>&& requires (std::ranges::range_value_t<T> v)
    {
        { v.sequence().length() } -> std::convertible_to<std::size_t>;
    }
    std::vector<std::ranges::range_value_t<T>>
        findAncestors(const T& proteins,
            const mipfinder::Mipfinder::RunParameters& run_params)
    {
        const std::size_t minimum_allowed_ancestor_length = run_params.minimum_ancestor_length;
        const std::size_t maximum_allowed_ancestor_length = run_params.maximum_ancestor_length;

        auto ancestors = filterProteinsByLength(proteins, minimum_allowed_ancestor_length, maximum_allowed_ancestor_length);

        if (std::ranges::distance(ancestors) == 0) {
            throw std::runtime_error("Could not find any ancestors, aborting processing");
        }
        return ancestors;
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
    /**
     * @brief  Creates a folder for mipfinder results.
     * @param  organism_identifier  Organism on which mipfinder is being run.
     * @throws  std::runtime_error  If results folder could not be created.
     * 
     * The function creates a new folder in the current working directory
     * prefixed with current date and the string "_results_" followed by
     * @a organism_identifier.
     */
    std::filesystem::path createMipfinderRunResultsFolder(const std::string organism_identifier)
    {
        LOG(DEBUG) << "Creating miPFinder run results folder";
        auto current_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        const std::string date_format = "%Y_%m_%d_%H_%M_%S"; //Get current date in the YYYY_MM_DD_Hour_Min_Sec format
        std::stringstream formatted_date;
        formatted_date << std::put_time(std::localtime(&current_time), date_format.c_str());
        std::filesystem::path results_folder{ formatted_date.str() + "_results_" + organism_identifier };
        LOG(DEBUG) << "Creating " << std::filesystem::current_path() / results_folder;
        std::error_code e;
        if (!std::filesystem::create_directory(results_folder, e)) {
            throw(std::runtime_error("Failed to create the results folder : \"" + e.message() + "\""));
        }
        return results_folder;
    }

    //Find which proteins are potential ancestors of the given single-copy microproteins, i.e.
    //which microproteins may have evolved from the larger proteins.
    //
    //@Params
    //single_copy_microproteins - A container 
    template <typename T, typename U>
    requires std::ranges::range<T>&& std::ranges::range<U>
        void findAncestorsOfSingleCopyMicroproteins(const T& single_copy_microproteins,
            const U& potential_ancestors,
            const mipfinder::Mipfinder::HmmerParameters& parameters,
            const std::filesystem::path& homology_search_output)
    {
        LOG(INFO) << "Finding ancestors of single-copy microProteins";

        if (std::ranges::size(single_copy_microproteins) == 0) {
            LOG(INFO) << "No single-copy microProteins found, aborting homology search";
            return;
        }

        if (std::ranges::size(potential_ancestors) == 0) {
            LOG(INFO) << "No ancestors found, aborting homology search";
            return;
        }

        const std::filesystem::path results_path = homology_search_output.parent_path();
        const std::filesystem::path query_file_location = results_path / "single_copy_microproteins.fasta";
        const std::filesystem::path database_location = results_path / "ancestors.fasta";
        mipfinder::protein::proteinToFasta(single_copy_microproteins, query_file_location);
        mipfinder::protein::proteinToFasta(single_copy_microproteins, database_location);

        const auto extra_param = "--mx " + parameters.scoring_matrix;
        const std::string extra_phmmer_parameters = "--popen " + std::to_string(parameters.gap_open_probability)
            + " --pextend " + std::to_string(parameters.gap_extension_probability)
            + extra_param;
        mipfinder::homology::phmmer(query_file_location, database_location, homology_search_output, extra_param);
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
    void filterOutTooShortAncestors(const T& proteome, const U& homology_search_results, const std::size_t minimum_length_difference)
    {
        U filtered_homology_results{};
        for (const auto& homology_hit : homology_search_results) {
            auto contains_protein = [&homology_hit](const auto& protein)
            {
                return (protein.identifier() == homology_hit.query && protein.identifier() == homology_hit.target);
            };

            if (!std::ranges::find_if(proteome, contains_protein) != std::ranges::end(proteome)) {
                continue;
            }
        }
    }



    template <typename T, typename U>
    requires std::ranges::range<T>&& std::ranges::range<U>
        auto filterAncestorHomologySearchResults(const T& proteome, const U& microprotein_ancestor_homology_results, const mipfinder::Mipfinder::RunParameters parameters)
    {
        //Remove homologous results with bitscores outside the acceptable range
        auto homologues = mipfinder::homology::parseHmmerTabularResults(microprotein_ancestor_homology_results);
        const double minimum_allowed_homology_bitscore = parameters.ancestor_bitscore_cutoff;
        const double maximum_allowed_homology_bitscore = 120;
        auto high_confidence_ancestors = mipfinder::homology::filterByBitscore(homologues, minimum_allowed_homology_bitscore, maximum_allowed_homology_bitscore);


        //Keep the top x ancestors for each MIP to ensure we don't pick up EVERY
        //protein with a similar domain. This would be a problem for very common
        //domains such as kinases, zinc fingers etc.
        const auto maximum_homologous_ancestors_per_microprotein_to_keep = parameters.maximum_homologues_per_microprotein;
        auto top_homologues_only = mipfinder::homology::keepTopHomologues(high_confidence_ancestors, maximum_homologous_ancestors_per_microprotein_to_keep);

        //Filter out ancestors that are within x a.a of the cMIP. If a protein is chosen as an ancestor
        //but is only slightly larger than the microprotein, then the ancestor is highly unlikely
        //to contain another domain and cannot function as an ancestor of a microProtein because it
        //does not have another effector domain.
        const auto minimum_length_difference = parameters.minimum_length_difference;






        //if (proteome)
        //filtered_homology_results

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
    void createHmmprofiles(const T& homologous_microproteins, const U& homology_relationship_table, const std::filesystem::path& hmmprofile_output_folder)
    {
        if (!std::filesystem::exists(hmmprofile_output_folder))
        {
            std::filesystem::create_directories(hmmprofile_output_folder);
        }

        LOG(INFO) << "Creating HMMER profiles from homologous microProtein sequences";
        std::size_t profiles_built_counter = 0;
        //Group all unique homologues together into one sequence based on homology
        //for (const auto& [protein_identifier, homologue_identifiers] : homology_relationship_table) {
        //    std::unordered_set<std::ranges::range_value_t<std::remove_const_t<T>>> grouped_homologous_proteins;
        //    if (auto found_protein = std::ranges::find_if(homologous_microproteins, [&protein_identifier](const auto& protein)
        //    {
        //        return protein.identifier() == protein_identifier;
        //    }); found_protein != std::ranges::end(homologous_microproteins)) {
        //        grouped_homologous_proteins.insert(*found_protein);
        //    }
        //    else {
        //        continue;
        //    };

        //    for (const auto& homologue_identifier : homologue_identifiers) {
        //        if (auto found_protein = std::ranges::find_if(homologous_microproteins, [&homologue_identifier](const auto& protein)
        //        {
        //            return protein.identifier() == homologue_identifier.target;
        //        }); found_protein != std::ranges::end(homologous_microproteins)) {
        //            grouped_homologous_proteins.insert(*found_protein);
        //        }
        //    }

        //    //Create files of unaligned sequences as input to an alignment program
        //    const std::filesystem::path unaligned_sequence_dir = hmmprofile_output_folder / "unaligned";
        //    std::filesystem::create_directory(unaligned_sequence_dir);
        //    const std::filesystem::path unaligned_grouped_sequences_fasta = unaligned_sequence_dir / (protein_identifier + ".fasta");
        //    mipfinder::proteinToFasta(grouped_homologous_proteins, unaligned_grouped_sequences_fasta);

        //    //Align the grouped sequences
        //    const std::filesystem::path multiple_sequence_alignment_dir = hmmprofile_output_folder / "msa";
        //    std::filesystem::create_directory(multiple_sequence_alignment_dir);
        //    const std::filesystem::path output_msa_file = multiple_sequence_alignment_dir / (protein_identifier + ".msa");
        //    detail::createMultipleSequenceAlignment(unaligned_grouped_sequences_fasta, output_msa_file);

        //    //Create a hmmprofile from the aligned sequences
        //    const auto hmmprofile_output_file = hmmprofile_output_folder / (protein_identifier + ".hmmprofile");

        //    const auto hmmprofile_name_command = "-n " + protein_identifier; //Names the MSA profile as the protein identifier
        //    mipfinder::homology::buildHmmerProfile(output_msa_file, hmmprofile_output_file, hmmprofile_name_command);
        //    ++profiles_built_counter;
        //}
        LOG(INFO) << "Done creating HMMER profiles";
        LOG(INFO) << "Created " << profiles_built_counter << " profiles";
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
        mipfinder::protein::proteinToFasta(potential_ancestors, ancestor_fasta_file);
        LOG(INFO) << "Performing hmmmsearch";
        mipfinder::homology::hmmsearch(merged_profile_file, ancestor_fasta_file, homology_search_output_file);
    }
}



mipfinder::Mipfinder::Mipfinder(const std::filesystem::path& configuration_file)
{
    //Set up configuration parameters from file
    LOG(DEBUG) << "Setting up miPFinder run configuration";
    //extractConfiguration(configuration_file);
    Configuration config{ configuration_file };

    LOG(DEBUG) << "Setting HMMER parameters";
    m_hmmer_parameters = HmmerParameters{
        .gap_open_probability = std::stod(config.value("HMMER", "gap_open_probability")),
        .gap_extension_probability = std::stod(config.value("HMMER", "gap_extend_probability")),
        .scoring_matrix = config.value("HMMER", "matrix")
    };

    LOG(DEBUG) << "Setting file parameters";
    m_file_parameters = FileParamaters{
        .input_proteome = config.value("TARGET", "input_proteome"),
        .known_microprotein_list = config.value("MIP", "known_microprotein_list"),
        .interpro_database = config.value("INTERPRO", "interpro_database"),
        .gene_ontology_database = config.value("GO", "go_database"),
        .uniprot_to_intepro_id_conversion_file = config.value("INTERPRO", "uniprot_to_interpro_id_conversion"),
        .uniprot_to_go_id_conversion_file = config.value("GO", "uniprot_to_go_id_conversion")
    };
    LOG(DEBUG) << "Setting run parameters";
    m_run_parameters = RunParameters{
        .minimum_microprotein_length = std::stoul(config.value("MIP", "minimum_microprotein_length")),
        .maximum_microprotein_length = std::stoul(config.value("MIP", "maximum_microprotein_length")),
        .minimum_ancestor_length = std::stoul(config.value("MIP", "minimum_ancestor_length")),
        .maximum_ancestor_length = std::stoul(config.value("MIP", "maximum_ancestor_length")),
        .maximum_homologues_per_microprotein = std::stoul(config.value("MIP", "maximum_homologues_per_microprotein")),
        .minimum_length_difference = std::stoul(config.value("MIP", "minimum_length_difference")),
        .maximum_ancestor_count = std::stoul(config.value("MIP", "maximum_ancestor_count")),
        .maximum_protein_existence_level = std::stoul(config.value("TARGET", "maximum_protein_existence_level")),
        .microprotein_homologue_bitscore_cutoff = std::stod(config.value("HMMER", "microprotein_homology_cutoff")),
        .ancestor_bitscore_cutoff = std::stod(config.value("HMMER", "ancestor_bitscore_cutoff")),
        .output_format = config.value("REPORT", "format"),
        .organism_identifier = config.value("TARGET", "organism_identifier")
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
        const auto proteome = detail::loadProteins(m_file_parameters.input_proteome);

        LOG(INFO) << "Performing intial filtering steps...";
        const std::size_t maximum_allowed_existence_level = m_run_parameters.maximum_protein_existence_level;
        auto potential_microproteins = detail::filterByExistenceLevel(proteome, maximum_allowed_existence_level);

        /**
         * Optional processing steps. These only get executed if the required
         * files have been provided in the configuration file.
         */
        if (!m_file_parameters.uniprot_to_intepro_id_conversion_file.empty() && !m_file_parameters.interpro_database.empty()) {
            LOG(INFO) << "Incorporating InterPro data into mipfinder";
            
            auto interpro_database = mipfinder::Interpro(m_file_parameters.interpro_database);

            /*
             * According to the current model, microProteins by definition only have
             * one domain. Therefore to be considered a microProtein, we have to
             * filter out every potential microProtein that has been annotated with
             * more than one InterPro entry of domain type.
             */
            auto uniprot_to_interpro_conversion_table = detail::parseProteinDomainList(m_file_parameters.uniprot_to_intepro_id_conversion_file);
            //TODO (09/09/21): Put the domain count parameters into the configuration file.
            constexpr std::size_t minimum_domains_per_microprotein = 1;
            constexpr std::size_t maximum_domains_per_microprotein = 1;
            potential_microproteins = detail::filterByDomainCount(potential_microproteins,
                minimum_domains_per_microprotein,
                maximum_domains_per_microprotein,
                interpro_database,
                uniprot_to_interpro_conversion_table);
        }
        /**
         * End of optional processing steps
         */

        LOG(INFO) << "Searching for all microProteins in the proteome";
        const std::filesystem::path classified_microproteins = m_results_folder / "all_microproteins_vs_microproteins.txt";
        auto categorised_microproteins = detail::findMicroproteins(potential_microproteins, m_run_parameters, m_hmmer_parameters, m_file_parameters, classified_microproteins);
        LOG(INFO) << "Found " << categorised_microproteins.single_copy.size() << " single-copy microProteins";
        LOG(INFO) << "Found " << categorised_microproteins.homologous.size() << " homologous microProteins";

        LOG(INFO) << "Finding ancestors";
        auto all_potential_ancestors = detail::findAncestors(proteome, m_run_parameters);
        LOG(INFO) << "Found " << all_potential_ancestors.size() << " potential ancestors";

        if (std::ranges::size(categorised_microproteins.single_copy) != 0) {
            //Deal with single copy cMIPS
            //Take all single-copy cMIPS and compare them against large proteins to find their ancestors
            auto unique_microproteins_vs_ancestors = m_results_folder / "unique_vs_ancestor.txt";
            detail::findAncestorsOfSingleCopyMicroproteins(categorised_microproteins.single_copy, all_potential_ancestors, m_hmmer_parameters, unique_microproteins_vs_ancestors);
            detail::filterAncestorHomologySearchResults(proteome, unique_microproteins_vs_ancestors, m_run_parameters);
        }

        if (std::ranges::size(categorised_microproteins.homologous) != 0) {
            //Deal with homologous cMIPS
            auto homologous_microproteins_vs_ancestors = m_results_folder / "homologous_vs_ancestor.txt";
            detail::findAncestorsOfHomologousMicroproteins(categorised_microproteins.homologous, all_potential_ancestors, categorised_microproteins.homology_table, m_hmmer_parameters, homologous_microproteins_vs_ancestors);
            detail::filterAncestorHomologySearchResults(proteome, homologous_microproteins_vs_ancestors, m_run_parameters);
        }

        LOG(INFO) << "Finished";

        //---------------------------------
        //Optional processing steps. These only get executed if the required files have been provided in the configuration
        //file.
        //---------------------------------
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

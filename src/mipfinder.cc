#include <algorithm>
#include <concepts>
#include <cmath>
#include <filesystem>
#include <functional>
#include <numeric>
#include <set>
#include <unordered_set>
#include <unordered_map>
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
    /**
     *  @brief  Checks that all required parameters have been set.
     *  @throws  std::runtime_error if any of the required parameters does not have a value.
     */
    void checkRequiredParameters(const mipfinder::Configuration& config)
    {
        //TODO: This needs to be updated with the new configuration settings
        static const std::unordered_map<std::string, std::vector<std::string>> required_parameters
        {
            { "TARGET", {
                "organism_identifier",
                "organism_proteome_fasta"
                }
            },
            { "MIP", {
                "max_mip_length",
                "min_ancestor_length",
                "maximum_homologues",
                "min_length_difference",
                "max_ancestor_count",
                }
            },
            { "HMMER", {
                "homologue_bitscore_cutoff",
                "ancestor_bitscore_cutoff",
                "gap_open_probability",
                "gap_extend_probability",
                "matrix"
                }
            },
            { "REPORT", {
                "format"
                }
            }
        };

        std::vector<std::string> missing_parameters;
        for (const auto& [header, parameters] : required_parameters) {
            for (const auto& parameter_name : parameters) {
                if (!config.contains(header, parameter_name)) {
                    missing_parameters.push_back(parameter_name);
                }
            }
        }

        if (missing_parameters.size() == 0) {
            std::string missing_parameter_info = join(' ', missing_parameters);
            throw std::runtime_error("The following required parameters are not set in the configuration file: \"" + missing_parameter_info + "\".");
        }
    }

    /**
     * @brief  Ensures that all files required by mipfinder are found.
     * @throw  std::out_of_range if any dependencies are missing.
     */
    void checkRequiredFiles(const mipfinder::Configuration& config)
    {
        //Required files
        if (!std::filesystem::is_regular_file(config.value("TARGET", "organism_proteins_fasta"))) {
            throw std::runtime_error("File specified by \"organism_proteins_fasta\" cannot be found.");
        }
    }

    /**
     * @brief  Verifies that all settings needed for mipfinder have been set.
     * @throws  std::runtime_error if the configuration is invalid.
     */
    void verifyConfiguration(const mipfinder::Configuration& config)
    {
        checkRequiredParameters(config);
        checkRequiredFiles(config);
    }
}


namespace detail
{
    //Find the corresponding proteins from the homology search results
    template <typename Cont>
    requires std::ranges::range<Cont>
        Cont findCorrespondingProteins(const mipfinder::homology::Results& results,
            const Cont& proteome)
    {
        /* Lookup table for fast searching */
        std::unordered_map<std::string, mipfinder::protein::Protein> lookup_table;
        for (const auto& protein : proteome) {
            lookup_table.insert(std::make_pair(protein.identifier(), protein));
        }

        Cont found_proteins{};
        for (const auto& result : results) {
            if (lookup_table.contains(result.query)) {
                found_proteins.push_back(lookup_table.at(result.query));
            }
        }

        //Remove duplicates
        std::sort(std::begin(found_proteins), std::end(found_proteins));
        auto new_last_element = std::unique(std::begin(found_proteins), std::end(found_proteins));
        found_proteins.erase(new_last_element, found_proteins.end());
        return found_proteins;
    }


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
    template <typename T, typename U>
    void compareMicroproteinsToMicroproteins(const T& potential_microproteins,
        const U& homology_search_parameters,
        const std::filesystem::path& results_output)
    {
        LOG(DEBUG) << "Finding homologous relationship between potential microproteins";
        const std::filesystem::path results_path = results_output.parent_path();
        const std::filesystem::path query_file_location = results_path / "microprotein_query.fasta";
        mipfinder::protein::proteinToFasta(potential_microproteins, query_file_location);

        mipfinder::homology::phmmer(query_file_location, query_file_location, results_output, homology_search_parameters);
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
    mipfinder::protein::ProteinList filterByExistenceLevel(T proteins, const unsigned int max_allowed_existence_level)
    {
        LOG(DEBUG) << "Removing proteins with existence level equal to or less than " << max_allowed_existence_level;
        auto isAboveAllowedExistenceLevel = [](const auto& elem) {
            return (elem.existenceLevel() > max_allowed_existence_level);
        };
        std::erase_if(proteins, isAboveAllowedExistenceLevel);
        return proteins;
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
                if (auto interpro_entry = mipfinder::interpro::find(interpro_database, interpro_entry_accession); interpro_entry != interpro_database.cend()) {
                    if (interpro_entry->type != mipfinder::interpro::Database::Entry::Type::domain_type) {
                        continue;
                    }
                    filtered.push_back(microprotein);
                }
            }
        }
        return filtered;
    }

    ///**
    // * @brief  Finds all proteins that satisfy the length criteria for a
    // *         microprotein.
    // * @param  proteins  A collection of proteins to filter.
    // * @param  minimum_protein_length  Minimum length of a protein to keep,
    // *                                 inclusive.
    // * @param  maximum_protein_length  Maximum length of a protein to keep,
    // *                                 inclusive.
    // * @return  A collection of proteins of [@a minimum_protein_length,
    // *          @a maximum_protein_length].
    // */
    //template <typename T>
    //    requires std::ranges::range<T>
    // && requires (typename std::ranges::range_value_t<T> t)
    //{
    //    { t.sequence().length() } -> std::convertible_to<std::size_t>;
    //}
    //T filterProteinsByLength(T proteins,
    //                         std::size_t minimum_protein_length,
    //                         std::size_t maximum_protein_length)
    //{
    //    auto isTooShortOrLong = [](const auto& elem) {
    //        return (elem.sequence().length() < minimum_protein_length ||
    //            elem.sequence().length() > maximum_protein_length);
    //    };
    //    std::erase_if(proteins, isTooShortOrLong);
    //    return proteins;
    //}

    /**
     * @brief  Find all proteins that meet the criteria for a microprotein
     *         ancestor.
     * @param  proteins  A collection of proteins to find ancestors from.
     * @return  A list of proteins that qualify as ancestors.
     */
    template <typename T>
    requires std::ranges::range<T> && requires (std::ranges::range_value_t<T> v)
    {
        { v.sequence().length() } -> std::convertible_to<std::size_t>;
    }
    std::vector<std::ranges::range_value_t<T>> findAncestors(const T& proteins,
                                                             const std::size_t min_ancestor_length,
                                                             const std::size_t max_ancestor_length)
    {
        auto ancestors = filterProteinsByLength(proteins, min_ancestor_length, max_ancestor_length);
        if (std::ranges::distance(ancestors) == 0) {
            throw std::runtime_error("Could not find any potential microprotein ancestors in the proteome, stopping mipfinder.");
        }
        return ancestors;
    }


    /**
     * @brief  Creates a folder for mipfinder results.
     * @param  organism_identifier  Organism on which mipfinder is being run.
     * @throws  std::runtime_error  If results folder could not be created.
     * 
     * The function creates a new folder in the current working directory
     * prefixed with current date and the string "_results_" followed by
     * @a organism_identifier.
     */
    std::filesystem::path createResultsFolder(const std::string& organism_identifier)
    {
        LOG(DEBUG) << "Creating mipfinder run results folder";
        auto current_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        constexpr auto date_format = "%Y_%m_%d_%H_%M_%S"; //Get current date in the YYYY_MM_DD_Hour_Min_Sec format
        std::stringstream formatted_date;
        formatted_date << std::put_time(std::localtime(&current_time), date_format);

        std::filesystem::path results_folder{ formatted_date.str() + "_results_" + organism_identifier };
        LOG(DEBUG) << "Creating " << std::filesystem::current_path() / results_folder;
        std::error_code e;
        if (!std::filesystem::create_directory(results_folder, e)) {
            throw(std::runtime_error("Failed to create the results folder : \"" + e.message() + "\""));
        }
        return results_folder;
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
            throw std::runtime_error("Could not find \"clustalo\", please ensure that it is installed and available in the PATH variable");
        }
    }

    /**
     * @brief Saves a copy of the configuration file used for the mipfinder in the results folder.
     */
    void saveRunConfiguration(const std::filesystem::path& results_folder, const std::filesystem::path& configuration_file)
    {
        std::filesystem::path config_file_copy = results_folder / configuration_file;
        std::filesystem::copy(configuration_file, config_file_copy, std::filesystem::copy_options::overwrite_existing);
    }

    //std::vector<std::string> assembleHmmerParametersFromConfiguration(const mipfinder::Configuration& configuration)
    //{
    //    auto homology_matrix = "--mx " + configuration.value("HMMER", "matrix");
    //    auto gap_open_probability = "--popen " + configuration.value("HMMER", "gap_open_probability");
    //    auto gap_extend_probability = "--pextend " + configuration.value("HMMER", "gap_extend_probability");
    //    return std::vector<std::string>{homology_matrix, gap_extend_probability, gap_open_probability};
    //}
}

namespace
{
    mipfinder::Mipfinder::HomologyParameters setHomologyParameters(const mipfinder::Configuration& configuration)
    {
        return mipfinder::Mipfinder::HomologyParameters{
            .homologue_bitscore_cutoff = configuration.parseAs<double>("HOMOLOGY", "homologue_bitscore_cutoff"),
            .minimum_ancestor_bitscore = configuration.parseAs<double>("HOMOLOGY", "minimum_ancestor_bitscore"),
            .maximum_ancestor_bitscore = configuration.parseAs<double>("HOMOLOGY", "maximum_ancestor_bitscore"),
            .gap_open_probability = configuration.parseAs<double>("HOMOLOGY", "gap_open_probability"),
            .gap_extend_probability = configuration.parseAs<double>("HOMOLOGY", "gap_extend_probability"),
            .matrix = configuration.parseAs<std::string>("HOMOLOGY", "matrix"),
        };
    }

    mipfinder::Mipfinder::FileParameters setFileParameters(const mipfinder::Configuration& configuration)
    {
        return mipfinder::Mipfinder::FileParameters{
            .proteome = configuration.parseAs<std::string>("ORGANISM", "proteome_file"),
            .interpro_database = configuration.contains("INTERPRO", "interpro_database")
                ? std::optional<std::filesystem::path>{configuration.parseAs<std::string>("INTERPRO", "interpro_database")} : std::nullopt,
            .uniprot_to_intepro = configuration.contains("INTERPRO", "uniprot_to_interpro")
                ? std::optional<std::filesystem::path>{configuration.parseAs<std::string>("INTERPRO", "uniprot_to_interpro")} : std::nullopt,
            .go_database = configuration.contains("GO", "go_database")
                ? std::optional<std::filesystem::path>{configuration.parseAs<std::string>("GO", "go_database")} : std::nullopt,
            .uniprot_to_go = configuration.contains("GO", "uniprot_to_go")
                ? std::optional<std::filesystem::path>{configuration.parseAs<std::string>("GO", "uniprot_to_go")} : std::nullopt,
            .known_microprotein_identifiers = configuration.contains("MIP", "known_microprotein_identifiers") 
                ? std::optional<std::filesystem::path>{configuration.parseAs<std::string>("MIP", "known_microprotein_identifiers")} : std::nullopt
        };
    }

    mipfinder::Mipfinder::GeneralParameters setGeneralParameters(const mipfinder::Configuration& configuration)
    {
        return mipfinder::Mipfinder::GeneralParameters{
            .organism_identifier = configuration.parseAs<std::string>("ORGANISM", "identifier")
        };
    }

    mipfinder::Mipfinder::MicroproteinParameters setMicroproteinParameters(const mipfinder::Configuration& configuration)
    {
        return mipfinder::Mipfinder::MicroproteinParameters{
            .minimum_ancestor_length = configuration.parseAs<unsigned long>("MICROPROTEIN", "minimum_ancestor_length"),
            .maximum_ancestor_length = configuration.parseAs<unsigned long>("MICROPROTEIN", "maximum_ancestor_length"),
            .minimum_length_difference = configuration.parseAs<unsigned long>("MICROPROTEIN", "minimum_length_difference"),
            .maximum_ancestor_homologues = configuration.parseAs<unsigned long>("MICROPROTEIN", "maximum_ancestor_homologues"),
            .minimum_domains_per_microprotein = configuration.parseAs<unsigned long>("MICROPROTEIN", "minimum_domains_per_microprotein"),
            .maximum_domains_per_microprotein = configuration.parseAs<unsigned long>("MICROPROTEIN", "maximum_domains_per_microprotein"),
            .minimum_microprotein_length = configuration.parseAs<unsigned long>("MICROPROTEIN", "minimum_microprotein_length"),
            .maximum_microprotein_length = configuration.parseAs<unsigned long>("MICROPROTEIN", "maximum_microprotein_length"),
            .maximum_allowed_protein_existence = configuration.contains("MICROPROTEIN", "maximum_allowed_protein_existence")
                ? std::optional<uint16_t>{configuration.parseAs<unsigned long>("MICROPROTEIN", "maximum_allowed_protein_existence")} : std::nullopt,
        };
    }

    /**
     * @brief Saves a copy of the configuration file used for the mipfinder in the results folder.
     */
    void saveRunConfiguration(const std::filesystem::path& results_folder, const std::filesystem::path& configuration_file)
    {
        std::filesystem::path config_file_copy = results_folder / configuration_file;
        std::filesystem::copy(configuration_file, config_file_copy, std::filesystem::copy_options::overwrite_existing);
    }
}






namespace mipfinder
{
    Mipfinder::Mipfinder(const std::filesystem::path& configuration_file) : configuration_file(configuration_file)
    {
        LOG(DEBUG) << "Verifying mipfinder run configuration";
        setConfiguration(configuration_file);
    }

    void Mipfinder::setConfiguration(const std::filesystem::path& configuration_file)
    {
        Configuration configuration{ configuration_file };
        verifyConfiguration(configuration);

        homology_parameters = setHomologyParameters(configuration);
        file_parameters = setFileParameters(configuration);
        general_parameters = setGeneralParameters(configuration);
        microprotein_parameters = setMicroproteinParameters(configuration);
    }


    void Mipfinder::run()
    {
        results_folder = detail::createResultsFolder(general_parameters.organism_identifier);
        detail::saveRunConfiguration(results_folder, configuration_file);

        LOG(INFO) << "Starting mipfinder v2.0";
        const auto proteome = detail::loadProteins(file_parameters.proteome);

        auto candidate_microproteins = findCandidateMicroproteins(proteome);
        if (std::ranges::empty(candidate_microproteins)) {
            LOG(ERROR) << "No candidate microproteins found in the proteome, exiting.";
            return;
        }

        auto candidate_ancestors = findCandidateAncestors(proteome);
        if (std::ranges::empty(candidate_ancestors)) {
            LOG(ERROR) << "No candidate ancestors found in the proteome, exiting.";
            return;
        }

        LOG(INFO) << "Searching for all microProteins in the proteome";
        auto categorised_microproteins = classifyCandidateMicroproteins(candidate_microproteins);
        LOG(INFO) << "Found " << categorised_microproteins.single_copy.size() << " single-copy microProteins";
        LOG(INFO) << "Found " << categorised_microproteins.homologous.size() << " homologous microProteins";


        mipfinder::homology::Results single_copy_vs_ancestor_homology_results;
        if (std::ranges::empty(categorised_microproteins.single_copy)) {
            LOG(ERROR) << "No single copy microproteins detected, skipping analysis.";
        }
        else {
            single_copy_vs_ancestor_homology_results = findAncestorsOfSingleCopyMicroproteins(categorised_microproteins.single_copy, candidate_ancestors);
        }

        mipfinder::homology::Results homologous_vs_ancestor_homology_results;
        if (std::ranges::empty(categorised_microproteins.homologous)) {
            LOG(ERROR) << "No homologous microproteins detected, skipping analysis.";
        }
        else {
            findAncestorsOfHomologousMicroproteins(categorised_microproteins.homologous,
                candidate_ancestors,
                categorised_microproteins.homology_table);
        }

        //TODO: Filter ancestor results by length difference, i.e. an ancestor who is less than X amino acids longer than
        //the microprotein is unlikely to be a real ancestor due to not containing another domain like ancestors must do.
        //filterByLength();
        //filterByLength();

        //filterByDomainCount();

        LOG(INFO) << "Finished";

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

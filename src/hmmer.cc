#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <unordered_set>
#include <ranges>

#include "configuration.h"
#include "easylogging++.h"

#include "hmmer.h"
#include "interpro.h"
#include "protein.h"

namespace mipfinder::homology
{
    bool operator==(const Result& lhs, const Result& rhs)
    {
        /* Since bitscores only ever have one digit, comparing them for equality
         * should be safe */
        return lhs.query == rhs.query &&
            lhs.target == rhs.target &&
            lhs.bitscore == rhs.bitscore;
    }


    void phmmer(const std::filesystem::path& query_file,
        const std::filesystem::path& database_file,
        const std::filesystem::path& results_file,
        const std::string& extra_parameters)
    {
        std::string phmmer_output = "hmmer_output.fasta";
        std::string phmmer_output_table = "hmmer_output_table.fasta";

        std::vector<std::string> command_tokens{ std::string{"phmmer"},
                                                std::string{"-o /dev/null"},
                                                std::string{"--tblout " + results_file.string()},
                                                extra_parameters,
                                                query_file.string(),
                                                database_file.string() };

        std::string phmmer_command;
        for (const auto& token : command_tokens) {
            const char token_separator{ ' ' };
            phmmer_command += token + token_separator;
        }

        LOG(DEBUG) << "Starting phmmer with the following command: \"" << phmmer_command << "\"";
        int sys_call_result = std::system(phmmer_command.c_str());
        if (sys_call_result != 0) {
            throw std::runtime_error("Failed to find phmmer. Please ensure that the HMMER package is installed and phmmer executable location is set in the PATH variable");
        }
    }

    void buildHmmerProfile(const std::filesystem::path& msa_file,
        const std::filesystem::path& output_file,
        const std::string& extra_parameters)
    {
        const std::vector<std::string> command_tokens{ std::string{"hmmbuild"},
                                                      std::string{"-o /dev/null"},
                                                      extra_parameters,
                                                      output_file.string(),
                                                      msa_file.string() };

        std::string hmmbuild_command;
        for (const auto& token : command_tokens) {
            hmmbuild_command += token + " "; //Space to separate the tokens
        }
        LOG(DEBUG) << "Calling hmmbuild with" << hmmbuild_command;
        int sys_call_result = std::system(hmmbuild_command.c_str()); //std::system expects a C-style string
        if (sys_call_result != 0) {
            throw std::runtime_error("Failed to find hmmbuild. Please ensure that the HMMER package is installed and phmmer executable location is set in the PATH variable");
        }
    }

    void hmmsearch(const std::filesystem::path& profile_file,
        const mipfinder::ProteinList& database,
        const std::filesystem::path& results_file,
        const std::string& extra_parameters)
    {
        const std::filesystem::path results_path = results_file.parent_path();

        std::filesystem::path database_file{ "hmmsearch_database.txt" };
        std::filesystem::path database_file_location = results_path / database_file;

        mipfinder::proteinToFasta(database, database_file_location);

        hmmsearch(profile_file,
            database_file_location,
            results_file,
            extra_parameters);
    }

    void hmmsearch(const std::filesystem::path& profile_file,
        const std::filesystem::path& database_file,
        const std::filesystem::path& results_file,
        const std::string& extra_parameters)
    {
        const std::vector<std::string> command_tokens{ std::string{"hmmsearch"},
                                                      std::string{"-o /dev/null"},
                                                      std::string{"--tblout"},
                                                      results_file.string(),
                                                      profile_file.string(),
                                                      database_file.string(),
                                                      extra_parameters };

        std::string hmmsearch_command;
        for (const auto& token : command_tokens) {
            hmmsearch_command += token + " "; //Space to separate the tokens
        }

        LOG(INFO) << "Starting hmmsearch with the following command: \"" << hmmsearch_command << "\"";
        int sys_call_result = std::system(hmmsearch_command.c_str()); //std::system expects a C-style string
        if (sys_call_result != 0) {
            throw std::runtime_error("Failed to find hmmsearch. Please ensure that the HMMER package is installed and phmmer executable location is set in the PATH variable");
        }
    }

    mipfinder::homology::Results parseHmmerTabularResults(const std::filesystem::path& results_file)
    {
        LOG(DEBUG) << "Parsing tabular HMMER homology search results";
        std::ifstream file{ results_file };
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open " + results_file.string() + ", aborting...");
        }

        mipfinder::homology::Results results;
        std::string line;
        while (std::getline(file, line)) {
            if (line.front() == '#') { //'#' lines are comments
                continue;
            }
            const auto tokens = mipfinder::tokenise(line, ' ');
            const std::string& query = tokens[2];
            const std::string& target = tokens[0];
            const double& bitscore = stod(tokens[5]);

            results.push_back(mipfinder::homology::Result{ .query = query, .target = target, .bitscore = bitscore });
        }

        //Sort the targets alphabetically first, and then by the bitscore to maintain the function return conditions.
        auto sorting_fn = [](const Result& lhs, const Result& rhs) {
            return lhs.query < rhs.query || lhs.query == rhs.query && lhs.bitscore > rhs.bitscore;
        };
        std::sort(std::begin(results), std::end(results), sorting_fn);

        LOG(DEBUG) << "Done parsing homology search results";
        return results;
    }

    mipfinder::homology::Results keepTopHomologues(const mipfinder::homology::Results& homology_search_results,
        const std::size_t maximum_homologues_allowed)
    {
        std::string currently_processed_query = homology_search_results[0].query; //Initialise to the first element to save a check for empty in the loop
        std::size_t proccesed_per_entry = 0;
        mipfinder::homology::Results filtered_results;
        for (const auto& result : homology_search_results) {
            if (currently_processed_query == result.query && proccesed_per_entry != maximum_homologues_allowed) {
                filtered_results.push_back(result);
                ++proccesed_per_entry;
            }
            else {
                currently_processed_query == result.query;
                filtered_results.push_back(result);
                proccesed_per_entry = 1;
            }
        }
        return filtered_results;
    }

    mipfinder::homology::Results filterByBitscore(const mipfinder::homology::Results& homology_results,
        const double minimum_bitscore,
        const double maximum_bitscore)
    {
        LOG(DEBUG) << "Filtering homology result by bitscore cutoffs";
        mipfinder::homology::Results filtered_results;
        for (const auto& result : homology_results)
        {
            if (result.bitscore < minimum_bitscore || result.bitscore > maximum_bitscore) {
                continue;
            }
            filtered_results.push_back(result);
        }
        return filtered_results;
    }

    mipfinder::homology::Results removeSelfHits(const mipfinder::homology::Results& homology_search_results)
    {
        mipfinder::homology::Results filtered;
        for (const auto& result : homology_search_results) {
            if (result.query == result.target) {
                continue;
            }
            filtered.push_back(result);
        }
        return filtered;
    }

}

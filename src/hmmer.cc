#include <fstream>
#include <stdexcept>
#include <unordered_set>

#include "hmmer.h"
#include "helpers.h"

namespace mipfinder::homology
{
    bool operator==(const Result& lhs, const Result& rhs)
    {
        //TODO: Bitscores are double and comparing them for equality needs special treatment.
        return lhs.query == rhs.query
            && lhs.target == rhs.target
            && lhs.bitscore == rhs.bitscore;
    }

    Results::iterator Results::begin()
    {
        return results.begin();
    }

    Results::const_iterator Results::begin() const
    {
        return results.begin();
    }

    Results::const_iterator Results::cbegin() const
    {
        return results.cbegin();
    }

    Results::iterator Results::end()
    {
        return results.end();
    }

    Results::const_iterator Results::end() const
    {
        return results.end();
    }

    Results::const_iterator Results::cend() const
    {
        return results.cend();
    }

    void Results::add(Result result)
    {
        results.insert(result);
    }

    Results parseResultsFile(const std::filesystem::path& results_file)
    {
        LOG(DEBUG) << "Parsing tabular HMMER homology search results";
        std::ifstream file{ results_file };
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open " + results_file.string() + ", aborting...");
        }

        Results results;
        std::string line;
        while (std::getline(file, line)) {
            if (line.front() == '#') { //'#' lines are comments
                continue;
            }
            const auto tokens = mipfinder::tokenise(line, ' ');
            const std::string& query = tokens[2];
            const std::string& target = tokens[0];
            const double& bitscore = stod(tokens[5]);

            results.add(mipfinder::homology::Result{ .bitscore = bitscore, .query = query, .target = target, });
        }

        LOG(DEBUG) << "Done parsing homology search results";
        return results;
    }

    Results keepTopHomologues(const Results& homology_search_results,
                              const std::size_t maximum_homologues_allowed)
    {
        std::size_t proccesed_per_entry = 0;
        std::unordered_set<std::string> processed_entries;
        Results filtered_results;
        for (const auto& result : homology_search_results) {
            if (processed_entries.contains(result.query) && proccesed_per_entry != maximum_homologues_allowed) {
                filtered_results.add(result);
                ++proccesed_per_entry;
            }
            else {
                processed_entries.insert(result.query);
                filtered_results.add(result);
                proccesed_per_entry = 1;
            }
        }
        return filtered_results;
    }

    Results filterByBitscore(const Results& homology_search_results,
                             const double minimum_bitscore,
                             const double maximum_bitscore)
    {
        LOG(DEBUG) << "Filtering homology result by bitscore cutoffs";
        Results filtered_results;
        for (const auto& result : homology_search_results)
        {
            if (result.bitscore >= minimum_bitscore && result.bitscore <= maximum_bitscore) {
                filtered_results.add(result);
            }
        }
        return filtered_results;
    }

    Results removeSelfHits(const Results& homology_search_results)
    {
        LOG(DEBUG) << "Filtering self-hits from homology results";
        Results filtered;
        for (const auto& result : homology_search_results) {
            if (result.query != result.target) {
                filtered.add(result);
            }
        }
        return filtered;
    }
}

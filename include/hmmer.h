#ifndef MIPFINDER_HMMER_H
#define MIPFINDER_HMMER_H

#include <filesystem>
#include <functional>
#include <limits>
#include <ranges>

#include "aliases.h"
#include "interpro.h"
#include "helpers.h"
#include "protein.h"

namespace mipfinder::homology
{
    struct Result
    {
        std::string query;
        std::string target;
        double bitscore;
    };

    struct Homologue
    {
        std::string target;
        double bitscore;
    };

    //Invariant: The homology search results are always ordered by bitscore for a given query, highest to lowest.
    using Results = std::unordered_map<std::string, std::vector<mipfinder::homology::Homologue>>;

    void phmmer(const std::filesystem::path& query_file,
        const std::filesystem::path& database_file,
        const std::filesystem::path& results_file,
        const std::string& extra_parameters);


    void buildHmmerProfile(const std::filesystem::path& msa_file,
        const std::filesystem::path& output_file,
        const std::string& extra_parameters = "");

    void hmmsearch(const std::filesystem::path& profile_file,
        const ProteinSet& database,
        const std::filesystem::path& output_file,
        const std::string& extra_parameters = "");

    void hmmsearch(const std::filesystem::path& profile_file,
        const std::filesystem::path& database_file,
        const std::filesystem::path& output_file,
        const std::string& extra_parameters = "");

    /* Reads in a results file in table format (`--tblout` from HMMMER) */
    mipfinder::homology::Results parseResults(const std::filesystem::path& results_file);

    //Filter the homology search results to only contain those key-value pairs
    //that have equal or less than 'maximum_homologues_allowed' entries.
    mipfinder::homology::Results keepTopHomologues(mipfinder::homology::Results& homology_search_results,
                                                   const std::size_t maximum_homologues_allowed);

    //Find the corresponding proteins from the homology search results
    template <typename Cont>
    requires std::ranges::range<Cont>
        Cont findCorrespondingProteins(const mipfinder::homology::Results& results,
            const Cont& proteome)
    {
        /* Lookup table for fast searching */
        std::unordered_map<std::string, mipfinder::Protein> lookup_table;
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


    template <typename T>
    requires std::ranges::range<T>/*&& requires (std::ranges::range_value_t<T> t)
    {
        { t.bitscore } -> std::convertible_to<double>;
    }*/
    mipfinder::homology::Results filterByBitscore(const T& homology_results,
        const double minimum_bitscore = (std::numeric_limits<double>::min)(), //Min and max have to be wrapped in 
        const double maximum_bitscore = (std::numeric_limits<double>::max)()) //parentheses due to unwanted macro expansion
    {
        LOG(DEBUG) << "Filtering homology result by bitscore cutoffs";
        mipfinder::homology::Results filtered_results;
        for (const auto& [protein_identifier, homologues] : homology_results)
        {
            for (const auto& homologue : homologues) {
                if (homologue.bitscore < minimum_bitscore || homologue.bitscore > maximum_bitscore) {
                    continue;
                }
                filtered_results[protein_identifier].push_back(homologue);
            }
        }
        return filtered_results;
    }

    ///* Keeps up to @hits_to_keep of best-scoring HMMER results for each query */
    //mipfinder::homology::Results
    //	keepTopHits(const mipfinder::homology::Results& results, std::size_t hits_to_keep);

    ///* Removes all entries from HMMER results files where the query is the same
    // * as the target */
    //mipfinder::homology::Results
    //	removeSelfHits(const mipfinder::homology::Results& results);

    ///* Returns a list of proteins from @proteins that match the identifiers found
    // * in HMMER @results */
    //mipfinder::ProteinSet
    //	convertToProtein(const mipfinder::homology::Results& results,
    //					 const mipfinder::ProteinSet& proteins);

    ////Creates HMMER profiles from all cMIP Multiple Sequence Alignments. 
    ////Creates one HMMER profile per one MSA.
    //void createHmmProfiles(const std::filesystem::path& msa_dir,
    //					   const std::filesystem::path& output_dir);

    ////Creates one large HMMER profile file from all homologous cMIP HMMER
    ////profiles. This simplifies downstream processing.
    //std::filesystem::path
    //	concatenateHmmProfiles(const std::filesystem::path& hmmer_profile_dir,
    //						   const std::filesystem::path& output_file);

    ///* Returns all HmmerResults queries where the target length - query length >=
    // * @min_length_difference */
    //mipfinder::homology::Results
    //	filterByLengthDifference(const mipfinder::homology::Results& results,
    //							 const mipfinder::ProteinSet& proteins,
    //							 unsigned int min_length_difference);

    ///* Takes HmmerResults and filters out all targets that do not contain the same
    // * protein family identifier as the query */
    //mipfinder::homology::Results
    //	filterByProteinFamilyAndDomain(const mipfinder::homology::Results& results,
    //								   const mipfinder::Proteome& proteome);

    ///* Returns all HmmerResults queries that have less or equal than
    // * @maximum_homologues targets */
    //mipfinder::homology::Results
    //	filterByHomologueCount(const mipfinder::homology::Results& results,
    //						   unsigned int maximum_homologues);

}

namespace std
{
    template <>
    struct hash<mipfinder::homology::Result>
    {
        std::size_t operator()(const mipfinder::homology::Result& k) const
        {
            using std::hash;

            return ((hash<std::string>()(k.query)
                ^ (hash<std::string>()(k.target) << 1)) >> 1);
        }
    };
}
#endif

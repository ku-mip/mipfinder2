#ifndef MIPFINDER_HMMER_H
#define MIPFINDER_HMMER_H

#include <ranges>
#include <filesystem>
#include <functional>
#include <limits>
#include <ranges>
#include <string>

#include "interpro.h"
#include "helpers.h"
#include "protein.h"

namespace mipfinder::homology
{
    //Call phmmer with the given `query_file` and `database_file`, and direct the output in a tabular format into `results_file`.
    //This corresponds to calling "phmmer [options] --tblout results_file query_file database_file
    template <typename Options>
        requires std::ranges::range<Options> &&
                 std::convertible_to<std::ranges::range_value_t<Options>, std::string>
    void phmmer(const std::filesystem::path& query_file,
                const std::filesystem::path& database_file,
                const std::filesystem::path& results_file,
                const Options& options)
    {
        std::string phmmer_options;
        static constexpr auto option_separator = ' ';
        for (const auto& option : options)
        {
            phmmer_options += option += option_separator;
        }
        std::string program_name = "phmmer";
        std::string phmmer_command = program_name + option_separator + phmmer_options + query_file.string() + database_file.string();

        LOG(DEBUG) << "Starting phmmer with the following command: \"" << phmmer_command << "\"";
        int sys_call_result = std::system(phmmer_command.c_str());
        if (sys_call_result != 0) {
            throw std::runtime_error("Failed to find phmmer. Please ensure that the HMMER package is installed and phmmer executable location is set in the PATH variable");
        }

    }

    template <typename Options>
        requires std::ranges::range<Options>
    void buildHmmerProfile(const std::filesystem::path& msa_file,
                           const std::filesystem::path& output_file,
                           const Options& options);

    //void hmmsearch(const std::filesystem::path& profile_file,
    //    const ProteinSet& database,
    //    const std::filesystem::path& output_file,
    //    const std::string& extra_parameters = "");

    void hmmsearch(const std::filesystem::path& profile_file,
        const std::filesystem::path& database_file,
        const std::filesystem::path& output_file,
        const std::string& options = "");

    struct Result
    {
        double bitscore;
        std::string query;
        std::string target;
    };

    using Results = std::vector<Result>;
    /* Parse a HMMER homology search result file that was written using the --tblout specifier
     *
     * Return a vector of Result objects whose relative ordering is indeterminate, except that every Result
     * object with the same 'query' value will be adjacent to each other, and within these Result objects the
     * relative ordering is based on their bitscore value in non-ascending order. In other words, while the queries
     * can appear in any order relative to each other, all homology search results related to a specific query are
     * grouped together with the highest-scoring (highest bitscore) appearing first.
     */
    mipfinder::homology::Results parseHmmerTabularResults(const std::filesystem::path& results_file);

    //Convenience functions

    /* Keep 'maximum_homologues_allowed' per unique query in the 'homology_search_results', discarding the rest. */
    mipfinder::homology::Results keepTopHomologues(const mipfinder::homology::Results& homology_search_results,
        std::size_t maximum_homologues_allowed);

    /* Remove all homology results where the homology bitscore is less than 'minimum_bitscore' 
     * or more than 'maximum_bitscore'. */
    mipfinder::homology::Results filterByBitscore(const mipfinder::homology::Results& homology_results,
        double minimum_bitscore = (std::numeric_limits<double>::min)(),
        double maximum_bitscore = (std::numeric_limits<double>::max)());

    /* Remove any homology search result where the query and the target are the same. */
    mipfinder::homology::Results removeSelfHits(const mipfinder::homology::Results& homology_search_results);



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

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


namespace
{


    /**
     * @brief  Join tokens in an iterable collection into a string.
     * @return  A string consisting of individual elements in the @a collection 
     *          separated by @a delimiter. If @a collection is empty, returns an
     *          empty string.
     */
    template <typename T>
    requires std::ranges::range<T>
          && std::convertible_to<std::ranges::range_value_t<T>, std::string>
    std::string join(const char delimiter, const T& collection)
    {
        std::string joined_elements;
        for (const auto& elem : collection)
        {
            joined_elements += elem += delimiter;
        }
        return joined_elements;
    }

    template <typename T>
    std::string join(char delimiter, T elem)
    {
        return std::string{ elem };
    }

    template <typename T, typename... Args>
    std::string join(char delimiter, T elem, Args... args)
    {
        return std::string{ elem } + delimiter + join(delimiter, args...);
    }
}

namespace mipfinder::homology
{
    template <typename T>
    concept IsOption = std::ranges::range<T>
                    && std::convertible_to<std::ranges::range_value_t<T>, std::string>;

    template <typename T>
        requires IsOption<T>
    void phmmer(const std::filesystem::path& input_query_file,
                const std::filesystem::path& input_database_file,
                const T& options)
    {
        constexpr auto phmmer_options = join(' ', options);
        constexpr auto program_name = "phmmer";
        constexpr char delimiter = ' ';
        auto phmmer_command = join(delimiter, program_name, phmmer_options, input_query_file.string(), input_database_file.string());

        LOG(DEBUG) << "Starting phmmer with the following command: \"" << phmmer_command << "\"";
        int sys_call_result = std::system(phmmer_command.c_str());
        if (sys_call_result < 0) {
            throw std::runtime_error("Failed to find phmmer. Please ensure that the HMMER package is installed and phmmer executable location is set in the PATH variable");
        }
    }

    template <typename T>
        requires IsOption<T>
    void buildHmmerProfile(const std::filesystem::path& input_msa_file,
                           const std::filesystem::path& output_file,
                           const T& options)
    {
        std::string hmmbuild_options = join(' ', options);
        std::string program_name = "hmmbuild";
        constexpr char delimiter = ' ';
        auto hmmbuild_command = join(delimiter, program_name, hmmbuild_options, output_file.string(), input_msa_file.string());

        LOG(DEBUG) << "Calling hmmbuild with" << hmmbuild_command;
        int sys_call_result = std::system(hmmbuild_command.c_str()); //std::system expects a C-style string
        if (sys_call_result < 0) {
            throw std::runtime_error("Failed to find hmmbuild. Please ensure that the HMMER package is installed and hmmbuild executable location is set in the PATH variable");
        }
    }

    template <typename T>
        requires IsOption<T>
    void hmmsearch(const std::filesystem::path& hmmer_profile_file,
                   const std::filesystem::path& sequence_database_file,
                   const T& options)
    {
        std::string hmmsearch_options = join(' ', options);
        std::string program_name = "hmmsearch";
        constexpr char delimiter = ' ';
        auto hmmsearch_command = join(delimiter, program_name, hmmsearch_options, hmmer_profile_file.string(), sequence_database_file.string());

        LOG(INFO) << "Calling hmmsearch with " << hmmsearch_command;
        int sys_call_result = std::system(hmmsearch_command.c_str()); //std::system expects a C-style string
        if (sys_call_result < 0) {
            throw std::runtime_error("Failed to find hmmsearch. Please ensure that the HMMER package is installed and phmmer executable location is set in the PATH variable");
        }
    }


    struct Result
    {
        double bitscore;
        std::string query;
        std::string target;
    };

    //class Results
    //{
    //public:
    //    Results(std::filesystem::path& homology_search_results_file);
    //private:
    //    std::set<Result> results;
    //};

    auto ResultComparator = [](const Result& lhs, const Result& rhs)
    {
        return (lhs.query < rhs.query)
            || (lhs.query == rhs.query && lhs.bitscore > rhs.bitscore);
    };

    using Results = std::set<Result, decltype(ResultComparator)>;

    /* Parse a HMMER homology search result file that was written using the --tblout specifier
     *
     * Return a vector of Result objects whose relative ordering is indeterminate, except that every Result
     * object with the same 'query' value will be adjacent to each other, and within these Result objects the
     * relative ordering is based on their bitscore value in non-ascending order. In other words, while the queries
     * can appear in any order relative to each other, all homology search results related to a specific query are
     * grouped together with the highest-scoring (highest bitscore) appearing first.
     */
    Results parseHmmerTabularResults(const std::filesystem::path& results_file);

    //Convenience functions

    /* Keep 'maximum_homologues_allowed' per unique query in the 'homology_search_results', discarding the rest. */
    Results keepTopHomologues(const Results& homology_search_results,
        std::size_t maximum_homologues_allowed);

    /* Remove all homology results where the homology bitscore is less than 'minimum_bitscore' 
     * or more than 'maximum_bitscore'. */
    Results filterByBitscore(const Results& homology_results,
        double minimum_bitscore = (std::numeric_limits<double>::min)(),
        double maximum_bitscore = (std::numeric_limits<double>::max)());

    /* Remove any homology search result where the query and the target are the same. */
    Results removeSelfHits(const Results& homology_search_results);



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

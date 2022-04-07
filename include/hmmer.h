#ifndef MIPFINDER_HMMER_H
#define MIPFINDER_HMMER_H

#include <filesystem>
#include <limits>
#include <ranges>
#include <set>
#include <string>

#include "easylogging++.h"

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

    auto ResultComparator = [](const Result& lhs, const Result& rhs)
    {
        return (lhs.query < rhs.query)
            || (lhs.query == rhs.query && lhs.bitscore > rhs.bitscore);
    };

    class Results
    {
    private:
        using ResultsContainer = std::set<Result, decltype(ResultComparator)>;
        ResultsContainer results;
    public:
        using iterator = ResultsContainer::iterator;
        using const_iterator = ResultsContainer::const_iterator;
        using reference = ResultsContainer::reference;
        using const_reference = ResultsContainer::const_reference;

        iterator begin();
        const_iterator begin() const;
        const_iterator cbegin() const;

        iterator end();
        const_iterator end() const;
        const_iterator cend() const;

        void add(Result result);
        iterator erase(const_iterator pos);
        iterator erase(const_iterator first, const_iterator last);
    };

    ///* Parse a HMMER homology search result file that was written using the --tblout specifier
    // *
    // * Return a Result objects whose relative ordering is indeterminate, except that every Result
    // * object with the same 'query' value will be adjacent to each other, and within these Result objects the
    // * relative ordering is based on their bitscore value in descending order. In other words, while the queries
    // * can appear in any order relative to each other, all homology search results related to a specific query are
    // * grouped together with the highest-scoring (highest bitscore) appearing first.
    // */
    Results parseResultsFile(const std::filesystem::path& results_file);


    /**
     *  @brief  Extracts top n homologues from the results.
     *  @pre  Each entry with the same query identifier in @a homology_search_results must be sorted
     *        in descending order based on their bitscore.
     */
    Results keepTopHomologues(const Results& homology_search_results,
                              std::size_t maximum_homologues_allowed);

    /**
     *  @brief  Extracts all results where @a minimum_bitscore <= bitscore <= @a maximum_bitscore
     */
    Results filterByBitscore(const Results& homology_results,
                             double minimum_bitscore = (std::numeric_limits<double>::min)(),
                             double maximum_bitscore = (std::numeric_limits<double>::max)());

    /**
     *  @brief  Removes any homology search result where the query and the target are the same.
     */
    Results removeSelfHits(const Results& homology_search_results);


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

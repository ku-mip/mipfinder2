#ifndef MIPFINDER_FASTA_H
#define MIPFINDER_FASTA_H

#include <array>
#include <filesystem>
#include <string>
#include <unordered_map>

namespace mipfinder::fasta
{
    /**
     * @brief  Denotes an individual FASTA record.
     */
    struct Entry
    {
        std::string header;
        std::string sequence;
    };

    using Entries = std::vector<mipfinder::fasta::Entry>;

    /**
     * @brief  Extracts all entries from a FASTA file.
     * @param  file  Path to the FASTA file.
     * @throw  std::runtime_error  If @a file cannot be opened.
     * @return  A collection of FASTA entries.
     */
    Entries parse(const std::filesystem::path& file);
    Entries parse(std::ifstream& stream);

    std::array<std::string, 4> extractUniprotHeader(const std::string& header);
}
#endif


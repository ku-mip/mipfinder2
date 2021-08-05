#ifndef MIPFINDER_FASTA_H
#define MIPFINDER_FASTA_H

#include <array>
#include <filesystem>
#include <string>
#include <unordered_map>


namespace mipfinder::fasta
{
    typedef std::unordered_map<std::string, std::string> Records;

    Records extractRecords(const std::filesystem::path& file);
    Records extractRecords(std::ifstream& stream);

    std::array<std::string, 4> extractUniprotHeader(const std::string& header);
}
#endif


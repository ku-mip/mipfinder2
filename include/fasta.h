#ifndef MIPFINDER_FASTA_H
#define MIPFINDER_FASTA_H

#include <array>
#include <filesystem>
#include <string>
#include <unordered_map>

#include "aliases.h"

namespace mipfinder::fasta
{
  FastaRecords extractRecords(const std::filesystem::path& file);
  FastaRecords extractRecords(std::ifstream& stream);

  std::array<std::string, 4> extractUniprotHeader(const std::string& header);
}
#endif


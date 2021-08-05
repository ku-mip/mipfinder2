#include <iostream>
#include <filesystem>
#include <fstream>
#include <regex>

#include "fasta.h"
#include "file.h"
#include "protein.h"

/* Takes a UniProt database FASTA header as input and returns an array with the
 * following elements:
 * 1: UniProt accession name
 * 2: Entry sequence version
 * 3: Entry description
 * 4: Entry protein existence level
 *
 * If the header cannot be successfully parsed, returns an empty array.
 */
std::array<std::string, 4>
mipfinder::fasta::extractUniprotHeader(const std::string& header)
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
        return std::array<std::string, 4>{};
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
  }

  std::array<std::string, 4> header_contents{accession_name,
                                             sequence_version,
                                             description,
                                             existence_level};
  return header_contents;
}

mipfinder::FastaRecords
mipfinder::fasta::extractRecords(const std::filesystem::path& file)
{
  auto stream = mipfinder::file::open(file);
  return mipfinder::fasta::extractRecords(stream);
}

mipfinder::FastaRecords mipfinder::fasta::extractRecords(std::ifstream& stream)
{
  /* This ensures we are only parsing lines that are part of a FASTA record */
  bool fasta_record_found = false;
  
  std::string fasta_header{};
  std::string fasta_sequence{};

  /* Map containing FASTA headers as keys and sequences as values */
  mipfinder::FastaRecords all_fasta_records;

  std::string line{""};
  while(std::getline(stream, line)) {
    if (line.empty()) {
        continue;
    }

    //If the FASTA header contains no information other than the '>' character marking it as a header,
    //ignore the line.
    if (line.front() == '>' && line.length() == 1) {
        continue;
    }

    if ((line.front() == '>') && (!fasta_record_found)) {
      fasta_record_found = true;
      fasta_header = line;
    }

    else if (!(line.front() == '>') && fasta_record_found) {
      fasta_sequence += line;
    }

    else if ((line.front() == '>') && fasta_record_found) {
      /* If we have found a new record but the current record did not have a 
       * sequence, ignore it */
      if (fasta_sequence.empty()) {
        fasta_sequence.clear();
        fasta_header = line;
        continue;
      }
      else {
        all_fasta_records.emplace(fasta_header, fasta_sequence);

        /* Reset the variable contents */
        fasta_sequence.clear();
        fasta_header = line;
      }
    }

    else if (!fasta_record_found) {
      /* If no record is found or being processed, the data is bad */
      continue;
    }
  }

  /* If final sequence is empty, ignore the header */
  if (fasta_sequence.empty()) {
    return all_fasta_records;
  }

  all_fasta_records.emplace(fasta_header, fasta_sequence);
  return all_fasta_records;
}

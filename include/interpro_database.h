#ifndef MIPFINDER_INTERPRO_H
#define MIPFINDER_INTERPRO_H

#include <unordered_set>
#include <string>

#include "aliases.h"

namespace mipfinder {
  class Configuration;
  class Proteome;
}

namespace mipfinder::interpro {

  typedef std::unordered_set<std::string> Identifiers;
  typedef std::vector<mipfinder::interpro::Result> Database;

  /* Returns a list of all InterPro IDs that map to Uniprot IDs.
   *
   * @data - A handle to a tsv file of an InterPro database file. Column 1
   * contains the UniProt Accession number, column 2 contains the sequence 
   * version of the specific protein, and column three is a comma-separated
   * list of all InterPro identifiers associated with that sequence version */
  mipfinder::interpro::Database parse(std::ifstream& database);

  void assignScores(const mipfinder::hmmer::Results& hmmer_results, 
                    const mipfinder::InterProResults& interpro_results, 
                    const mipfinder::Proteome& proteome,
                    const mipfinder::Configuration& config);

  /* Returns a list of InterPro identifiers that are of type @type 
   * 
   * @interpro_to_type - A tsv file where the first column is the InterPro ID 
   *                     and the second column is its @type */
  Identifiers
  findTypeIdentifiers(const std::filesystem::path interpro_id_types,
                      const std::string& type);
  Identifiers
  findTypeIdentifiers(std::ifstream& interpro_id_types,
                      const std::string& type);
}
#endif

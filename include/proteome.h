#ifndef MIPFINDER_PROTEOME_H
#define MIPFINDER_PROTEOME_H

#include <filesystem>
#include <limits>
#include <memory>
#include <string>

#include "aliases.h"
#include "protein.h"

namespace mipfinder {
  class Proteome {
  public:
    Proteome(const std::filesystem::path& fasta_file);

    //Returns the number of proteins in the proteome
    std::size_t size() const;

    /* Returns all the proteins present in the proteome */
    ProteinSet data() const;

    /* Returns a handle to a protein with a given identifier. If no such protein
     * can be found, returns a nullptr */
    Protein* find(std::string identifier) const;

  private:
    typedef std::unordered_map<std::string, std::unique_ptr<Protein>> Entries;
    Entries proteins_;
  };

  /* Associate GO identifiers with proteins in the @proteome */
  void addGoIdentifiers(const mipfinder::Proteome& proteome,
                        const mipfinder::Go& go_database,
                        const std::filesystem::path& uniprot_to_go);
                        
  /* Associate InterPro identifiers with proteins in the @proteome */
  void addInterproIdentifiers(const mipfinder::Proteome& proteome,
                        const mipfinder::Interpro& interpro_database,
                        const std::filesystem::path& uniprot_to_go);

  static constexpr std::size_t SEQUENCE_MAX_LENGTH = (std::numeric_limits<std::size_t>::max)();

  /* Returns a list of proteins between [min_length, max_length] */
  ProteinSet filterByLength(const ProteinSet& proteins, std::size_t min_length = 0, 
      std::size_t max_length = SEQUENCE_MAX_LENGTH);

  /* Returns a list of proteins between [min_length, max_length] */
  ProteinSet filterByLength(const Proteome& proteome, std::size_t min_length = 0, 
      std::size_t max_length = SEQUENCE_MAX_LENGTH);

  /* Returns  list of proteins with a protein existence level <= @level_cutoff */
  ProteinSet filterByExistence(const ProteinSet& proteome, int level_cutoff);
}
#endif
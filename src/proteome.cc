#include <algorithm>
#include <cassert>

#include "fasta.h"
#include "file.h"
#include "helpers.h"
#include "protein.h"
#include "proteome.h"

#include <iostream>

namespace
{
  /* Parses UniProt to GO terms mapping file. Returns a map where the keys are
   * Protein IDs (ENTRY_NAME.SEQUENCE_VERSION), and the values are a list of
   * GO identifiers */
  std::unordered_map<std::string, std::vector<std::string>>
  parseUniProtToGo(const std::filesystem::path& uniprot_to_go)
  {
    std::unordered_map<std::string, std::vector<std::string>> results;
    std::ifstream f;
    f.open(uniprot_to_go);
    std::string line;
    while (getline(f, line)) {
      const auto tokens = mipfinder::tokenise(line, '\t');
      /* If a UniProt entry name does not have any GO IDs associated with it */
      if (tokens.size() != 3) {
        continue;
      }

      const std::string uniprot_accession_name = tokens[0];
      const std::string sequence_version = tokens[1];
      const std::string go_identifiers = tokens[2];

      /* mipfinder stores Protein IDs with their sequence version */
      const std::string protein_id = uniprot_accession_name + "." + sequence_version;

      auto go_id_tokens = mipfinder::tokenise(go_identifiers, ';');

      /* The GO IDs come separated by `; `, so we remove all spaces from the
       * tokens */
      for (auto& token : go_id_tokens) {
        token.erase(std::remove_if(token.begin(), token.end(), isspace), token.end());
        results[protein_id].push_back(token);
      }
    }
    return results;
  }

  /* Parses UniProt to InterPro identifier mapping file. Returns a map where the 
   * keys are Protein IDs (ENTRY_NAME.SEQUENCE_VERSION), and the values are a 
   * list of InterPro identifiers */
  std::unordered_map<std::string, std::vector<std::string>>
  parseUniProtToInterpro(const std::filesystem::path& uniprot_to_interpro)
  {
    std::unordered_map<std::string, std::vector<std::string>> results;
    auto stream = mipfinder::file::open(uniprot_to_interpro);
    std::string line;
    std::getline(stream, line); /* Skips the header in database file */
    while(std::getline(stream, line)) {
      const auto tokens = mipfinder::tokenise(line, '\t');

      /* Correctly formatted file has three columns of data */
      if (tokens.size() != 3) {
        continue;
      }

      const auto uniprot_id = tokens[0];
      const auto sequence_version = tokens[1];
      const auto full_protein_id = uniprot_id + "." + sequence_version;

      const auto interpro_id_list = tokens[2];
      auto interpro_id_tokens = mipfinder::tokenise(interpro_id_list, ';');

      /* The InterPro IDs come separated by `; `, so we remove all spaces from the
        * tokens */
      for (auto& interpro_id : interpro_id_tokens) {
        interpro_id.erase(std::remove_if(interpro_id.begin(), interpro_id.end(), isspace),
                          interpro_id.end());
                          
        assert(!uniprot_id.empty());
        assert(!interpro_id.empty());

        results[full_protein_id].push_back(interpro_id);
      }
    }
    return results;
  }

}

namespace mipfinder
{ 
  /* Initialise a proteome into a functional state. This includes associating
   * all Protein objects with GO and InterPro identifiers to set them in a 
   * valid state */
  Proteome::Proteome(const std::filesystem::path& fasta_file)
  {
    auto file = file::open(fasta_file);
    const FastaRecords proteome_fasta_records = fasta::extractRecords(file);

    for (const auto& [header, sequence] : proteome_fasta_records) {
      const auto& [protein_id, sequence_version, description, existence_level] 
          = fasta::extractUniprotHeader(header);


      const std::string identifier = protein_id + "." + sequence_version;
      auto p = std::make_unique<Protein>(identifier,
                                         sequence,
                                         description,
                                         std::stoi(existence_level));
      proteins_[identifier] = std::move(p);
    }
  }

  mipfinder::ProteinSet Proteome::data() const
  {
    mipfinder::ProteinSet data;
    for (const auto& [identifier, protein_unique_ptr] : proteins_) {
      data.insert(protein_unique_ptr.get());
    }
    return data;
  }

  std::size_t Proteome::size() const
  {
    return proteins_.size();
  }

  Protein* Proteome::find(std::string identifier) const
  {
    toupper(identifier);

    if (proteins_.count(identifier) == 1) {
      return (proteins_.at(identifier)).get();
    }
    else {
      return nullptr;
    }
  }

  /* Free functions */
  
  //void addGoIdentifiers(const mipfinder::Proteome& proteome,
  //                      const mipfinder::Go& go_database,
  //                      const std::filesystem::path& uniprot_to_go)
  //{
  //  const auto lookup_table = parseUniProtToGo(uniprot_to_go);

  //  for (const auto& [protein_id, go_ids] : lookup_table) {
  //    const auto protein = proteome.find(protein_id);

  //    if (protein == nullptr) {
  //      continue;
  //    }

  //    for (const auto& go_id : go_ids) {
  //      const auto go_entry = go_database.find(go_id);
  //      if (go_entry) {
  //        protein->addGoEntry(go_entry.value());
  //      }
  //    }

  //  }
  //}

  //void addInterproIdentifiers(const mipfinder::Proteome& proteome,
  //                            const mipfinder::Interpro& interpro_database,
  //                            const std::filesystem::path& uniprot_to_go)
  //{
  //  const auto lookup_table = parseUniProtToInterpro(uniprot_to_go);

  //  for (const auto& [protein_id, interpro_ids] : lookup_table) {
  //    const auto protein = proteome.find(protein_id);

  //    if (protein == nullptr) {
  //      continue;
  //    }

  //    assert(protein != nullptr);
  //    for (const auto& interpro_id : interpro_ids) {
  //      const auto interpro_entry = interpro_database.find(interpro_id);
  //      if (interpro_entry) {
  //        protein->addInterproEntry(interpro_entry.value());
  //      }
  //    }
  //  }
  //}

  mipfinder::ProteinSet
  filterByLength(const ProteinSet& proteins,
                 const std::size_t min_length,
                 const std::size_t max_length)
  {
    if (min_length > max_length) {
      throw std::logic_error("Minimum length cannot be larger than maximum length");
    }
    
    ProteinSet filtered_proteins;
    for (auto& protein : proteins) {
      if ((protein->length() >= min_length) && (protein->length() <= max_length)) {
        filtered_proteins.insert(protein);
      }
    }
    return filtered_proteins;
  }

  mipfinder::ProteinSet
  filterByLength(const Proteome& proteome,
                 const std::size_t min_length,
                 const std::size_t max_length)
  {
    if (min_length > max_length) {
      throw std::logic_error("Minimum length cannot be larger than maximum length");
    }
    
    ProteinSet filtered_proteins;
    for (auto& protein : proteome.data()) {
      if ((protein->length() >= min_length) && (protein->length() <= max_length)) {
        filtered_proteins.insert(protein);
      }
    }
    return filtered_proteins;
  }

  mipfinder::ProteinSet
  filterByExistence(const ProteinSet& proteins, int level_cutoff)
  {
    ProteinSet filtered_proteins;
    for (auto& protein : proteins) {
      if (protein->existenceLevel() <= level_cutoff) {
        filtered_proteins.insert(protein);
      }
    }
    return filtered_proteins;
  }


}
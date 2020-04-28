#ifndef MIPFINDER_PROTEIN_H
#define MIPFINDER_PROTEIN_H

#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>

#include "ancestor.h"
#include "aliases.h"
#include "homologue.h"
#include "go.h"
#include "interpro.h"

namespace mipfinder {

  //Protein holds information about one specific protein isoform. The Protein identifier is composed
  //of UniProt accession name plus the sequence version.
  class Protein {
  public:
    enum class Type {
      UNKNOWN,
      CMIP,
      SINGLE_COPY,
      HOMOLOGOUS,
      ANCESTOR,
      KNOWN_MIP
    };

    Protein() = delete;
    Protein(const std::string& identifier,
            const std::string& sequence,
            const std::string& description,
            int protein_existence);
    ~Protein() = default;

    std::string identifier() const;
    std::string sequence() const;
    std::string description() const;

    std::vector<mipfinder::Ancestor> ancestors() const;
    void addAncestor(const Ancestor& ancestor);

    std::vector<mipfinder::Homologue> homologues() const;
    void addHomologue(const Homologue& homologue);

    Type type() const;
    std::string type_to_string() const;
    void setType(const mipfinder::Protein::Type& type);

    mipfinder::Go::Entries goEntries() const;
    void addGoEntry(const mipfinder::Go::Entry& entry);

    mipfinder::Interpro::Entries interproEntries() const;
    void addInterproEntry(const mipfinder::Interpro::Entry& entry);

    int existenceLevel() const;

    std::size_t length() const;
    double score() const;
    void changeScore(double score);
  private:
    std::string identifier_;
    std::string sequence_;
    std::string description_;
    int existence_level_;
    double score_;

    mipfinder::Protein::Type type_;
    std::vector<mipfinder::Ancestor> ancestors_;
    std::vector<mipfinder::Homologue> homologues_;

    mipfinder::Interpro::Entries interpro_entries_;
    mipfinder::Go::Entries go_entries_;
  };

  /* Returns the instability index of the sequence */
  double instability_index(const std::string& sequence);

  bool operator==(const Protein& lhs, const Protein& rhs);

  //Creates a FASTA file from existing Proteins
  template <typename Container>
  void proteinToFasta(const Container& proteins, 
                      const std::filesystem::path& output)
  {
    std::ofstream f;
    f.open(output);
    for (const auto& protein : proteins) {
      f << ">" << protein->identifier() << "\n" << protein->sequence() << "\n";
    }
  }

}
#endif
#ifndef MIPFINEDR_PRINTER_H
#define MIPFINDER_PRINTER_H

#include <filesystem>

namespace {
  class Protein;
  struct ProteinCmp
  {
    bool operator()(const Protein* lhs, const Protein* rhs) const;
  };
}

namespace mipfinder::printer
{
  /* @format specifies what does the final output will look like 
   * %acc - cMIP UniProt accession name
   * %length - cMIP length
   * %score - cMIP score
   * %ancestors - cMIP ancestors
   * %homologues - cMIP homologues
   * %description - cMIP description
   * %gene_ontology - cMIP GO terms
   * %type - cMIP type (unique or homologous)
   * %gene - cMIP gene name
   * %known_mip - Prints "Y" if the detected cMIP is in the known_microproteins.fasta
   * %ancestor_nr - Prints the number of detected ancestors, up to "max_allowed_ancetors"
   * %ancestor_domains - Prints ancestor InterPro domain descriptions
   * %ancestor_go - Prints ancestor GO annotations
   */
  void createReport(const std::string& format,
                    char delimiter,
                    const std::vector<Protein*>& proteins,
                    const std::filesystem::path& output_file);
}

#endif
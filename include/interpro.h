#ifndef MIPFINDER_INTERPRO_H
#define MIPFINDER_INTERPRO_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>
#include <optional>

namespace mipfinder
{
  class Interpro {
  public:
    /* All InterPro identifiers fall into one of these categories */
    enum class Type {
      ACTIVE_SITE,
      BINDING_SITE,
      CONSERVED_SITE,
	  DOMAIN_TYPE,
      FAMILY,
      HOMOLOGOUS_SUPERFAMILY,
      PTM,
      REPEAT
    };

    struct Entry {
      std::string interpro_id;
      std::string entry_name;
      Interpro::Type type;
      auto operator<=>(const Interpro::Entry&) const = default;
    };

    typedef std::vector<Entry> Entries;
    typedef std::unordered_map<std::string, Entry> Data;

    /* @database_file is a tsv-file that specifies which type each InterPro id is
     * Column 1: InterPro identifier
     * Column 2: Identifier type
     * Column 3: Identifier description
     */
    Interpro(const std::filesystem::path& database_file);

    /* Returns the entry corresponding to the @interpro_identifier. If no such 
     * entry exists, std::optional::value() returns false */
    std::optional<Interpro::Entry> find(const std::string& interpro_identifier) const;
  private:
    Data interpro_entries_;
  };
}
#endif

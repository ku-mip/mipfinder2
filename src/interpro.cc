#include <algorithm>
#include <cassert>
#include <filesystem>
#include <iostream>

#include "file.h"
#include "helpers.h"
#include "interpro.h"

namespace
{
  mipfinder::Interpro::Data
  parseDatabase(const std::filesystem::path& interpro_database)
  {
    mipfinder::Interpro::Data results;
    auto stream = mipfinder::file::open(interpro_database);
    std::string line;
    std::getline(stream, line); /* Skips the header in database file */
    while(std::getline(stream, line)) {
      const auto tokens = mipfinder::tokenise(line, '\t');

      /* Correctly formatted file has three columns of data */
      if (tokens.size() != 3) {
        continue;
      }

      const auto interpro_id = tokens[0];
      const auto id_type = tokens[1];
      const auto entry_name = tokens[2];

      mipfinder::Interpro::Type type;
      if (id_type == "Active_site") {
        type = mipfinder::Interpro::Type::ACTIVE_SITE;
      }
      else if (id_type == "Binding_site") {
        type = mipfinder::Interpro::Type::BINDING_SITE;
      }
      else if (id_type == "Conserved_site") {
        type = mipfinder::Interpro::Type::CONSERVED_SITE;
      }
      else if (id_type == "Domain") {
        type = mipfinder::Interpro::Type::DOMAIN_TYPE;
      }
      else if (id_type == "Family") {
        type = mipfinder::Interpro::Type::FAMILY;
      }
      else if (id_type == "Homologous_superfamily") {
        type = mipfinder::Interpro::Type::HOMOLOGOUS_SUPERFAMILY;
      }
      else if (id_type == "PTM") {
        type = mipfinder::Interpro::Type::PTM;
      }
      else if (id_type == "Repeat") {
        type = mipfinder::Interpro::Type::REPEAT;
      }
      else {
        continue;
      }
       
      results[interpro_id] = mipfinder::Interpro::Entry{interpro_id,
                                                        entry_name,
                                                        type};
    }
    return results;
  }

}

namespace mipfinder
{
  Interpro::Interpro(const std::filesystem::path& database_file)
                     : interpro_entries_(parseDatabase(database_file)) {}

  std::optional<Interpro::Entry> Interpro::find(const std::string& interpro_identifier) const
  {
    if (interpro_entries_.count(interpro_identifier) == 1) {
      return interpro_entries_.at(interpro_identifier);
    }
    return {};
  }

  bool operator==(const Interpro::Entry& lhs, const Interpro::Entry& rhs)
  {
    return lhs.interpro_id == rhs.interpro_id &&
           lhs.entry_name == rhs.entry_name &&
           lhs.type == rhs.type;
  }
}
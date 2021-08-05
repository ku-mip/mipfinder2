#include <filesystem>
#include <fstream>
#include <optional>
#include <string>

#include "easylogging++.h"
#include "helpers.h"
#include "go.h"

#include <iostream>

namespace {
  /* Reads in the GO database and creates GO entries. */
  mipfinder::Go::Data parseDatabase(const std::filesystem::path& database_file)
  {
    /* Store relevant terms while parsing a GO record */
    std::string go_id;
    std::string go_process_name;
    std::string go_process_type;

    mipfinder::Go::Data go_entries;
    bool go_record_found = false;
    std::ifstream f{database_file};
    std::string line;
    while (getline(f, line)) {
      /* If we have reached a new record while parsing current record */
      if (go_record_found && line == "[Term]") {
        /* Construct a new GO object and reset variables. Also add it to the
         * lookup table for quick searching */
        mipfinder::Go::Entry go_entry{go_id, go_process_name, go_process_type};
        go_entries[go_id] = go_entry;

        go_id.clear();
        go_process_name.clear();
        go_process_type.clear();
        continue;
      }

      else if (line == "[Term]") {
        go_record_found = true;
        continue;
      }

      else if (!go_record_found) {
        continue;
      }

      /* Separate `parameter_name: value` into `parameter_name` and `value` */
      const auto separator_pos = line.find_first_of(':');

      if (separator_pos == std::string::npos) {
        continue;
      }

      const std::string parameter = line.substr(0, separator_pos);
      /* After separator there is always a space, so adding + 2 gets the first 
       * character after the space */
      const std::string value = line.substr(separator_pos + 2, line.size());

      if (parameter == "id") {
        go_id = value;
      }
      else if (parameter == "name") {
        go_process_name = value;
      }
      else if (parameter == "namespace") {
        go_process_type = value;
      }
      else {
        continue;
      }
    }
    LOG(DEBUG) << "Created " << go_entries.size() << " GO records";
    return go_entries;
  }
}

namespace mipfinder
{
  Go::Go(const std::filesystem::path& go_database)
    : go_entries_(parseDatabase(go_database)) { }

  std::optional<Go::Entry> Go::find(const std::string& go_identifier) const
  {
    if (go_entries_.count(go_identifier) == 1) {
      return go_entries_.at(go_identifier);
    }
    return {};
  }
}
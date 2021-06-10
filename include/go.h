#ifndef MIPFINDER_GO_H
#define MIPFINDER_GO_H

#include <filesystem>
#include <string>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace mipfinder
{
  class Go {
  public:
    struct Entry {
      std::string identifier;
      std::string process_name;
      std::string process_type;
      auto operator<=>(const Entry&) const = default;
    };

    typedef std::vector<Entry> Entries;
    typedef std::unordered_map<std::string, Entry> Data;

    /* Initialises the Go database */
    Go(const std::filesystem::path& go_database);

    /* Returns the entry corresponding to the @go_identifier. If no such 
     * entry exists, std::optional::value() returns false */
    std::optional<Go::Entry> find(const std::string& go_identifier) const;
  private:
    Data go_entries_;
  };
}
#endif
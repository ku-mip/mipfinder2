#include <cassert>
#include <filesystem>

#include "aliases.h"
#include "configuration.h"
#include "easylogging++.h"
#include "file.h"
#include "helpers.h"
#include "interpro_database.h"
#include "proteome.h"
#include "protein.h"

namespace mipfinder::interpro
{
  /* Returns the InterPro identifiers for @type types */
  Identifiers
  findTypeIdentifiers(const std::filesystem::path interpro_to_type,
                      const std::string& type)
  {
    Identifiers lookup_table;

    std::ifstream f;
    f.open(interpro_to_type);
    std::string line;
    std::getline(f, line); /* Skips the header in database file */

    while(std::getline(f, line)) {
      const auto tokens = mipfinder::tokenise(line, '\t');
      const auto entry_type = tokens[1];

      if (entry_type == type) {
        const auto interpro_id = tokens[0];
        lookup_table.insert(interpro_id);
      }
    }
    return lookup_table;
  }

  mipfinder::interpro::Database
  parse(std::filesystem::path& database)
  {
    auto stream = mipfinder::file::open(database);

  }

  mipfinder::interpro::Database
  parse(std::ifstream& database)
  {
    mipfinder::InterProResults results;
    std::string line;
    std::getline(database, line); /* Skips the header in database file */
    while(std::getline(database, line)) {
      const auto tokens = mipfinder::tokenise(line, '\t');

      /* Correctly formatted file has three columns of data */
      if (tokens.size() != 3) {
        continue;
      }

      const auto interpro_id_list = tokens[2];
      auto interpro_id_tokens = mipfinder::tokenise(interpro_id_list, ';');

      /* The InterPro IDs come separated by `; `, so we remove all spaces from the
        * tokens */
      for (auto& interpro_id : interpro_id_tokens) {
        interpro_id.erase(remove_if(interpro_id.begin(), interpro_id.end(), isspace),
                          interpro_id.end());
                          
        const auto uniprot_id = tokens[0];
        const auto sequence_version = tokens[1];
        const auto full_id = uniprot_id + "." + sequence_version;

        assert(!uniprot_id.empty());
        assert(!interpro_id.empty());
        results.push_back(mipfinder::interpro::Result{full_id, interpro_id});
      }
    }
    LOG(DEBUG) << "Finished parsing InterPro database";
    LOG(DEBUG) << "Created " << results.size() << " entries\n";
    return results;
  }

  void assignScores(const mipfinder::HmmerResults& hmmer_results, 
                    const mipfinder::InterProResults& interpro_results, 
                    const mipfinder::Proteome& proteome,
                    const mipfinder::Configuration& config)
  {

    if (hmmer_results.size() == 0) {
      LOG(ERROR) << "No records found in HMMER results, aborting...";
      return;
    }
    if (interpro_results.size() == 0) {
      LOG(ERROR) << "No records found in InterPro results, aborting...";
      return;
    }

    //Read InterPro results into a map for quick lookup
    std::unordered_map<std::string, std::vector<std::string>> interpro_lookup;
    LOG(DEBUG) << "Constructing InterPro lookup table";
    for (const auto& i : interpro_results) {
      interpro_lookup[i.uniprot_id].push_back(i.interpro_id);
    }
    LOG(DEBUG) << "Finished";

    static constexpr double INTERPRO_DOMAIN_MATCH_SCORE = 100;

    for (const auto& h_result : hmmer_results) {
      //Don't add scores to self
      if (h_result.query == h_result.target) {
        continue;
      }

      if (interpro_lookup.count(h_result.query) == 0) {
        continue;
      }

      const auto hmmer_query_domains = interpro_lookup[h_result.query];
      const auto hmmer_target_domains = interpro_lookup[h_result.query];  

      for (const auto& domain_id : hmmer_query_domains) {
        const auto search_result = std::find(hmmer_target_domains.begin(), hmmer_target_domains.end(), domain_id);
        if (search_result != hmmer_target_domains.end()) {
          const auto p = proteome.find(h_result.query);
          if (p != nullptr) {
            p->changeScore(INTERPRO_DOMAIN_MATCH_SCORE);
          }
        }
      }
    }
  }

}
#ifndef MIPFINDER_ALIASES_H
#define MIPFINDER_ALIASES_H

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "ancestor.h"
#include "go_record.h"
#include "homologue.h"

namespace mipfinder {
  class Protein;
  class Interpro;

  typedef std::unordered_map<std::string, std::string> FastaRecords;

  typedef std::unordered_map<Protein*, std::unordered_set<Protein*>> ProteinMap;
  typedef std::unordered_set<Protein*> ProteinSet;



  typedef std::vector<mipfinder::Ancestor> Ancestors;
  //typedef std::vector<mipfinder::Result> Homologues;

  typedef std::vector<GoRecord> GoRecords;
}
#endif

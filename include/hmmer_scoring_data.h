#ifndef MIPFINDER_HMMER_SCORING_DATA
#define MIPFINDER_HMMER_SCORING_DATA

#include "protein.h"

namespace mipfinder
{
  struct HmmerScoringData {
    const mipfinder::protein::Protein* query;
    const mipfinder::protein::Protein* target;
    const double bitscore;
  };
}
#endif

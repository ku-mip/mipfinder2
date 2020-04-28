#ifndef MIPFINDER_HMMER_SCORING_DATA
#define MIPFINDER_HMMER_SCORING_DATA

#include "protein.h"

namespace mipfinder
{
  struct HmmerScoringData {
    const Protein* query;
    const Protein* target;
    const double bitscore;
  };
}
#endif

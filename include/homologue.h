#ifndef MIPFINDER_HOMOLOGUE_H
#define MIPFINDER_HOMOLOGUE_H

namespace mipfinder 
{
  class Protein;

  struct Homologue {
    mipfinder::Protein* protein;
    double bitscore;
  };
}
#endif
#ifndef MIPFINDER_ANCESTOR_H
#define MIPFINDER_ANCESTOR_H

namespace mipfinder 
{
  class Protein;

  struct Ancestor {
    mipfinder::Protein* protein;
    double bitscore;
  };
}
#endif
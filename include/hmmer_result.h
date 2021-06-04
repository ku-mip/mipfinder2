#ifndef MIPFINDER_HMMER_RESULT
#define MIPFINDER_HMMER_RESULT

#include <string>

namespace mipfinder::homology {
  struct Result {
    std::string query;
    std::string target;
    double bitscore;
  };

  inline bool operator==(const Result& lhs, const Result& rhs)
  {
    /* Since bitscores only ever have one digit, comparing them for equality 
     * should be safe */
    return lhs.query == rhs.query &&
           lhs.target == rhs.target &&
           lhs.bitscore == rhs.bitscore;
  }

}
#endif

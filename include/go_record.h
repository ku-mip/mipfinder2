#ifndef MIPFINDER_GO_RECORD_H
#define MIPFINDER_GO_RECORD_H

#include <string>

namespace mipfinder 
{
  struct GoRecord {
    std::string identifier;
    std::string process_name;
    std::string process_type;
  };
}
#endif


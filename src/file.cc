#include "file.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace mipfinder::file
{
  std::ifstream open(const std::filesystem::path& filename)
  {
    mipfinder::file::exists(filename);

    std::ifstream f;
    f.open(filename);
    return f;
  }

  void exists(const std::filesystem::path& filename) {
    if (!std::filesystem::exists(filename)) { 
      const std::string error_msg = filename.string() + " could not be found";
      throw std::runtime_error(error_msg);
    }
  }

}
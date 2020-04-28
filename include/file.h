#ifndef MIPFINDER_FILE_H
#define MIPFINDER_FILE_H

#include <filesystem>
#include <fstream>

namespace mipfinder::file
{
  /* Returns a input filestream to the file provided. If the file cannot be
   * found, throws std::runtime_error */
  std::ifstream open(const std::filesystem::path& file);

  /* If file does not exist, throw an std::runtime_error with a message
  * containing the filename */
  void exists(const std::filesystem::path& filename);
}
#endif

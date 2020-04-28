#include <filesystem>
#include <iostream>

#include "easylogging++.h"
#include "mipfinder.h"

int main(int argc, char* argv[]) {

  if (argc != 2) {
    std::cerr << "Unknown arguments deteced.\n";
    std::cerr << "Usage: " << argv[0] << " [configuration_file]\n";
    return 1;
  }

  std::filesystem::path configuration_file = argv[1];
  mipfinder::Mipfinder mpf = mipfinder::Mipfinder(configuration_file);
  mpf.run();
}

#ifndef MIPFINDER_ARGPARSE_H
#define MIPFINDER_ARGPARSE_H

#include <string>
#include <unordered_map>

namespace mipfinder
{
    using Args = std::unordered_map<std::string, std::string>;
    std::vector<std::string> split(std::string str, std::string delimiter);
    Args parse_args(int argc, char* args[]);
    Args parse_args(const std::string& args);
}

#endif // !MIPFINDER_ARGPARSE_H

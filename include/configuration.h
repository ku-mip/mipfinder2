#ifndef MIPFINDER_CONFIGURATION_H
#define MIPFINDER_CONFIGURATION_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace mipfinder {

  /* Parameter class is is a wrapper around an std::unordered_map to overload
   * the operator[] */
  class Parameters {
  public:
    Parameters() = default;
    void insert(const std::string& parameter, const std::string& value);
    const std::string& operator[](const std::string& parameter) const;
    std::size_t size() const;
  private:
    std::unordered_map<std::string, std::string> params_;
  };

  /* Configuration is a class to parse a configuration file and provide an uniform access to global
   * parameters set in a file. Every configuration object can have its own configuration file
   * associated with it. This way it is possible to have multiple configurations in one program, if
   * so desired.
   *
   * The rules of a properly formed configuration file are as follows:
   * A header is enclosed in square brackets (`[]`) and creates a scope. Variable names must
   * not be duplicated in the header scope but variables may have the same names under different 
   * headers.
   * Parameter and it's value are set as `parameter=value`. Any extra spaces are ignored, e.g. it is
   * also allowed to write it as `parameter = value`. However extra spaces at the beginning of the
   * line are not allowed.
   * Optional parameters are preceded by `:optional:`. These do not have to be set in the configuration
   * file.
   * Any parameters NOT set as :optional: must have a value. If a compulsory paramater does not have
   * a value, the class throws std::invalid_argument */
  class Configuration {
  public:
    Configuration() = delete;
    Configuration(const std::filesystem::path& configuration_file);
    const Parameters& operator[](const std::string& header) const;
    std::size_t size() const;
  private:
    void parse(const std::filesystem::path& config_file);
    //A configuration file is represented by a map of maps, where the keys of the outer map are the
    //headers and the keys of the inner map are the parameters. The values of the inner map are the
    //parameters values.
    std::unordered_map<std::string, Parameters> config_;
    std::filesystem::path config_file_;
  };

}
#endif

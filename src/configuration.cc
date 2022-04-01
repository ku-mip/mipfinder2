#include <fstream>
#include <regex>

#include "configuration.h"


namespace mipfinder
{
    //struct Configuration;
    /**
     *  @brief  Extracts the parameter = value from configuration file and saves them to the Configuration object.
     */
    Configuration parse(const std::filesystem::path& configuration_file)
    {
        std::ifstream f{ configuration_file };
        if (!f.is_open()) {
            throw std::runtime_error("Could not open the configuration file \"" + configuration_file.string() + "\" for parsing");
        }

        Configuration configuration;
        Configuration::Header current_header;
        std::string line;
        while (getline(f, line)) {
            static constexpr auto comment_line_token = '#';
            if (line.empty()) { continue; }
            if (line.front() == comment_line_token) { continue; }

            static const std::regex header_regex(R"((?:\[)(\w+)(?:\]$))");
            static const std::regex parameter_regex(R"((\w+)(?:\W*=\W*)(\w+))");

            if (std::smatch matches; std::regex_search(line, matches, header_regex)) {
                current_header = matches[1];
                continue;
            }

            if (std::smatch matches; std::regex_search(line, matches, parameter_regex)) {
                Configuration::Parameter::Name name = matches[1];
                Configuration::Parameter::Value value = matches[2];
                if (!configuration.contains(current_header, name)) {
                    configuration.insert(current_header, mipfinder::Configuration::Parameter{ name, value });
                }
                else {
                    throw std::logic_error("Disallowed duplicate parameter \"" + name + "\" detected.");
                }
            }
        }
        return configuration;
    }


    Configuration::Configuration(const std::filesystem::path& configuration_file) : configuration_file(configuration_file)
    {
        parse(configuration_file);
    }


    void Configuration::insert(const Header& header, const Parameter& parameter)
    {
        parameters[header].push_back(parameter);
    };


    /**
     *  @return true if @parameter_name exists under a given @a header, false otherwise.
     */
    bool Configuration::contains(const Header& header, const Parameter::Name& parameter_name) const
    {
        if (!parameters.contains(header)) {
            return false;
        }

        for (const auto& param : parameters.at(header)) {
            if (parameter_name == param.name) {
                return true;
            }
        }
        return false;
    }


    /**
     *  @return  Value of a @parameter_name specified under the @a header section.
     *  @throw  std::out_of_range if @header or @parameter_name does not represent an existing
     *          parameter.
     */
    Configuration::Parameter::Value& Configuration::value(const Header& header, const Parameter::Name& parameter_name)
    {
        return const_cast<Configuration::Parameter::Value&>(static_cast<const Configuration&>(*this).value(header, parameter_name));
    }


    const Configuration::Parameter::Value& Configuration::value(const Header& header, const Parameter::Name& parameter_name) const
    {
        if (!parameters.contains(header)) {
            throw std::out_of_range(header + " was not detected in configuration file");
        }

        for (auto& param : parameters.at(header)) {
            if (parameter_name == param.name) {
                return param.value;
            }
        }
        throw std::out_of_range("No parameter under header \"" + header + "\" with name \"" + parameter_name + "\" found");
    };


    std::size_t Configuration::size() const noexcept
    {
        return parameters.size();
    }

}
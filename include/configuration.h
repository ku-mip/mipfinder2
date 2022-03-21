#ifndef MIPFINDER_CONFIGURATION_H
#define MIPFINDER_CONFIGURATION_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>
#include <string>

#include "easylogging++.h"

namespace mipfinder
{
    template <typename Parameter>
    struct Configuration;

    //Extracts the parameter=value from configuration file and saves them to the Configuration object.
    //uses regex to ensure that any user-formatting errors do not affect the parsing.
    template <typename Parameter>
    Configuration<Parameter> parse(const std::filesystem::path& configuration_file)
    {
        LOG(DEBUG) << "Parsing configuration file (" + configuration_file.string() + ")";
        std::ifstream f{ configuration_file };
        if (!f.is_open()) {
            throw std::runtime_error("Could not open the configuration file \"" + configuration_file.string() + "\" for parsing");
        }

        std::string line;
        //unsigned int config_line_nr = 0;
        Configuration configuration;
        while (getline(f, line)) {
            //    ++config_line_nr;

            static constexpr auto comment_line_token = '#';
            if (line.empty()) { continue; }
            if (line.front() == comment_line_token) { continue; }

            static const std::regex header_regex(R"((?:\[)(\w+)(?:\]$))");
            static const std::regex parameter_regex(R"((\w+)(?:\W*=\W*)(\w+))");

            Configuration<Parameter>::Header current_header;
            if (std::std::match_results<Parameter::name_type> matches; std::regex_search(line, matches, header_regex)) {
                current_header = matches[1];
                LOG(DEBUG) << "A valid header with value \"" + header << +"\" was found";
                continue;
            }

            if (std::smatch matches; std::regex_search(line, matches, parameter_regex)) {
                std::string parameter = matches[1];
                std::string value = matches[2];

                if (value.empty() && !is_optional) {
                    //If parameter is not set as optional, value cannot be empty! Leaving it empty is a logic
                    //error.
                    throw std::logic_error(configuration_file.string() + ":" + std::to_string(config_line_nr) + ": Parameter ("
                        + parameter + ") is not marked as optional and must have a value");
                }

                LOG(DEBUG) << "A valid parameter \"" + parameter + "\" with value \"" + value + "\" was found";
                if (added_params.find(parameter) == added_params.end()) {
                    config[header].insert(parameter, value);
                    added_params.insert(parameter);
                    LOG(DEBUG) << "Added " + parameter + " with " + value + " to " + header;
                    continue;
                }
                else {
                    throw std::logic_error(configuration_file.string() + ":" + std::to_string(config_line_nr) + ": Duplicate"
                        " parameters (" + parameter + ") are not allowed");
                }

        }


        //    //The logic is as follows:
        //    //First check for the header line. Due to the way the regexes were written, the parameter=value 
        //    //regex also matches the header line, therefore we have to check for that first.
        //    //If there is no match, check whether the line contains the word ":optional:". This signals that
        //    //the parameter is optional. The flag is set and the line is then matched against the
        //    //parameter=value regex, extracting the information.
        //    static const std::regex header_regex(R"((?:^\[)(\S+)(?:\]))");
        //    static const std::regex param_value_regex(R"((?:^[\s]+)*(?:\:optional\:)?([\S]+)(?:[=])([\S]+)*)");
        //    std::smatch matches;



        //    bool is_optional = false;
        //    if (line.find(":optional:") != std::string::npos) {
        //        LOG(DEBUG) << "Optional token was found";
        //        is_optional = true;
        //    }

        //    if (std::regex_search(line, matches, param_value_regex)) {
        //        std::string parameter = matches[1];
        //        std::string value = matches[2];

        //        if (value.empty() && !is_optional) {
        //            //If parameter is not set as optional, value cannot be empty! Leaving it empty is a logic
        //            //error.
        //            throw std::logic_error(configuration_file.string() + ":" + std::to_string(config_line_nr) + ": Parameter ("
        //                + parameter + ") is not marked as optional and must have a value");
        //        }

        //        LOG(DEBUG) << "A valid parameter \"" + parameter + "\" with value \"" + value + "\" was found";
        //        if (added_params.find(parameter) == added_params.end()) {
        //            config[header].insert(parameter, value);
        //            added_params.insert(parameter);
        //            LOG(DEBUG) << "Added " + parameter + " with " + value + " to " + header;
        //            continue;
        //        }
        //        else {
        //            throw std::logic_error(configuration_file.string() + ":" + std::to_string(config_line_nr) + ": Duplicate"
        //                " parameters (" + parameter + ") are not allowed");
        //        }
        //    }
        //}
        //LOG(DEBUG) << "Finised parsing configuration file";
        //return config;
    }
}


namespace mipfinder {

    /* Parameter class is is a wrapper around an std::unordered_map to overload
     * the operator[] */


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
    template <typename Name, typename Value>
    struct Basic_Parameter
    {
        using name_type = Name;
        using value_type = Value;
        Name name;
        Value value;
    };

    using Parameter = Basic_Parameter<std::string, std::string>;

    template <typename Parameter = Parameter>
    class Configuration
    {
    public:
        using Header = std::string;
        using ParameterName = Parameter::name_type;
        using ParameterValue = Parameter::value_type;

        Configuration() = default;
        Configuration(std::filesystem::path configuration_file)
        {
            parse<Parameter>(configuration_file);
        }

        void insert(const Header& header, const Parameter& parameter)
        {
            parameters[header].push_back(parameter);
        };

        void insert(const Header& header)
        {
            parameters[header];
        }

        bool contains(const Header& header, const ParameterName& parameter_name)
        {
            for (const auto& param : parameters[header]) {
                if (parameter_name == param.name) {
                    return true;
                }
            }
            return false;
        }

        ParameterValue& value(const Header& header, const ParameterName& parameter_name)
        {
            for (const auto& param : parameters[header]) {
                if (parameter_name == param.name) {
                    return param.value;
                }
            }
            throw std::out_of_range("No parameter under header \"" + header + "\" with name \"" + parameter_name + "\" found");
        };

        std::size_t size() const noexcept
        {
            return parameters.size();
        }

    private:
        std::unordered_map<Header, std::vector<Parameter>> parameters;
        std::filesystem::path configuration_file;
    };

}
#endif

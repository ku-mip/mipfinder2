#ifndef MIPFINDER_CONFIGURATION_H
#define MIPFINDER_CONFIGURATION_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace mipfinder
{
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
     */
    class Configuration
    {
    public:
        struct Parameter
        {
            using Name = std::string;
            using Value = std::string;
            Name name;
            Value value;
        };

        using Header = std::string;

        Configuration() = default;
        Configuration(const std::filesystem::path& configuration_file);

        void insert(const Header& header, const Parameter& parameter);

        /**
         *  @return true if @parameter_name exists under a given @a header, false otherwise.
         */
        bool contains(const Header& header, const Parameter::Name& parameter_name) const;

        /**
         *  @return  Value of a @parameter_name specified under the @a header section.
         *  @throw  std::out_of_range if @header or @parameter_name does not represent an existing
         *          parameter.
         */
        Parameter::Value& value(const Header& header, const Parameter::Name& parameter_name);

        std::size_t size() const noexcept;

    private:
        std::unordered_map<Header, std::vector<Parameter>> parameters;
        std::filesystem::path configuration_file;
    };

}
#endif

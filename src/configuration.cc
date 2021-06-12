#include <fstream>
#include <regex>
#include <unordered_set>

#include "configuration.h"
#include "easylogging++.h"

//////////////////////////
//                      //
//   Class Parameters   //
//                      //
//////////////////////////

void mipfinder::Parameters::insert(const std::string& parameter, const std::string& value)
{
	params_.emplace(parameter, value);
}

const std::string& mipfinder::Parameters::operator[](const std::string& parameter) const
{
	return params_.at(parameter);
}

/////////////////////////////
//                         //
//   Class Configuration   //
//                         //
/////////////////////////////

mipfinder::Configuration::Configuration(const std::filesystem::path& configuration_file)
	: config_file_(configuration_file)
{
	parse(config_file_);
}

const mipfinder::Parameters& mipfinder::Configuration::operator[](const std::string& header) const
{
	return config_.at(header);
}

//Extracts the parameter=value from configuration file and saves them to the Configuration object.
//uses regex to ensure that any user-formatting errors do not affect the parsing.
void mipfinder::Configuration::parse(const std::filesystem::path& configuration_file)
{
	LOG(DEBUG) << "Parsing configuration file (" + configuration_file.string() + ")";
	std::ifstream f;
	f.open(configuration_file);
	std::string line;
	std::string header;
	std::unordered_set<std::string> added_params;
	unsigned int config_line_nr = 0;
	while (getline(f, line)) {
		//LOG(DEBUG) << "Parsing \"" << line << "\"";
		++config_line_nr;
		//Empty lines are used for formatting the configuration file
		if (line.empty()) {
			LOG(DEBUG) << "Ignoring empty line";
			continue;
		}

		if (line.front() == '#') { //Lines marked with # are comments in the configuration file
			LOG(DEBUG) << "Ignoring comment line";
			continue;
		}

		//The logic is as follows:
		//First check for the header line. Due to the way the regexes were written, the parameter=value 
		//regex also matches the header line, therefore we have to check for that first.
		//If there is no match, check whether the line contains the word ":optional:". This signals that
		//the parameter is optional. The flag is set and the line is then matched against the
		//parameter=value regex, extracting the information.
		static const std::regex header_regex(R"((?:[\[])(\S+)(?:[\]]))");
		static const std::regex param_value_regex(R"((?:^[\s]+)*(?:\:optional\:)?([\S]+)(?:[=])([\S]+)*)");
		std::smatch matches;
		if (std::regex_search(line, matches, header_regex)) {
			header = matches[1];
			//If first (or new) header is found, reset the "added_params" set since parameters are allowed
			//to be duplicated within separate headers.
			added_params.clear();
			config_[header] = Parameters();
			LOG(DEBUG) << "A valid header with value \"" + header << + "\" was found";
			continue;
		}

		bool is_optional = false;
		if (line.find(":optional:") != std::string::npos) {
			LOG(DEBUG) << "Optional token was found";
			is_optional = true;
		}

		if (std::regex_search(line, matches, param_value_regex)) {
			std::string parameter = matches[1];
			std::string value = matches[2];

			if (value.empty() && !is_optional) {
				//If parameter is not set as optional, value cannot be empty! Leaving it empty is a logic
				//error.
				throw std::logic_error(config_file_.string() + ":" + std::to_string(config_line_nr) + ": Parameter ("
									   + parameter + ") is not marked as optional and must have a value");
			}

			LOG(DEBUG) << "A valid parameter \"" + parameter + "\" with value \"" + value + "\" was found";
			if (added_params.find(parameter) == added_params.end()) {
				config_[header].insert(parameter, value);
				added_params.insert(parameter);
				LOG(DEBUG) << "Added " + parameter + " with " + value + " to " + header;
				continue;
			}
			else {
				throw std::logic_error(config_file_.string() + ":" + std::to_string(config_line_nr) + ": Duplicate"
									   " parameters (" + parameter + ") are not allowed");
			}
		}
	}
	LOG(DEBUG) << "Finised parsing configuration file";
	return;
}

std::size_t mipfinder::Configuration::size() const
{
	return config_.size();
}

std::size_t mipfinder::Parameters::size() const
{
	return params_.size();
}
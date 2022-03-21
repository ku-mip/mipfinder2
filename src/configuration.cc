#include <fstream>
#include <regex>
#include <unordered_set>

#include "configuration.h"
#include "easylogging++.h"


namespace
{
	////Extracts the parameter=value from configuration file and saves them to the Configuration object.
	////uses regex to ensure that any user-formatting errors do not affect the parsing.
	//mipfinder::Configuration parse(const std::filesystem::path& configuration_file)
	//{
	//	LOG(DEBUG) << "Parsing configuration file (" + configuration_file.string() + ")";
	//	std::ifstream f;
	//	f.open(configuration_file);
	//	if (!f.is_open()) {
	//		throw std::runtime_error("Could not open " + configuration_file.string() + " for parsing");
	//	}

	//	std::string line;
	//	std::string header;
	//	std::unordered_set<std::string> added_params;
	//	unsigned int config_line_nr = 0;
	//	mipfinder::Configuration config;
	//	while (getline(f, line)) {
	//		++config_line_nr;

	//		const auto comment_line_token = '#';
	//		if (line.empty()) {	continue; }
	//		if (line.front() == comment_line_token) { continue; }

	//		//The logic is as follows:
	//		//First check for the header line. Due to the way the regexes were written, the parameter=value 
	//		//regex also matches the header line, therefore we have to check for that first.
	//		//If there is no match, check whether the line contains the word ":optional:". This signals that
	//		//the parameter is optional. The flag is set and the line is then matched against the
	//		//parameter=value regex, extracting the information.
	//		static const std::regex header_regex(R"((?:^\[)(\S+)(?:\]))");
	//		static const std::regex param_value_regex(R"((?:^[\s]+)*(?:\:optional\:)?([\S]+)(?:[=])([\S]+)*)");
	//		std::smatch matches;

	//		if (std::regex_search(line, matches, header_regex)) {
	//			header = matches[1];
	//			//If a header is found, reset the "added_params" set since parameters are allowed
	//			//to be duplicated within separate headers.
	//			added_params.clear();
	//			config[header] = mipfinder::Parameters{};
	//			LOG(DEBUG) << "A valid header with value \"" + header << +"\" was found";
	//			continue;
	//		}

	//		bool is_optional = false;
	//		if (line.find(":optional:") != std::string::npos) {
	//			LOG(DEBUG) << "Optional token was found";
	//			is_optional = true;
	//		}

	//		if (std::regex_search(line, matches, param_value_regex)) {
	//			std::string parameter = matches[1];
	//			std::string value = matches[2];

	//			if (value.empty() && !is_optional) {
	//				//If parameter is not set as optional, value cannot be empty! Leaving it empty is a logic
	//				//error.
	//				throw std::logic_error(configuration_file.string() + ":" + std::to_string(config_line_nr) + ": Parameter ("
	//					+ parameter + ") is not marked as optional and must have a value");
	//			}

	//			LOG(DEBUG) << "A valid parameter \"" + parameter + "\" with value \"" + value + "\" was found";
	//			if (added_params.find(parameter) == added_params.end()) {
	//				config[header].insert(parameter, value);
	//				added_params.insert(parameter);
	//				LOG(DEBUG) << "Added " + parameter + " with " + value + " to " + header;
	//				continue;
	//			}
	//			else {
	//				throw std::logic_error(configuration_file.string() + ":" + std::to_string(config_line_nr) + ": Duplicate"
	//					" parameters (" + parameter + ") are not allowed");
	//			}
	//		}
	//	}
	//	LOG(DEBUG) << "Finised parsing configuration file";
	//	return config;
	//}
}

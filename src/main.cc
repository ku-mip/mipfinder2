#include <cstdlib>
#include <filesystem>
#include <iostream>

#include "easylogging++.h"
#include "mipfinder.h"
#include "argparse.h"

INITIALIZE_EASYLOGGINGPP

namespace
{
	void configureLogger(bool debug_enabled)
	{
		el::Configurations default_logging;
		default_logging.setToDefault();
		if (debug_enabled) {
			default_logging.set(el::Level::Debug, el::ConfigurationType::Enabled, "true");
			default_logging.set(el::Level::Debug, el::ConfigurationType::Filename, "debug.log");
		}
		default_logging.set(el::Level::Info, el::ConfigurationType::Filename, "mipfinder.log");
		default_logging.set(el::Level::Info, el::ConfigurationType::Format, "%datetime %level %msg");
		el::Loggers::reconfigureLogger("default", default_logging);
	}

	struct CommandLineArguments
	{
		std::string configuration_file = "";
		bool debug_mode = false;
	};

	CommandLineArguments extractArguments(const mipfinder::Args& args)
	{
		CommandLineArguments cmd_args;
		if (args.contains("-d") || args.contains("-debug")) {
			cmd_args.debug_mode = true;
		}
		if (args.contains("-c")) {
			cmd_args.configuration_file = args.at("-i");
		}
		else if (args.contains("--configuration")) {
			cmd_args.configuration_file = args.at("--input");
		}
		return cmd_args;
	}
}

int main(int argc, char* argv[])
{
	auto cmd_args = mipfinder::parse_args(argc, argv);
	auto parsed_args = extractArguments(cmd_args);
	configureLogger(parsed_args.debug_mode);

	if (parsed_args.configuration_file.empty()) {
		std::cerr << "No configuration file supplied. Please use \"-c [configuration_file]\"\n";
		return EXIT_FAILURE;
	}

	mipfinder::Mipfinder mpf = mipfinder::Mipfinder(parsed_args.configuration_file);
	mpf.run();

	return EXIT_SUCCESS;
}

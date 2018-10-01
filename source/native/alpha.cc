// This is t./he main function for the NATIVE version of this project.

#include <iostream>

#include "config/command_line.h"
#include "config/ArgManager.h"

#include "../hp_config.h"
#include "../Experiment.h"



int main(int argc, char* argv[])
{
	// Read configs.
	std::string config_fname = "configs.cfg";
	auto args = emp::cl::ArgManager(argc, argv);
	HPConfig config;
	config.Read(config_fname);

	if (args.ProcessConfigOptions(config, std::cout, config_fname, "../hp_config.h") == false)
		exit(0);
	if (args.TestUnknown() == false)
		exit(0);

	Experiment e(config);
	e.Run();

	return 0;
}

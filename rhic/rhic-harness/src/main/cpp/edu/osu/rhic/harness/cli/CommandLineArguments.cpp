/*
 * CommandLineArguments.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include "edu/osu/rhic/harness/cli/CommandLineArguments.h"

namespace rhic
{

	cli_arguments::cli_arguments()
	{

		options =
		{{
			{"test",  't', "RUN_TEST", OPTION_ARG_OPTIONAL, "Run software tests"},
			{"hydro",  'h', "RUN_HYDRO", OPTION_ARG_OPTIONAL, "Run hydrodynamic simulation"},
			{"output",  'o', "OUTPUT_DIRECTORY", 0, "Path to output directory"},
			{"config", 'c', "CONFIG_DIRECTORY", 0, "Path to configuration directory"},
			{0}
		}};
			
		args_doc = "ARG1 ARG2";
			
		doc = "Run -- A program to run a single viscous hydrodynamic simulation of a relativistic heavy ion collision";
			
		argp = {options.data(), parse_opt, args_doc.c_str(), doc.c_str()};
	}

	/*
	   PARSER. Field 2 in ARGP.
	   Order of parameters: KEY, ARG, STATE.
		 Now uses internal cli_args state, keeps argp_state for compat
	*/
	error_t cli_arguments::parse_opt(int key, char *arg, struct argp_state *state)
	{
		switch (key) {
		case 't':
			cli_args.run_test = true;
			break;
		case 'h':
			cli_args.run_hydro = true;
			break;
		case 'o':
			cli_args.output_dir = std::string(arg);
			break;
		case 'c':
			cli_args.config_dir = std::string(arg);
			break;
	//	case ARGP_KEY_ARG:
	//		if (state->arg_num >= 2) {
	//			argp_usage(state);
	//		}
	//		arguments->args[state->arg_num] = arg;
	//		break;
	//	case ARGP_KEY_END:
	//		if (state->arg_num < 2) {
	//			argp_usage(state);
	//		}
	//		break;
		default:
			return ARGP_ERR_UNKNOWN;
		}
		return 0;
	}

	error_t cli_arguments::load_cli_args(int argc, char **argv)
	{
		argp_parse(&argp, argc, argv, 0, 0, &cli_args); 
		return 0;
	}

	error_t cli_arguments::print_arg_vals()
	{
		std::cout <<   "config_dir = " << cli_args.config_dir
		          << "\noutput_dir = " << cli_args.output_dir
		          << "\nrun_hydro  = " << cli_args.run_hydro
		          << "\nrun_test   = " << cli_args.run_test << std::endl;
		return 0;
	}

}

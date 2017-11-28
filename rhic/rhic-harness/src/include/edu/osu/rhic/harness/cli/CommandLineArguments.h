/*
 * CommandLineArguments.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef COMMANDLINEARGUMENTS_H_
#define COMMANDLINEARGUMENTS_H_

#include <stdbool.h>
#include <argp.h>

namespace rhic
{
	class cli_arguments
	{
		public:
			cli_arguments();
			/*
			   PARSER. Field 2 in ARGP.
			   Order of parameters: KEY, ARG, STATE.
			*/
			error_t parse_opt(int key, char *arg, struct argp_state *state);
			/*
				 The main function.
				 Notice how now the only function call needed to process
				 all command-line options and arguments nicely
				 is argp_parse.
			*/
			error_t load_cli_args(int argc, char **argv);

		private:
			struct cli_args
			{
				std::array<std::string, 2> args;            /* ARG1 and ARG2 */
				bool run_test = false;
				bool run_hydro = false;
				std::string config_dir;              /* The -v flag */
				std::string output_dir;            /* Argument for -o */
			};

			std::string argp_program_ver;
			std::string argp_program_bug_addr;

			/*
				 OPTIONS.  Field 1 in ARGP.
				 Order of fields: {NAME, KEY, ARG, FLAGS, DOC}.
			*/
			std::array<struct argp_option, 5> options;

			/*
				 ARGS_DOC. Field 3 in ARGP.
				 A description of the non-option command-line arguments
					 that we accept.
			*/
			std::string args_doc;

			/*
				DOC.  Field 4 in ARGP.
				Program documentation.
			*/
			std::string doc;

			/*
				 The ARGP structure itself.
			*/
			struct argp argp; 

			struct cli_args cli_args;
	}
}

#endif /* COMMANDLINEARGUMENTS_H_ */

/*
  ReadCommandLine.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 9/2006

  Reads command line input for Concentrations.c.  The argument is the
  either the prefix for the files that contain the input data, e.g.,
  prefix.cx and prefix.con, or it is the full file name of the input,
  e.g., prefix.ox-epairs, in which case the concentration information
  is if prefix.con.  Uses the package getopt.h to retrieve option
  flags.

  Justin Bois, Caltech, 2 September 2006
*/

#include "ReadCommandLine.h"
#include "constants.h"

/* ******************************************************************************** */
void ReadCommandLine(int nargs, char **args, char *cxFile, char *conFile, 
		     char *logFile, char *eqFile,
		     int *SortOutput, int *MaxIters, double *tol, double *kT,
		     int *MaxNoStep, int *MaxTrial, double *PerturbScale, int *quiet,
		     int *WriteLogFile, int *Toverride, int *NoPermID, 
		     unsigned long *seed, int * NUPACK_VALIDATE) {

  int options;  // Counters used in getting flags
  int ShowHelp; // ShowHelp = 1 if help option flag is selected
  char prefix[MAXLINE]; // The prefix for the input and output files
  char InputStr[MAXLINE]; // Dummy string for storing flag options from command line
  FILE *fp; // The cx file, used to check if we can open it.

  if (nargs == 1) {
    printf("For instructions on running this program, run it with the ");
    printf("-help flag.\n\nExiting....\n\n");
    exit(ERR_NOINPUT);
  }

  *SortOutput = 1;  // Default is to sort the output by concentration
  *MaxIters = 10000; // Default maxiters
  *tol = 0.0000001; // Default tolerance is 0.00001% of minimum of the
                    // minimum count among the single-strands
  *kT = kB*(37.0 + ZERO_C_IN_KELVIN); // Default temperature is 37 deg. C
  *quiet = 0; // Default is to show messages on the screen
  *MaxNoStep = 50; // Default is 50 iterations with no step
  *MaxTrial = 100000; // Default is maximum of 100,000 trials
  *PerturbScale = 100; // Default is a scale of 100 on the perturb scale
  *WriteLogFile = 0; // Default is not to write a log file
  *Toverride = 0; // Default is to either use T = 37 or that specified in input file
  *NoPermID = 0; // Default is to use .ocx file => permutation IDs in file
  *seed = 0; // Default is to seed off the clock.
  double cutoff = 0.001; // Default cutoff
  *NUPACK_VALIDATE = 0;
  ShowHelp = 0;

  SetExecutionPath(nargs, args);
  
  // Get the option flags
  while (1)
    {
      static struct option long_options [] =
	      {
	        {"T",             required_argument,  0, 'd'},
	        {"quiet",         no_argument,        0, 'e'},
	        {"help",          no_argument,        0, 'h'},
	        {"perturbscale",  required_argument,  0, 'i'},
	        {"writelogfile",  no_argument,        0, 'j'},
	        {"seed",          required_argument,  0, 'm'},
          {"validate",      no_argument,        0, 'o'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;


      options = getopt_long_only (nargs, args, 
				  "d:ehi:jm:o", long_options, 
				  &option_index);

      // Detect the end of the options.
      if (options == -1)
        break;

      switch (options)
        {
	      case 'd':
	        strcpy(InputStr,optarg);
	        *kT = kB*(str2double(InputStr) + ZERO_C_IN_KELVIN);
	        *Toverride = 1;
	        break;

	      case 'e':
	        *quiet = 1;
	        break;

        case 'h':
	        ShowHelp = 1;
          break;

        case 'i':
	        strcpy(InputStr,optarg);
	        (*PerturbScale) = str2double(InputStr);
          break;

	      case 'j':
	        *WriteLogFile = 1;
	        break;

	      case 'm':
	        strcpy(InputStr,optarg);
	        (*seed) = (unsigned long)atol(InputStr);
	        break;

        case 'o':
          *NUPACK_VALIDATE = 1;
          *tol = 0.0000000000001;
	        *SortOutput = 3;
	        // NoSortOutputOption = 0; // Record that we've selected a sorting option
          cutoff = 0;
          break;

        case '?':
          // getopt_long already printed an error message.
          break;
        
        default:
          abort();
        }
    }

  if (ShowHelp) {
    DisplayHelpConc();
    exit(ERR_HELP);
  }

  if (*SortOutput > 4) {
    if (*quiet == 0) {
      printf("Sorting option is an integer from 0 to 4.  Output will be sorted by\n");
      printf("concentration.\n\n");
    }
    *SortOutput = 1;
  }

  // Get the the input file
  if (optind == nargs) { // There's no input from the user
    if (*quiet == 0) {
      printf("You must have a prefix or an input file on the command line.\n");
      printf("For instructions on running this program, run it with the ");
      printf("-help flag.\n\nExiting....\n\n");
    }
    exit(ERR_NOINPUT);
  }
  else {
    strcpy(prefix,args[optind]);
  }

  // Name the files
  strcpy(cxFile,prefix);
  strcpy(conFile,prefix);
  strcpy(logFile,prefix);
  strcpy(eqFile,prefix);
  strcat(cxFile,".ocx");
  strcat(conFile,".con");
  strcat(logFile,".log");
  strcat(eqFile,".eq");

  // Do a quick check to make sure the cx file exists before we proceed
  if ((fp = fopen(cxFile,"r")) == NULL) {
    if (*quiet == 0) {
      printf("Error opening %s!\n\nExiting....\n",cxFile);
    }
    exit(ERR_CX);
  }
  fclose(fp);

} 
/* ******************************************************************************** */


/* ******************************************************************************** */
void DisplayHelpConc() {
  /*
    Displays the contents of the Concentrations help file.
  */

  printf("Please read the NUPACK User Guide for detailed instructions.\n");
  printf("Usage: concentrations [OPTIONS] PREFIX\n");
  printf("Calculate concentrations for each complex specified\n");
  printf("Options:\n");
  printf(" -quiet               suppress output to the screen\n");
  printf("\n");

}


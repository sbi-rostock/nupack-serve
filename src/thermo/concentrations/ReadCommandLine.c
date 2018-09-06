/*
 * ReadCommandLine.c is part of the NUPACK software suite
 * Copyright (c) 2007 Caltech. All rights reserved.
 * Coded by: Justin Bois 9/2006
 */

#include "ReadCommandLine.h"
#include "constants.h"


void ReadCommandLine(int nargs, char **args, int *SortOutput, int *MaxIters,
        double *tol, double *kT, int *MaxNoStep, int *MaxTrial,
        double *PerturbScale, int *quiet, int *Toverride, int *NoPermID,
        unsigned long *seed, int *NUPACK_VALIDATE){

  int options;
  int ShowHelp; // 1 if help option flag is selected
  char prefix[MAXLINE]; // prefix for the input and output files
  char InputStr[MAXLINE]; // dummy string for storing flag options from command line

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
  *quiet = 1; // Default is not to show messages on the screen
  *MaxNoStep = 50; // Default is 50 iterations with no step
  *MaxTrial = 100000; // Default is maximum of 100,000 trials
  *PerturbScale = 100; // Default is a scale of 100 on the perturb scale
  *Toverride = 0; // Default is to either use T = 37 or that specified in input file
  *NoPermID = 0; // Default is to use .ocx file => permutation IDs in file
  *seed = 0; // Default is to seed off the clock.
  double cutoff = 0.001; // Default cutoff
  *NUPACK_VALIDATE = 0;
  ShowHelp = 0;

  SetExecutionPath(nargs, args);

  // Get the option flags
  while (1){
    static struct option long_options[] = {
        {"T",             required_argument,  0, 'd'},
        {"help",          no_argument,        0, 'h'},
        {"validate",      no_argument,        0, 'o'},
        {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;


    options = getopt_long_only (nargs, args, "d:ho", long_options,
        &option_index);

    // detect the end of the options
    if (options == -1){
      break;
    }

    switch(options){
      case 'd':
        strcpy(InputStr,optarg);
        *kT = kB*(str2double(InputStr) + ZERO_C_IN_KELVIN);
        *Toverride = 1;
        break;

      case 'h':
        ShowHelp = 1;
        break;

      case 'o':
        *NUPACK_VALIDATE = 1;
        *tol = 0.0000000000001;
        *SortOutput = 3;
        cutoff = 0;
        break;

      default:
        abort();
    }
  }

  if(ShowHelp){
    DisplayHelpConc();
    exit(ERR_HELP);
  }

  if(*SortOutput > 4){
    *SortOutput = 1;
  }

  // Get the the input file
  if(optind == nargs){ // no input from the user
    exit(ERR_NOINPUT);
  }
  else{
    strcpy(prefix,args[optind]);
  }
}



/* Displays the contents of the Concentrations help file
 */
void DisplayHelpConc(){

  printf("Please read the NUPACK User Guide for detailed instructions.\n");
  printf("Usage: concentrations [OPTIONS] PREFIX\n");

}


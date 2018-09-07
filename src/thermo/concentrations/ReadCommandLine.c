/*
 * ReadCommandLine.c is part of the NUPACK software suite
 * Copyright (c) 2007 Caltech. All rights reserved.
 * Coded by: Justin Bois 9/2006
 */

#include "ReadCommandLine.h"
#include "constants.h"


void ReadCommandLine(int nargs, char** args, int* SortOutput, int* MaxIters,
        double* tol, double* kT, int* MaxNoStep, int* MaxTrial,
        double* PerturbScale, int* Toverride, unsigned long* seed,
        double* cutoff, int* NUPACK_VALIDATE){

  int options;
  int option_index = 0;
  char InputStr[MAXLINE];


  // default settings
  *NUPACK_VALIDATE = 0;
  *kT = kB*(37.0 + ZERO_C_IN_KELVIN); // temperature
  *PerturbScale = 100; // perturbation scale
  *MaxIters = 10000;  // max iterations
  *MaxTrial = 100000; // maximum number of trials
  *seed = 0; // seed off the clock
  *Toverride = 0;   // do not override default temperature value
  *cutoff = 0.001;  // cutoff
  *MaxNoStep = 50;  // iterations
  *SortOutput = 1;  // sort output by concentration
  *tol = 0.0000001; // tolerance percentage of the minimum count among
                    // single-strands

  SetExecutionPath(nargs, args);

  // get the option flags
  while (1){
    static struct option long_options[] = {
        {"T",             required_argument,  0, 'd'},
        {"validate",      no_argument,        0, 'o'},
        {0, 0, 0, 0}
    };


    options = getopt_long_only (nargs, args, "d:o", long_options,
        &option_index);

    // detect the end of the options
    if (options == -1){
      break;
    }

    switch(options){
      case 'd':
        strcpy(InputStr, optarg);
        *kT = kB * (str2double(InputStr) + ZERO_C_IN_KELVIN);
        *Toverride = 1;
        break;

      case 'o':
        *NUPACK_VALIDATE = 1;
        *tol = 0.0000000000001;
        *SortOutput = 3;
        *cutoff = 0;
        break;

      default:
        abort();
    }
  }

  if(*SortOutput > 4){
    *SortOutput = 1;
  }
}


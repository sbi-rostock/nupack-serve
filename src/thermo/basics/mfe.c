/*
    mfe.c is part of the NUPACK software suite
    Copyright (c) 2007 Caltech. All rights reserved.
    Coded by: Robert Dirks, 6/2006 and Justin Bois 1/2007


    This function will calculate and print all mfe structures (if the
    -degenerate flag is selected) or one mfe structure (if not),
    taking into account symmetry corrections.  Consequently, if
    -degenerate is chosen, this could scale as poorly as exponential
    with regard to space and time.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <thermo/core.h>


int main(int argc, char *argv[]) {

  char seq[MAXSEQLENGTH];
  int seqNum[MAXSEQLENGTH+1];
  int isNicked[MAXSEQLENGTH];
  int nNicks = 0;


  int nicks[MAXSTRANDS];
  int nickIndex;
  int **etaN;
  int pf_ij;

  int complexity = 3;
  int length;
  int tmpLength;
  DBL_TYPE mfe;
  int vs;
  char inputFile[MAXLINE];
  int inputFileSpecified;

  // provenance blocks
  int len_header = 1000;
  int len_parameters = 1000;
  int len_structures = 1000;
  int len_provenance;

  dnaStructures mfeStructs = {NULL, 0, 0, 0, NAD_INFINITY};


  strcpy(inputFile, "");

  inputFileSpecified = ReadCommandLineNPK( argc, argv, inputFile);
  if(NupackShowHelp){
    printf("Usage: mfe [OPTIONS] PREFIX\n");
    printf("Compute and store the minimum free energy and the MFE\n");
    printf("secondary structure(s) of the input sequence.\n");
    printf("Example: mfe -multi -T 25 -material dna example\n");
    PrintNupackThermoHelp();
    PrintNupackUtilitiesHelp();
    exit(1);
  }

  getUserInput(seq, &vs, NULL, NULL);


  /* echo provenance header
   */

  // allocate provenance block
  char *header = malloc(sizeof(char) * len_header);
  if(!header){
    exit(1);
  }
  for(int y=0 ; y<len_header ; ++y){
    header[y] = 0;
  }

  // fill provenance block
  len_provenance = header2provenance(header, argc, argv);
  for(int y=0 ; y<len_provenance ; ++y){
    printf("%c", header[y]);
  }

  // free provenance block
  free(header);
  header = NULL;


  /* echo provenance setting informations
   */

  // allocate provenance block
  char *parameters = malloc(sizeof(char) * len_parameters);
  if(!parameters){
    exit(1);
  }
  for(int y=0 ; y<len_parameters ; ++y){
    parameters[y] = 0;
  }

  // fill provenance block
  len_provenance = parameters2provenance(parameters, argc, argv, seq,
    NULL, NULL);
  for(int y=0 ; y<len_provenance ; ++y){
    printf("%c", parameters[y]);
  }

  // free provenance block
  free(parameters);
  parameters = NULL;


  if(!DO_PSEUDOKNOTS){
    complexity = 3;
  }
  else{
    complexity = 5;
  }

  tmpLength = strlen(seq);
  convertSeq(seq, seqNum, tmpLength);

  mfe = mfeFullWithSym(seqNum, tmpLength, &mfeStructs, complexity,
    DNARNACOUNT, DANGLETYPE, TEMP_K - ZERO_C_IN_KELVIN, vs, ONLY_ONE_MFE,
    SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);


  //the rest is for printing purposes
  tmpLength = length = strlen(seq);


  for(int i=0 ; i<tmpLength ; ++i){
    isNicked[i] = 0;
    if(seq[i] == '+') {
      --length;
      isNicked[i - (++nNicks) - 1] = 1;
    }
  }

  //initialize nicks
  for(int i=0 ; i<MAXSTRANDS ; ++i){
    nicks[i] = -1;
  }

  nickIndex = 0;
  for(int i=0 ; i<length ; ++i){
    if(isNicked[i]){
      nicks[++nickIndex] = i;
    }
  }


  //overkill, but convenient
  etaN = (int**) malloc(sizeof(int*) * (length*(length+1)/2 + (length+1)));
  InitEtaN(etaN, nicks, length);


  /* echo provenance DNA structures
   */

  // allocate provenance blocks
  char *structures = malloc(sizeof(char) * len_structures);
  if(!structures){
    exit(1);
  }
  for(int y=0 ; y<len_structures ; ++y){
    structures[y] = 0;
  }

  // fill provenance block
  len_provenance = dnastructures2provenance(structures, &mfeStructs, etaN,
    nicks, vs);
  for(int y=0 ; y<len_provenance ; ++y){
    printf("%c", structures[y]);
  }

  // free provenance block
  free(structures);
  structures = NULL;


  clearDnaStructures(&mfeStructs);

  for(int i=0 ; i<=(length-1) ; ++i){
    for(int j=(i-1) ; j<=(length-1); ++j){
      pf_ij = pf_index(i, j, length);
      free(etaN[pf_ij]);
    }
  }
  free(etaN);

  return 0;
}


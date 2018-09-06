/*
 * concentrations.c is part of the NUPACK software suite
 * Copyright (c) 2007 Caltech. All rights reserved.
 * Coded by: Justin Bois 9/2006
 *
 * This program does thermodynamic analysis of interacting nucleic acid strands
 * in a box where the partition functions for the complexes are known.
 * The algorithm is described in Dirks, Bois, Schaeffer, Winfree, and Pierce,
 * "Thermodynamic Analysis of interacting nucleic acid strands", SIAM Review.
 * Variable names in the code should be referenced with variable names in that
 * paper.
 * The program calculates the equilibrium concentrations of possible complexes
 * in the thermodynamic limit, i.e., for a large system with many strands. The
 * trust region algorithm for solving the dual problem is that in Nocedal and
 * Wright, Numerical Optimization, 1999, page 68, with the dogleg method on
 * page 71. The subroutine used to do this calculation is CalcConc.c.
 *
 * For usage instructions, input and output formats, etc., see the associated
 * manual.
 */

#include "constants.h"
#include "CalcConc.h"
#include "FracPair.h"
#include "InputFileReader.h"
#include "OutputWriter.h"
#include "ReadCommandLine.h"


int main(int argc, char *argv[]) {

  unsigned long seed; // seed for random number generation
  int numSS;    // number of single-strand (monomer) types
  int numSS0;   // number of monomer types including those with zero concentration
  int numTotal; // total number of complexes
  int nTotal;   // total number of permutations
  int LargestCompID; // largest complex ID
  int MaxIters;   // maximum number of iterations in trust region method
  int SortOutput; // sorting options for output
  int NoPermID;   // 1 if there are no perumtation IDs in 2nd column of input file
  int CalcConcConverge; // 1 for convergence, 0 otherwise
  double tol;      // absolute tolerance is tol*(mininium monomer init. conc.)
  double deltaBar; // maximum allowed step size in trust region method
  double eta;      // eta parameter in trust region method, 0 < eta < 0.25
  double kT;       // thermal energy in kcal/mol
  double temperature;
  double MolesWaterPerLiter; // moles of water per liter
  int Toverride; // 1 if the user has enforced a temperature in the command line
  int MaxNoStep; // maximum number of iterations allowed without taking a step
  int MaxTrial;  // maximum number ot perturbations allowed in a calculation
  double PerturbScale; // multiplier on the random number for perturbations
  int **A;    // number of monomers of type i in complex j
  double *G;  // free energies of complexes
  double *x;  // the mole fractions
  double *x0; // total concentrations of single-species
  double *conc;
  int NUPACK_VALIDATE; // 1 if validation mode (14 digit printout)
  int *numPermsArray; // number of permutations of each species
  int *CompIDArray;   // complex IDs
  int *PermIDArray;   // permutation IDs

  // provenance blocks
  int len_header = 1000;
  int len_parameters = 1000;
  int len_concentrations = 2000;
  int len_provenance;


  eta = TRUST_REGION_ETA;
  deltaBar = TRUST_REGION_DELTABAR;


  // read command line arguments
  ReadCommandLine(argc, argv, &SortOutput, &MaxIters, &tol, &kT, &MaxNoStep,
        &MaxTrial, &PerturbScale, &Toverride, &NoPermID, &seed,
        &NUPACK_VALIDATE);


  /* echo provenance header starts
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
  len_provenance = concentrations_header(header, argc, argv);
  for(int y=0 ; y<len_provenance ; ++y){
    printf("%c", header[y]);
  }

  // free provenance block
  free(header);
  header = NULL;
  /*
   * echo provenance header ends */


  // get the system's size
  getSize(&numSS,&numTotal,&nTotal,&LargestCompID,&numPermsArray);


  // store input parameters
  struct InStruct* InputStruct = malloc(sizeof(InStruct) * nTotal);
  for(int j=0 ; j<nTotal; ++j){
    InputStruct[j].Aj = malloc (sizeof(int) * numSS);
  }

  // read input files
  MolesWaterPerLiter = ReadInputFiles(&A, &G, &CompIDArray, &PermIDArray, &x0,
        &conc, &numSS, &numSS0, &numTotal, numPermsArray, &kT, &temperature,
        Toverride, InputStruct);


  /* echo provenance parameters starts
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
  len_provenance = concentrations_parameters(parameters, numSS, conc,
        temperature);
  for(int y=0 ; y<len_provenance ; ++y){
    printf("%c", parameters[y]);
  }

  // free provenance block
  free(parameters);
  parameters = NULL;
  /*
   * echo provenance parameters ends */


  // compute convergence
  x = malloc (sizeof(double) * numTotal);
  CalcConcConverge = CalcConc(x, A, G, x0, numSS, numTotal, MaxIters, tol,
        deltaBar, eta, kT, MaxNoStep, MaxTrial, PerturbScale,
        MolesWaterPerLiter, seed);


  WriteOutput(x, G, CompIDArray, LargestCompID, numSS0, numTotal, nTotal, kT,
        SortOutput, MolesWaterPerLiter, NoPermID, NUPACK_VALIDATE, InputStruct);


  /* echo provenance concentrations starts
   */
  // allocate provenance block
  char *concentrations = malloc(sizeof(char) * len_concentrations);
  if(!concentrations){
    exit(1);
  }
  for(int y=0 ; y<len_concentrations ; ++y){
    concentrations[y] = 0;
  }

  // fill provenance block
  len_provenance = concentrations_results(concentrations, numSS, nTotal, kT,
        MolesWaterPerLiter, NoPermID, NUPACK_VALIDATE, InputStruct);

  for(int y=0 ; y<len_provenance ; ++y){
    printf("%c", concentrations[y]);
  }

  // free provenance block
  free(concentrations);
  concentrations = NULL;
  /*
   * echo provenance concentrations ends */


  // free memory allocations
  for(int i=0 ; i<numSS ; ++i){
    free(A[i]); // allocated in ReadInput
  }
  free(A); // allocated in ReadInput
  free(G); // allocated in ReadInput
  free(numPermsArray); // allocated in getSize
  free(CompIDArray);   // allocated in ReadInput
  free(PermIDArray);   // allocated in ReadInput
  free(x0); // allocated in ReadInput
  free(conc);
  free(x);

  for(int i=0 ; i<nTotal ; ++i){
    free(InputStruct[i].Aj);
  }
  free(InputStruct);


  // If didn't converge, give error message
  if (CalcConcConverge == 0) {
    exit(ERR_NOCONVERGE);
  }

  return 0;
}


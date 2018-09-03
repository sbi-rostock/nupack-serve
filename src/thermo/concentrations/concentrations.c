/*
  concentrations.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 9/2006

  This program does thermodynamic analysis of interacting nucleic acid
  strands in a box where the partition functions for the complexes are
  known.  The algorithm is described in Dirks, Bois, Schaeffer,
  Winfree, and Pierce, "Thermodynamic Analysis of interacting nucleic
  acid strands", SIAM Review, in press.  Variable names in the code
  should be referenced with variable names in that paper.

  The program calculates the equilibrium concentrations of possible
  complexes in the thermodynamic limit, i.e., for a large system with
  many strands.  The trust region algorithm for solving the dual
  problem is that in Nocedal and Wright, Numerical Optimization, 1999,
  page 68, with the dogleg method on page 71.  The subroutine used to
  do this calculation is CalcConc.c.
  
  For usage instructions, input and output formats, etc., see the
  associated manual.
  
  Justin Bois, Caltech, 2 September 2006
*/

#include "CalcConc.h"
#include "FracPair.h"
#include "InputFileReader.h"
#include "OutputWriter.h"
#include "ReadCommandLine.h"
#include "constants.h"

int main(int argc, char *argv[]) {

  int i; // Counter
  unsigned long seed; // Seed for random number generation
  char cxFile[MAXLINE]; // File containing complex ID's and free energies
  char conFile[MAXLINE]; // File containing initiatial monomer concentrations
  char logFile[MAXLINE]; // File containing data about the calculation
  char eqFile[MAXLINE];  // Name of file for equilibrium concentrations
  char fpairsFile[MAXLINE];  // Name of file for fraction base pair concentrations
  int numSS; // Number of single-strand (monomer) types. 
  int numSS0; // Number of monomer types including those with zero concentration
  int numTotal; // Total number of complexes
  int nTotal; // Total number of permutations
  int LargestCompID; // Largest complex ID
  int MaxIters; // Maximum number of iterations in trust region method
  int SortOutput; // Sorting options for output
  int quiet; // = 1 for no displays of results to screen
  int NoPermID; // = 1 if there are no perumtation IDs in 2nd column of input file
  int DoBPfracs; // = 1 if we need to calculation fraction of strands that for pairs
  int WriteLogFile; // Whether or not to write a log file
  int CalcConcConverge; // 1 is CalcConc converged and 0 otherwise
  double tol; // The absolute tolerance is tol*(mininium monomer init. conc.)
  double deltaBar; // Maximum allowed step size in trust region method
  double eta; // eta parameter in trust region method, 0 < eta < 0.25
  double kT; // The thermal energy in kcal/mol.
  double MolesWaterPerLiter; // Moles of water per liter
  int Toverride; // = 1 if the user has enforced a temperature in the command line
  int MaxNoStep; // The maximum number of iterations allowed without taking a step
  int MaxTrial; //  The maximum number ot perturbations allowed in a calculation
  double PerturbScale; // The multiplier on the random number for perturbations
  int **A; // A[i][j] is the number of monomers of type i in complex j
  double *G; // Free energies of complexes
  double *x; // The mole fractions
  double *x0; // Total concentrations of single-species
  int NUPACK_VALIDATE; // 1 if validation mode (14 digit printout)
  int *numPermsArray; // Number of permutations of each species
  int *CompIDArray; // The complex ID's
  int *PermIDArray; // Permutation ID's
  FILE *fpeq; // The .eq file, which contains the output of the file.

  // Read command line arguments
  ReadCommandLine(argc,argv,cxFile,conFile,eqFile,
		  &SortOutput,&MaxIters,&tol,&kT,&MaxNoStep,&MaxTrial,&PerturbScale,
		  &quiet,&WriteLogFile,&Toverride,&NoPermID,&seed,&NUPACK_VALIDATE);


  // Pull eta and deltaBar from global variables
  eta = TRUST_REGION_ETA;
  deltaBar = TRUST_REGION_DELTABAR;

  // Write information to .eq file
  if ((fpeq = fopen(eqFile,"w")) == NULL) {
    exit(ERR_EQ);
  }
  fprintf(fpeq,"%% NUPACK %s\n", CMAKE_NUPACK_VERSION);
  fprintf(fpeq,"%% This is %s, an output file generated for a \"concentrations\"\n",
	  eqFile);
  fprintf(fpeq,"%% calculation of equilibrium concentrations.\n");
  fprintf(fpeq,"%% For information on contents, see NUPACK manual.\n");
  fprintf(fpeq,"%% Program: concentrations\n");
  fprintf(fpeq,"%% Command: ");
  for (i = 0; i < argc; i++) {
    fprintf(fpeq,"%s ",argv[i]);
  }
  fprintf(fpeq,"\n");
  fprintf(fpeq,"%% Initial monomer concentrations:\n");
  fclose(fpeq);
  
  // Get the size of the system.
  getSize(&numSS,&numTotal,&nTotal,&LargestCompID,&numPermsArray,cxFile,
	  quiet);
  
  // Read input files and sort if necessary.
  // Note: A, G, and either x0 or m0 are all allocated in ReadInput
  MolesWaterPerLiter = ReadInputFiles(&A,&G,&CompIDArray,&PermIDArray,&x0,&numSS,
				      &numSS0,&numTotal,numPermsArray,cxFile,conFile,
				      &kT,Toverride,logFile,eqFile,
				      fpairsFile,quiet,WriteLogFile,DoBPfracs,
				      NoPermID);

  // Allocate memory for mole fractions
  x = (double *) malloc (numTotal * sizeof(double));
  
  CalcConcConverge = CalcConc(x,A,G,x0,numSS,numTotal,MaxIters,tol,deltaBar,eta,kT,
			      MaxNoStep,MaxTrial,PerturbScale,quiet,WriteLogFile,
			      logFile,MolesWaterPerLiter,seed);

  // Show warning in eq file if we failed to converge
  if (CalcConcConverge == 0) {
    if ((fpeq = fopen(eqFile,"a")) == NULL) {
      exit(ERR_EQ);
    }
    fprintf(fpeq,"%%\n");
    fprintf(fpeq,"%% ***************************************************************\n");
    fprintf(fpeq,"%%      TRUST REGION DID NOT CONVERGE DUE TO PRECISION ISSUES     \n");
    fprintf(fpeq,"%% ***************************************************************\n");
    fprintf(fpeq,"%%\n");
    fclose(fpeq);
  }

  WriteOutput(x,G,CompIDArray,LargestCompID,numSS0,numTotal,nTotal,kT,cxFile,
	      SortOutput,eqFile,MolesWaterPerLiter,quiet,NoPermID,NUPACK_VALIDATE);

  // Free memory
  for (i = 0; i < numSS; i++) {
    free(A[i]); // Allocated in ReadInput
  }
  free(A); // Allocated in ReadInput
  free(G); // Allocated in ReadInput
  free(numPermsArray); // Allocated in getSize
  free(CompIDArray); // Allocated in ReadInput
  free(PermIDArray); // Allocated in ReadInput
  free(x0); // Allocated in ReadInput
  free(x);  // Allocated in main

  // If didn't converge, give error message
  if (CalcConcConverge == 0) {
    exit(ERR_NOCONVERGE);
  }

  return 0; // Return

}

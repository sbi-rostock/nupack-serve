/*
  InputFileReader.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 9/2006

  For use with Concentrations.c. 

  Reads the input from the .cx and .con files.

  After the comment lines, which must begin with a % character, each
  row of the cxFile contains the complex ID and its free energy in
  kcal/mol.  E.g., if permutation 2 of complex ID number 17 is ABBD
  and has a free energy -18.64 kcal/mol, and there are four monomer
  types, the corresponding row in the cxFile would be:
      17  2  1  2  0  1  -18.62

  If the input NoPermID == 1, i.e., the second column above doesn't
  exist, the entry in the file looks like:
      17  1  2  0  1  -18.62

  The conFile has the initial concentrations of each of the monomer
  types in solution.  These are in units of MOLAR.

  For further info about input and output formatting, see associated
  manual.

  WARNING: This program does little format checking, so if there is an
  error in the input, it is likely to result in some strange error
  and/or seg fault.
*/


#include "InputFileReader.h" // File with important definitions
#include <float.h>
#include "constants.h"


// Structures for storing input and subsequent sorting
struct CompStruct { // Struct for complexes (used for output)
  int *Aj; // Array representing column j of A
  int numSS; // number of entries in Aj
  int CompID;
  double FreeEnergy; // Partition function for species
  char *AuxStr; // String containing auxillary information from input file
};



void getSize(int *numSS, int *numTotal, int *nTotal, int *LargestCompID,
             int **numPermsArray) {

  char newline;


  /* read number of complex IDs
   */
  printf("Enter number of complex IDs: ");
  scanf("%d%c", numTotal, &newline);
  *nTotal = *numTotal;
  *LargestCompID = *numTotal;


  /* read number of concentrations
   */
  printf("Enter number of different concentrations: ");
  scanf("%d%c", numSS, &newline);


  // Allocate memory and set array of number of permutations
  *numPermsArray = malloc(sizeof(int) * (*numTotal));
  for(int x=0 ; x<(*numTotal) ; ++x){
    (*numPermsArray)[x] = 1;
  }

  // Check to make sure the complex ID's are sequential
  for (int j=0; j<(*numTotal); ++j) {
    if ((*numPermsArray)[j] == 0) {
      exit(ERR_NONSEQUENTIAL);
    }
  }
}



/* ******************************************************************************** */
double ReadInputFiles(int ***A, double **G, int **CompIDArray,
        int **PermIDArray, double **x0, double** concentrations, int *numSS,
        int *numSS0, int *numTotal, int *numPermsArray, double *kT,
        double* temperature, int Toverride, int quiet){
  /*
    If one of the entries in the con file is zero, the problem is
    reformulated as if that strand does not exist.

    The input is stored in the arrays A, G, and x0.  A[i][j] is the
    number of monomers of type i in complex j.  G[j] is the free
    energy of complex j IN UNITS OF kT.  x0[i] is the initial mole
    fraction of monomer species i.

    The arrays CompIDArray and PermIDArray store the corresponding
    complex IDs and Permutation IDs for the entries loaded in the A
    and G.

    THE MEMORY FOR THESE ARRAYS IS ALLOCATED IN THIS FUNCTION AND MUST
    BE FREED OUTSIDE OF IT.
  */

  int i,j,k; // Counters
  struct CompStruct *InputStruct; // Struct we store the input in.
  char *tok; // Token
  char tokseps[] = " \t\n"; // Token separators
  int nSS; // Local number of single species
  int cTotal; // Local number of complexes
  int *nonzerox0; // Identities of strands that are not zero in ccon
  int *zerox0; // The identities of strands that are set to zero in ccon
  int **newA; // The matrix A for the reformulated problem with zero ccon's taken out
  int *newCompIDArray; // The comp ID's for the reformulated problem
  int *newPermIDArray; // The perm ID's for the reformulated problem
  double *newG; // Free energies for reformulated prob. with zero ccon's taken out
  double *newx0; // Mole fractions of single species with zero ccon's taken out
  long double *Q; // Partition functions for complexes
  double Gperm; // free energy of a given permutation
  int newnumTotal; // New number of complexes after zero ccon's are removed
  int newnumSS; // New number of single strands after zero ccon's are removed
  int notOK; // Whether or not an entry in A can be kept if there are zero ccon's
  int noPerms; // noPerms = 1 if permutations are not explicitly considered
  double MolesWaterPerLiter; // Moles of water per liter

  // Rename these just so we don't have to use cumbersome pointers
  nSS = *numSS;
  cTotal = *numTotal;

  // Record the number of monomer types including those with zero conc.
  *numSS0 = nSS;

  // Find out if we need to explicitly consider permutations
  if (sumint(numPermsArray,cTotal) > cTotal) {
    noPerms = 0;
  }
  else {
    noPerms = 1;
  }

  // Allocate memory for A, G, and x0
  // THESE ARE NOT FREED UNTIL THE END OF MAIN
  *A = (int **) malloc(nSS * sizeof(int *));
  for (i = 0; i < nSS; i++) {
    (*A)[i] = (int *) malloc(cTotal * sizeof(int));
  }
  *G = (double *) malloc(cTotal * sizeof(double));
  *CompIDArray = (int *) malloc(cTotal * sizeof(int));
  *PermIDArray = (int *) malloc(cTotal * sizeof(int));
  // For this PermIDArray is all zeros
  for (j = 0; j < cTotal; j++) {
    (*PermIDArray)[j] = 0;
  }

  *x0 = (double *) malloc(nSS * sizeof(double));
  *concentrations = (double *) malloc(nSS * sizeof(double));

  // Allocate memory for the struct
  InputStruct = (struct CompStruct *) malloc(cTotal * sizeof(struct CompStruct));
  for (j = 0; j < cTotal; j++) {
    InputStruct[j].Aj = (int *) malloc (nSS * sizeof(int));
  }

  // Allocate memory for the partition functions and initialize
  // We do this even if noPerms == 1 so the compiler doesn't give a warning when
  // optimization if turned on.
  Q = (long double *) malloc (cTotal * sizeof(long double));
  for (j = 0; j < cTotal; j++) {
    Q[j] = 0.0;
  }


  /* read concentrations
   */
  char nupack_concentration[MAXLINE];
  for(int x=0 ; x<nSS ; ++x){
    printf("Enter concentration %d: ", x+1);
    scanf("%s", nupack_concentration);
    tok = strtok(nupack_concentration, tokseps);
    (*x0)[x] = str2double(tok);
    (*concentrations)[x] = (*x0)[x];
  }


  /* read temperature
   */
  char nupack_temperature[MAXLINE];
  printf("Enter temperature: ");
  scanf("%s", nupack_temperature);
  tok = strtok(nupack_temperature, tokseps);
  *temperature = str2double(tok);
  *kT = kB*((*temperature) + ZERO_C_IN_KELVIN);


  /* HARDCODE OCX STARTS */
  char ocx_sep[] = ",";
  char* ocx = malloc(sizeof(char) * MAXLINE);
  for(int x=0 ; x<cTotal ; ++x){
    InputStruct[x].CompID = (x + 1);
    switch(x){
      case 0:
        strcpy(ocx, "1,0,0,-7.92078773e+00");
        break;
      case 1:
        strcpy(ocx, "0,1,0,-9.79502400e+00");
        break;
      case 2:
        strcpy(ocx, "0,0,1,-9.79502400e+00");
        break;
      case 3:
        strcpy(ocx, "1,1,0,-4.84277745e+01");
        break;
      case 4:
        strcpy(ocx, "1,0,1,-4.84277745e+01");
        break;
      case 5:
        strcpy(ocx, "1,1,1,-6.36285141e+01");
        break;
    }
    for(int y=0 ; y<nSS ; ++y){
      if(y < 1){
        InputStruct[x].Aj[y] = atoi(strtok(ocx, ocx_sep));
      }
      else{
        InputStruct[x].Aj[y] = atoi(strtok(NULL, ocx_sep));
      }
    }
    Gperm = str2double(strtok(NULL, ocx_sep))/(*kT);
    InputStruct[x].FreeEnergy = Gperm;
    InputStruct[x].numSS = nSS;
  }

  free(ocx);
  /* HARDCODE OCX ENDS */


  // Compute and enter free energies
  if (noPerms == 0) {
    for (j = 0; j < cTotal; j++) {
      InputStruct[j].FreeEnergy = -(double) logl(Q[j]);
    }
  }
  /* *************************************************************** */


  // Make the matrix A and free energy G and the complex ID list
  for (j = 0; j < cTotal; j++) {
    for (i = 0; i < nSS; i++) {
      (*A)[i][j] = InputStruct[j].Aj[i];
    }
    (*G)[j] = InputStruct[j].FreeEnergy;
    (*CompIDArray)[j] = InputStruct[j].CompID;
  }

  // Free the struct
  for (j = 0; j < cTotal; j++) {
    free(InputStruct[j].Aj);
  }
  free(InputStruct);

  // Free the partition functions
  free(Q);
  
  // Do a quick check of the free energies.  If any are > 0, it's likely there's
  // an input error.  Let the user know if this is the case.
  if (quiet == 0) {
    j = 0;
    while (j < cTotal && (*G)[j] <= 0.0001) {
      j++;
    }
    if (j < cTotal) {
      printf("\n\nWarning: At least one free energy is > 0. %lf\n", (*G)[j]);
      printf("It is likely there is an input error.\n");
      printf("If there is such an error, the the program will still run\n");
      printf("normally and give results, which may be nonsensical.\n");
      printf("If your input file suffix is .cx, .cx-epairs, or .cx-mfe,\n");
      printf("there should be no ordered complex identifier in the input file.\n");
      printf("\n\n");
    }
  }


  /* ************** BEGIN    REFORMATTING PROBLEM *********************** */
  /* 
     This section of the code reformats the problem if there are zero
     entries in the con file.  I.e., if there is a variable whose
     initial concentration is entered as zero, the problem is
     reformulated as if that strand doesn't exist and if a complex
     cannot be formed, the problem is reformulated as if it doesn't
     exist.  Arrays are reallocated to the adjusted size of the
     problem
  */
  
  zerox0 = (int *) malloc(nSS * sizeof(int));
  nonzerox0 = (int *) malloc(nSS * sizeof(int));

  // Check to see if any of the concentrations are zero.  If so, 
  // reformulate the problem accordingly.
  j = 0;  // How many entries are zero
  k = 0;  // How many entries are nonzero
  for (i = 0; i < nSS; i++) {
    if ((*x0)[i] <= DBL_MIN) {
      zerox0[j] = i;
      j++;
    }
    else {
      nonzerox0[k] = i;
      k++;
    }
  }

  if (j > 0) { // Have to reformulate
  
    newnumSS = nSS - j;

    // First count how many entries we have.  notOK = 1 if complex contains something
    // with zero initial concentration
    newnumTotal = 0;
    
    for (j = 0; j < cTotal; j++) {
      notOK = 0;
      for (i = 0; i < nSS-newnumSS; i++) {
        if ((*A)[zerox0[i]][j] > 0) {
          notOK = 1;
        }
      }
      if (notOK == 0) {
        newnumTotal++;
      }
    }

    // Allocate memory for new arrays
    newA = (int **) malloc(newnumSS * sizeof(int *));
    for (i = 0; i < newnumSS; i++) {
      newA[i] = (int *) malloc(newnumTotal * sizeof(int));
    }
    newG =  (double *) malloc(newnumTotal * sizeof(double));
    newx0 = (double *) malloc(newnumSS * sizeof(double));
    newCompIDArray = (int *) malloc(newnumTotal * sizeof(int));
    newPermIDArray = (int *) malloc(newnumTotal * sizeof(int));
    
    // Put in the new x0
    for (i = 0; i < newnumSS; i++) {
      newx0[i] = (*x0)[nonzerox0[i]];
    }
    
    // Go through and pick out the entries to keep
    k = 0;
    for (j = 0; j < cTotal; j++) {
      notOK = 0;
      for (i = 0; i < nSS-newnumSS; i++) {
        if ((*A)[zerox0[i]][j] > 0) {
          notOK = 1;
        }
      }
      if (notOK == 0) {
        for (i = 0; i < newnumSS; i++) {
          newA[i][k] = (*A)[nonzerox0[i]][j];
        }
        newG[k] = (*G)[j];
        newCompIDArray[k] = (*CompIDArray)[j];
        newPermIDArray[k] = (*PermIDArray)[j];
        k++;
      }
    }
    
    // Change names of "new" variables
    // Rename newA
    for (i = 0; i < nSS; i++) {
      free((*A)[i]);
    }
    free(*A);    
    *A = (int **) malloc(newnumSS * sizeof(int *));
    for (i = 0; i < newnumSS; i++) {
      (*A)[i] = (int *) malloc(newnumTotal * sizeof(int *));
      for (j = 0; j < newnumTotal; j++) {
        (*A)[i][j] = newA[i][j];
      }
      free(newA[i]);
    }
    free(newA);
    
    // Rename newG
    free(*G);
    *G = (double *) malloc(newnumTotal * sizeof(double));
    for (j = 0; j < newnumTotal; j++) {
      (*G)[j] = newG[j];
    }
    free(newG);
    
    // Rename CompIDArray
    free(*CompIDArray);
    free(*PermIDArray);
    *CompIDArray = (int *) malloc(newnumTotal * sizeof(int));
    *PermIDArray = (int *) malloc(newnumTotal * sizeof(int));
    for (j = 0; j < newnumTotal; j++) {
      (*CompIDArray)[j] = newCompIDArray[j];
      (*PermIDArray)[j] = newPermIDArray[j];
    }
    free(newCompIDArray);
    free(newPermIDArray);
    
    // Rename newx0
    free(*x0);
    (*x0) = (double *) malloc(newnumSS * sizeof(double));
    for (i = 0; i < newnumSS; i++) {
      (*x0)[i] = newx0[i];
    }
    free(newx0);
    
    // Rename numTotal and numSS
    *numTotal = newnumTotal;
    *numSS = newnumSS;
  
  }
  free(zerox0);
  free(nonzerox0);
  /* ************** FINISHED REFORMATTING PROBLEM *********************** */

  // Calculate molarity of water and convert appropriate quantities to the right units
  MolesWaterPerLiter = WaterDensity((*kT)/kB - ZERO_C_IN_KELVIN);
  for (i = 0; i < (*numSS); i++) {
    (*x0)[i] /= MolesWaterPerLiter;
  }

  return MolesWaterPerLiter;

}
/* ******************************************************************************** */


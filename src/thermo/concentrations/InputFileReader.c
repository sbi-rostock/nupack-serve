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


#include "constants.h"
#include "InputFileReader.h" // File with important definitions
#include <float.h>


/* Get the system's size
 */
void getSize(int *numSS, int *numTotal, int *nTotal, int *LargestCompID,
             int **numPermsArray) {

  char newline;


  // read number of complex IDs
  printf("Enter number of complex IDs: ");
  scanf("%d%c", numTotal, &newline);
  *nTotal = *numTotal;
  *LargestCompID = *numTotal;


  // read number of concentrations
  printf("Enter number of different concentrations: ");
  scanf("%d%c", numSS, &newline);


  // allocate memory and set array of number of permutations
  *numPermsArray = malloc(sizeof(int) * (*numTotal));
  for(int x=0 ; x<(*numTotal) ; ++x){
    (*numPermsArray)[x] = 1;
  }
}



/* If one of concentrations is zero, the problem is reformulated as if the
 * corresponding strand does not exist.
 * The input is stored as follows:
 * - A[i][j] is the number of monomers of type i in complex j
 * - G[j] is the free energy of complex j IN UNITS OF kT
 * - x0[i] is the initial mole fraction of monomer species i
 * CompIDArray and PermIDArray store the corresponding complex IDs and
 * Permutation IDs for the entries loaded in the A and G
 */
double ReadInputFiles(int ***A, double **G, int **CompIDArray,
        int **PermIDArray, double **x0, double** concentrations, int *numSS,
        int *numSS0, int *numTotal, int *numPermsArray, double *kT,
        double* temperature, int Toverride, struct InStruct* InputStruct){

  double MolesWaterPerLiter; // moles of water per liter
  char *tok;
  char separators[] = " \t\n";
  int nSS;     // number of single species
  int cTotal;  // number of complexes
  int noPerms; // 1 if permutations are not explicitly considered

  nSS = *numSS;
  cTotal = *numTotal;

  // record the number of monomer types including those with zero conc
  *numSS0 = nSS;

  // find out if we need to explicitly consider permutations
  if(sumint(numPermsArray,cTotal) > cTotal){
    noPerms = 0;
  }
  else{
    noPerms = 1;
  }


  /* memory allocations
   */

  // A
  *A = malloc(sizeof(int*) * nSS);
  for(int i=0 ; i<nSS ; ++i){
    (*A)[i] = malloc(sizeof(int) * cTotal);
  }

  // G
  *G = malloc(sizeof(double) * cTotal);

  // CompIDArray
  *CompIDArray = malloc(sizeof(int) * cTotal);

  // PermIDArray
  *PermIDArray = malloc(sizeof(int) * cTotal);
  for (int i=0 ; i<cTotal ; ++i){
    (*PermIDArray)[i] = 0;
  }

  // x0
  *x0 = malloc(sizeof(double) * nSS);

  // concentrations
  *concentrations = malloc(sizeof(double) * nSS);

  // partition function for complexes Q
  long double* Q = malloc (sizeof(long double) * cTotal);
  for(int i=0 ; i<cTotal ; ++i){
    Q[i] = 0.0;
  }


  /* read concentrations
   */
  char nupack_concentration[MAXLINE];
  for(int x=0 ; x<nSS ; ++x){
    printf("Enter concentration %d: ", x+1);
    scanf("%s", nupack_concentration);
    tok = strtok(nupack_concentration, separators);
    (*x0)[x] = str2double(tok);
    (*concentrations)[x] = (*x0)[x];
  }


  /* read temperature
   */
  char nupack_temperature[MAXLINE];
  printf("Enter temperature: ");
  scanf("%s", nupack_temperature);
  tok = strtok(nupack_temperature, separators);
  *temperature = str2double(tok);
  *kT = kB*((*temperature) + ZERO_C_IN_KELVIN);


  /* HARDCODE OCX STARTS
   */
  char ocx_sep[] = ",";
  char* ocx = malloc(sizeof(char) * MAXLINE);
  for(int x=0 ; x<cTotal ; ++x){
    InputStruct[x].CompID = (x + 1);
    switch(x){
      case 0:
        strcpy(ocx, "1,1,0,0,-7.92078773e+00");
        break;
      case 1:
        strcpy(ocx, "1,0,1,0,-9.79502400e+00");
        break;
      case 2:
        strcpy(ocx, "1,0,0,1,-9.79502400e+00");
        break;
      case 3:
        strcpy(ocx, "1,1,1,0,-4.84277745e+01");
        break;
      case 4:
        strcpy(ocx, "1,1,0,1,-4.84277745e+01");
        break;
      case 5:
        strcpy(ocx, "1,1,1,1,-6.36285141e+01");
        break;
    }
    InputStruct[x].PermID = atoi(strtok(ocx, ocx_sep));
    for(int y=0 ; y<nSS ; ++y){
      InputStruct[x].Aj[y] = atoi(strtok(NULL, ocx_sep));
    }
    double energy_of_permutation = str2double(strtok(NULL, ocx_sep))/(*kT);
    InputStruct[x].FreeEnergy = energy_of_permutation;
    InputStruct[x].numSS = nSS;
  }

  free(ocx);
  /*
   * HARDCODE OCX ENDS */


  // compute and enter free energies
  if(noPerms == 0){
    for(int i=0 ; i<cTotal ; ++i){
      InputStruct[i].FreeEnergy = -(double) logl(Q[i]);
    }
  }


  // make the matrix A and free energy G and the complex ID list
  for(int j=0 ; j<cTotal ; ++j){
    for (int i=0 ; i<nSS ; ++i){
      (*A)[i][j] = InputStruct[j].Aj[i];
    }
    (*G)[j] = InputStruct[j].FreeEnergy;
    (*CompIDArray)[j] = InputStruct[j].CompID;
  }


  // free allocated memory
  free(Q);


  // calculate molarity of water and convert appropriate quantities to the
  // right units
  MolesWaterPerLiter = WaterDensity((*kT)/kB - ZERO_C_IN_KELVIN);
  for(int i=0 ; i<(*numSS); ++i){
    (*x0)[i] /= MolesWaterPerLiter;
  }

  return MolesWaterPerLiter;
}


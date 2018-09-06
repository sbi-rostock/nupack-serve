/*
 * OutputWriter.c is part of the NUPACK software suite
 * Copyright (c) 2007 Caltech. All rights reserved.
 * Coded by: Justin Bois 9/2006
 *
 * For input and output formatting, see associated manual
 */


#include "OutputWriter.h" // Concentrations header file
#include <shared.h> // Concentrations header file
#include "constants.h" // Concentrations header file


// Structure for sorting output that includes permutations
struct PermSortStruct {
  int CompID; // The ID tag associated with a complex
  int PermID; // ID number for permutation
  int *Aj; // Array representing column j of A
  double FreeEnergy; // Free energy for permutation
  double xj; // The concentration of the permutation
  double xjc; // The complex concentration
  int numSS; // Number of entries in Aj (needed for sorting routine)
  char *AuxStr; // String containing auxillary information from input file
};


/* create a string containing the provenance's header informations:
 * - version
 * - command invocation
 * return the length of the generated string
 */
int concentrations_header(char *provenance, int argc, char **argv) {

  char PROVENANCE_STARTS[] = "{ ";
  char FIELD_VERSION[] = "\"version\": \"";
  char FIELD_COMMAND[] = "\"command\": \"";
  char FIELD_ENDS[] = "\"";
  char FIELD_NEXT[] = ", ";

  int len_entry_version;
  int len_entry_command;
  int len_provenance = 0;


  /* retrieve each entry's value, and store it as provenance
   */

  len_provenance += strlen(PROVENANCE_STARTS);
  provenance = strcpy(provenance, PROVENANCE_STARTS);


  /* version
   */

  // retrieve the version's length, and add it to the provenance
  len_entry_version = strlen(FIELD_VERSION) + strlen(NUPACK_VERSION)
        + strlen(FIELD_ENDS) + strlen(FIELD_NEXT);
  len_provenance += len_entry_version;

  // store the version
  provenance = strcat(
                strcat(
                  strcat(
                    strcat(provenance, FIELD_VERSION),
                    NUPACK_VERSION),
                  FIELD_ENDS),
                FIELD_NEXT);


  /* command
   */

  // retrieve the command's length, and add it to the provenance
  int len_nupack_command = 0;
  for(int x=0 ; x<argc ; ++x) {
    if(x<(argc-1)){
      len_nupack_command += strlen(argv[x]) + 1;
    }
    else{
      len_nupack_command += strlen(argv[x]);
    }
  }
  len_entry_command = strlen(FIELD_COMMAND) + len_nupack_command
        + strlen(FIELD_ENDS) + strlen(FIELD_NEXT);
  len_provenance += len_entry_command;

  // retrieve the command's value
  char nupack_command[len_nupack_command];
  char *n = strcpy(nupack_command, "");
  for(int x=0 ; x<argc ; ++x) {
    n = strcat(n, argv[x]);
    if(x < (argc-1)){
      n = strcat(n, " ");
    }
  }
  // store the command
  provenance = strcat(
                strcat(
                  strcat(
                    strcat(provenance, FIELD_COMMAND),
                    nupack_command),
                  FIELD_ENDS),
                FIELD_NEXT);

  return len_provenance;
}



/* create a string containing all provenance's parameter informations:
 * - temperature
 * - concentrations
 * return the length of the generated string
 */
int concentrations_parameters(char* provenance, int numSS,
        double* concentrations, double temperature){

  char FIELD_CONCENTRATIONS[] = "\"concentrations (M)\": ";
  char FIELD_TEMPERATURE[] = "\"temperature (C)\": ";
  char FIELD_NEXT[] = ", ";
  char LIST_STARTS[] = "[";
  char LIST_ENDS[]   = "]";
  char COMMA[] = ",";


  int len_entry_concentrations;
  int len_entry_temperature;
  int len_provenance = 0;


  /* retrieve each entry's value, and store it as provenance
   */

  provenance = strcpy(provenance, "");


  /* temperature
   */

  // retrieve the temperature's length, and add it to the provenance
  int len_nupack_temperature = (snprintf(NULL, 0, "%.1f", temperature) + 1);
  len_entry_temperature = strlen(FIELD_TEMPERATURE) + len_nupack_temperature
        + strlen(FIELD_NEXT);
  len_provenance += len_entry_temperature;

  // retrieve the temperature's value
  char nupack_temperature[len_nupack_temperature];
  snprintf(nupack_temperature, len_nupack_temperature, "%.1f", temperature);

  // store the temperature
  provenance = strcat(
                  strcat(
                    strcat(provenance, FIELD_TEMPERATURE),
                    nupack_temperature),
                  FIELD_NEXT);


  /* concentrations
   */

  // retrieve the concentrations' length, and add it to the provenance
  int len_nupack_concentrations = strlen(LIST_STARTS);
  for(int j=0 ; j<numSS ; ++j) {
    int len_nupack_concentration = (
        snprintf(NULL, 0, "%e", concentrations[j]) + 1);
    len_nupack_concentrations += len_nupack_concentration;
    if(j < (numSS - 1)){
      len_nupack_concentrations += strlen(COMMA);
    }
  }
  len_nupack_concentrations += strlen(LIST_ENDS);

  len_entry_concentrations = strlen(FIELD_CONCENTRATIONS)
        + len_nupack_concentrations + strlen(FIELD_NEXT);

  len_provenance += len_entry_concentrations;

  // retrieve the concentrations' value
  char *nupack_concentrations = malloc(sizeof(char) * len_nupack_concentrations);
  if(!nupack_concentrations){
    exit(1);
  }
  for(size_t x=0 ; x < len_nupack_concentrations ; ++x){
    nupack_concentrations[x] = 0;
  }
  nupack_concentrations = strcpy(nupack_concentrations, LIST_STARTS);

  for(int j=0 ; j<numSS ; ++j) {
    int len_nupack_concentration = (
        snprintf(NULL, 0, "%e", concentrations[j]) + 1);
    char nupack_concentration[len_nupack_concentration];
    snprintf(nupack_concentration, len_nupack_concentration, "%e",
        concentrations[j]);

    nupack_concentrations = strcat(nupack_concentrations, nupack_concentration);

    if(j < (numSS - 1)){
      nupack_concentrations = strcat(nupack_concentrations, COMMA);
    }
  }
  nupack_concentrations = strcat(nupack_concentrations, LIST_ENDS);

  // store the concentrations
  provenance = strcat(
                strcat(
                  strcat(provenance, FIELD_CONCENTRATIONS),
                  nupack_concentrations),
                FIELD_NEXT);


  // free all memory objects
  free(nupack_concentrations);
  nupack_concentrations = NULL;

  return len_provenance;
}



/* Writes the output of a concentrations calculation.
 * We assume no strand has been associated to a zero concentration, therefore
 * we avoid re-reading the input to check whether the problem has changed
 */
void WriteOutput(double *X, double *G, int *CompIDArray, int LargestCompID,
        int numSS, int numTotal, int nTotal, double kT,
        int SortOutput, double MolesWaterPerLiter, int quiet,
        int NoPermID, int NUPACK_VALIDATE){

  int *CompLookup;
  struct PermSortStruct *PermOutStruct; // output structure for sorting


  // Allocate memory for lookup table for complexes (to check if conc. is nonzero)
  CompLookup = (int *) malloc (LargestCompID * sizeof(int));
  // Initialize it to -1
  for (int j=0 ; j<LargestCompID ; ++j) {
    CompLookup[j] = -1;
  }
  // If included in problem change entry
  for (int j=0; j<numTotal ; ++j) {
    CompLookup[CompIDArray[j]-1] = j;
  }

  // Allocate memory for the output structure
  PermOutStruct = malloc(sizeof(struct PermSortStruct) * nTotal);
  for(int x=0 ; x<nTotal ; ++x){
    PermOutStruct[x].Aj = malloc(sizeof(struct PermSortStruct) * numSS);
  }


  /* HARDCODE OCX STARTS */
  char ocx_sep[] = ",";
  char* ocx = malloc(sizeof(char) * MAXLINE);
  for(int x=0 ; x<nTotal ; ++x){
    PermOutStruct[x].CompID = (x + 1);
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
    for(int y=0 ; y<numSS ; ++y){
      if(y < 1){
        PermOutStruct[x].Aj[y] = atoi(strtok(ocx, ocx_sep));
      }
      else{
        PermOutStruct[x].Aj[y] = atoi(strtok(NULL, ocx_sep));
      }
    }
    double Gperm = str2double(strtok(NULL, ocx_sep))/kT;
    PermOutStruct[x].FreeEnergy = Gperm;
    PermOutStruct[x].AuxStr = malloc(sizeof(char) * 1);
    PermOutStruct[x].AuxStr[0] = '\0';
    PermOutStruct[x].numSS = numSS;

    if (CompLookup[x] == -1) {
      PermOutStruct[x].xj = 0;
      PermOutStruct[x].xjc = 0;
    }
    else {
      PermOutStruct[x].xjc = X[CompLookup[x]];
      PermOutStruct[x].xj = (
          X[CompLookup[x]]
          * exp(
              G[CompLookup[x]] - PermOutStruct[x].FreeEnergy));
    }
  }

  free(ocx);
  /* HARDCODE OCX ENDS */


  // Sort the results (no sorting necessary for SortOutput == 0)
  if (SortOutput == 1) {
    qsort(PermOutStruct,nTotal,sizeof(struct PermSortStruct),Compare11);
  }
  else if (SortOutput == 2) {
    qsort(PermOutStruct,nTotal,sizeof(struct PermSortStruct),Compare12);
  }
  else if (SortOutput == 3) {
    qsort(PermOutStruct,nTotal,sizeof(struct PermSortStruct),Compare13);
  }
  else if (SortOutput == 4) {
    qsort(PermOutStruct,nTotal,sizeof(struct PermSortStruct),Compare14);
  }

  for (int k=0 ; k<nTotal ; ++k){
    // ComplexID
    printf("%d\t", PermOutStruct[k].CompID);
    if (NoPermID == 0) {
      printf("%d\t", PermOutStruct[k].PermID);
    }
    // The corresponding column in A
    for (int i=0 ; i<numSS ; ++i){
      printf("%d\t", PermOutStruct[k].Aj[i]);
    }
    // Free energy in kcal/mol
    if(!NUPACK_VALIDATE) {
      printf("%8.6e\t", PermOutStruct[k].FreeEnergy * kT);
      // Concentration in MOLAR
      printf("%8.6e\t", PermOutStruct[k].xj * MolesWaterPerLiter);
    } else {
      printf("%.14e\t", PermOutStruct[k].FreeEnergy * kT);
      // Concentration in MOLAR
      printf("%.14e\t", PermOutStruct[k].xj * MolesWaterPerLiter);
    }
    printf("%s\n", PermOutStruct[k].AuxStr);
  }

  // Free the allocated memory
  for (int k=0 ; k<nTotal ; ++k){
    free(PermOutStruct[k].Aj);
    free(PermOutStruct[k].AuxStr);
  }
  free(PermOutStruct);

  // Free CompLookup
  free(CompLookup);
}



/* ******************************************************************************** */
int Compare11(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.

     Whichever has highest concentration (xj) returns -1.

     If the concentrations are the same (usually only the case if they're both zero),
     they are sorted by complex ID and then permutation ID.
  */

  const struct PermSortStruct *ps1 = p1;  // Get the right type of pointer
  const struct PermSortStruct *ps2 = p2;

  if (ps1->xj < ps2->xj) {
    return 1;
  }
  else if (ps1->xj > ps2->xj) {
    return -1;
  }
  else { // Equal permutation concentration
    if (ps1->CompID < ps2->CompID) {
      return -1;
    }
    else if (ps1->CompID > ps2->CompID) {
      return 1;
    }
    else {  // same complex ID
      if (ps1->PermID < ps2->PermID) {
	return -1;
      }
      else if (ps1->PermID > ps2->PermID) {
	return 1;
      }
      else { // Shouldn't ever get here
	return 0;
      }
    }
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int Compare12(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.

     Used for sorting first by complex concentration then permutation concentration
  */

  const struct PermSortStruct *ps1 = p1;  // Get the right type of pointer
  const struct PermSortStruct *ps2 = p2;


  if (ps1->CompID == ps2->CompID) { // Same complex, sort by perm conc.
    if (ps1->xj < ps2->xj) {
      return 1;
    }
    else if (ps1->xj > ps2->xj) {
      return -1;
    }
    else { // Same permutation concentration (sort by Perm ID)
      if (ps1->PermID < ps2->PermID) {
	return -1;
      }
      else if (ps1-> PermID > ps2->PermID) {
	return 1;
      }
      else { // Same PermID
	return 0;
      }
    }
  }
  else { // Different complexes, sort by complex conc.
    if (ps1->xjc < ps2->xjc) {
      return 1;
    }
    else if (ps1->xjc > ps2->xjc) {
      return -1;
    }
    else { // Same complex concentration (sort by complex ID)
      if (ps1->CompID < ps2->CompID) {
	return -1;
      }
      else if (ps1-> CompID > ps2->CompID) {
	return 1;
      }
      else { // Shouldn't ever get here
	return 0;
      }
    }
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int Compare13(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.

     They are sorted by complex ID and then permutation ID.
  */

  const struct PermSortStruct *ps1 = p1;  // Get the right type of pointer
  const struct PermSortStruct *ps2 = p2;

  if (ps1->CompID < ps2->CompID) {
    return -1;
  }
  else if (ps1->CompID > ps2->CompID) {
    return 1;
  }
  else {  // same complex ID
    if (ps1->PermID < ps2->PermID) {
      return -1;
    }
    else if (ps1->PermID > ps2->PermID) {
      return 1;
    }
    else { // Shouldn't ever get here
      return 0;
    }
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int Compare14(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.

     Used for sorting first by number of strands in complex, and then
     alphabetically.  I.e., an order of complexes might be:
     A, B, AA, AB, BB

     Then it's sorted by permutation ID number.
  */

  int i; // Counter
  const struct PermSortStruct *ps1 = p1;  // Get the right type of pointer
  const struct PermSortStruct *ps2 = p2;
  int s1;
  int s2;

  if (ps1->CompID == ps2->CompID) {  // Same complex
    if (ps1->PermID < ps2->PermID) {
      return -1;
    }
    if (ps1->PermID > ps2->PermID) {
      return 1;
    }
    else { // Same Perm ID number
      return 0;
    }
  }
  else { // Different complexes
    s1 = sumint(ps1->Aj,ps1->numSS);
    s2 = sumint(ps2->Aj,ps2->numSS);
    if (s1 > s2) {
      return 1;
    }
    else if (s1 < s2) {
      return -1;
    }
    else { // Have same number of strands
      for (i = 0; i < ps1->numSS; i++) {
	if (ps1->Aj[i] > ps2->Aj[i]) {
	  return -1;
	}
	else if (ps1->Aj[i] < ps2->Aj[i]) {
	  return 1;
	}
      }
      return 0; // Shouldn't ever get here
    }
  }

}
/* ******************************************************************************** */



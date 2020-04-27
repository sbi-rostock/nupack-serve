/*
 * OutputWriter.c is part of the NUPACK software suite
 * Copyright (c) 2007 Caltech. All rights reserved.
 * Coded by: Justin Bois 9/2006
 */


#include "constants.h"
#include "OutputWriter.h"
#include <shared.h>


// structure for sorting output that includes permutations
struct PermSortStruct {
  int numSS; // number of entries in Aj (needed for sorting routine)
  int CompID; // ID tag associated with a complex
  int PermID; // ID number for permutation
  int *Aj; // array representing column j of A
  double FreeEnergy; // free energy for permutation
  double xj; // concentration of the permutation
  double xjc; // complex concentration
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

  char FIELD_CONCENTRATIONS[] = "\"monomer concentrations (M)\": ";
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


  // free allocated memory
  free(nupack_concentrations);
  nupack_concentrations = NULL;

  return len_provenance;
}



/* create a string containing all provenance's concentrations informations:
 * - complex ID
 * - permutation ID
 * - strand code (one per provided strand)
 * - concentrations
 * - free energy
 * return the length of the generated string
 */
int concentrations_results(char* provenance, int numSS, int nTotal, double kT,
        double MolesWaterPerLiter, int NUPACK_VALIDATE,
        struct InStruct* inStruct){

  char FIELD_CONCENTRATIONS[] = "\"complex concentrations\": ";
  char PROVENANCE_ENDS[] = " }\n";
  char LIST_STARTS[] = "[";
  char LIST_ENDS[]   = "]";
  char PAIR_STARTS[] = "[";
  char PAIR_ENDS[]   = "]";
  char COMMA[] = ",";


  int len_entry_concentrations = 0;
  int len_provenance = 0;


  /* retrieve each concentrations' entry values, and store it as provenance
   */


  // start provenance
  len_provenance += (strlen(FIELD_CONCENTRATIONS) + strlen(LIST_STARTS));
  provenance = strcat(
                strcpy(provenance, FIELD_CONCENTRATIONS),
                LIST_STARTS);


  for(int i=0 ; i<nTotal ; ++i){

    // start pair
    len_entry_concentrations += strlen(PAIR_STARTS);
    provenance = strcat(provenance, PAIR_STARTS);


    /* complex ID
     */

    // length
    int len_nupack_complex_id = (
        snprintf(NULL, 0, "%d", inStruct[i].CompID) + 1);
    len_entry_concentrations += (len_nupack_complex_id + strlen(COMMA));
    // value
    char nupack_complex_id[len_nupack_complex_id];
    snprintf(nupack_complex_id, len_nupack_complex_id, "%d",
        inStruct[i].CompID);
    // store
    provenance = strcat(
                    strcat(provenance, nupack_complex_id),
                    COMMA);


    /* permutation
     */

    // length
    int len_nupack_permutation_id = (
        snprintf(NULL, 0, "%d", inStruct[i].PermID) + 1);
    len_entry_concentrations += (len_nupack_permutation_id + strlen(COMMA));
    // value
    char nupack_permutation_id[len_nupack_permutation_id];
    snprintf(nupack_permutation_id, len_nupack_permutation_id, "%d",
        inStruct[i].PermID);
    // store
    provenance = strcat(
                    strcat(provenance, nupack_permutation_id),
                    COMMA);


    /* A
     */

    int len_nupack_set_number = 0;
    for(int j=0 ; j<numSS ; ++j){
      // length
      len_nupack_set_number = (
        snprintf(NULL, 0, "%d", inStruct[i].Aj[j]) + 1);
      len_entry_concentrations += (len_nupack_set_number + strlen(COMMA));
      //value
      char nupack_set_number[len_nupack_set_number];
      snprintf(nupack_set_number, len_nupack_set_number, "%d",
        inStruct[i].Aj[j]);
      // store
      provenance = strcat(
                    strcat(provenance, nupack_set_number),
                    COMMA);
    }


    /* energy (Kcal/mol)
     */

    int len_nupack_energy = 0;
    if(!NUPACK_VALIDATE){
      // length
      len_nupack_energy = (
        snprintf(NULL, 0, "%8.6e", ((inStruct[i].FreeEnergy) * kT)) + 1);
      len_entry_concentrations += (len_nupack_energy + strlen(COMMA));
      // value
      char nupack_energy[len_nupack_energy];
      snprintf(nupack_energy, len_nupack_energy, "%8.6e",
        ((inStruct[i].FreeEnergy) * kT));
      // store
      provenance = strcat(
                    strcat(provenance, nupack_energy),
                    COMMA);
    }else{
      // length
      len_nupack_energy = (
        snprintf(NULL, 0, "%.14e", ((inStruct[i].FreeEnergy) * kT)) + 1);
      len_entry_concentrations += (len_nupack_energy + strlen(COMMA));
      // value
      char nupack_energy[len_nupack_energy];
      snprintf(nupack_energy, len_nupack_energy, "%.14e",
        ((inStruct[i].FreeEnergy) * kT));
      // store
      provenance = strcat(
                    strcat(provenance, nupack_energy),
                    COMMA);
    }


    /* concentration (M)
     */

    int len_nupack_concentration = 0;
    if(!NUPACK_VALIDATE){
      // length
      len_nupack_concentration = (
        snprintf(NULL, 0, "%8.6e", ((inStruct[i].xj) * MolesWaterPerLiter)) + 1);
      len_entry_concentrations += (len_nupack_concentration + strlen(PAIR_ENDS));
      // value
      char nupack_concentration[len_nupack_concentration];
      snprintf(nupack_concentration, len_nupack_concentration, "%8.6e",
        ((inStruct[i].xj) * MolesWaterPerLiter));
      // store
      provenance = strcat(provenance, nupack_concentration);
    }else{
      // length
      len_nupack_concentration = (
        snprintf(NULL, 0, "%.14e", ((inStruct[i].xj) * MolesWaterPerLiter)) + 1);
      len_entry_concentrations += (len_nupack_concentration + strlen(PAIR_ENDS));
      // value
      char nupack_concentration[len_nupack_concentration];
      snprintf(nupack_concentration, len_nupack_concentration, "%.14e",
        ((inStruct[i].xj) * MolesWaterPerLiter));
      // store
      provenance = strcat(provenance, nupack_concentration);
    }

    // close pair
    if(i < (nTotal - 1)){
      len_entry_concentrations += (strlen(PAIR_ENDS) + strlen(COMMA));
      provenance = strcat(
                    strcat(provenance, PAIR_ENDS),
                    COMMA);
    }else{
      len_entry_concentrations += strlen(PAIR_ENDS);
      provenance = strcat(provenance, PAIR_ENDS);
    }
    len_provenance += len_entry_concentrations;
  }

  // end provenance
  len_provenance += (strlen(LIST_ENDS) + strlen(PROVENANCE_ENDS));
  provenance = strcat(
                strcat(provenance, LIST_ENDS),
                PROVENANCE_ENDS);

  return len_provenance;
}



/* Writes the output of a concentrations calculation.
 * We assume no strand has been associated to a zero concentration, therefore
 * we avoid re-reading the input to check whether the problem has changed
 */
void WriteOutput(double *X, double *G, int *CompIDArray, int LargestCompID,
        int numSS, int numTotal, int nTotal, double kT, int SortOutput,
        double MolesWaterPerLiter, int NUPACK_VALIDATE,
        struct InStruct* InputStruct){

  int *CompLookup;


  // allocate memory for lookup table for complexes (to check if conc. is nonzero)
  CompLookup = malloc (sizeof(int) * LargestCompID);
  for(int j=0 ; j<LargestCompID ; ++j){ // initialize it to -1
    CompLookup[j] = -1;
  }
  // if included in problem change entry
  for(int j=0; j<numTotal ; ++j){
    CompLookup[CompIDArray[j]-1] = j;
  }


  for(int x=0 ; x<nTotal ; ++x){
    if(CompLookup[x] == -1){
      InputStruct[x].xj = 0;
      InputStruct[x].xjc = 0;
    }
    else{
      InputStruct[x].xjc = X[CompLookup[x]];
      InputStruct[x].xj = (X[CompLookup[x]]
        * exp(G[CompLookup[x]] - InputStruct[x].FreeEnergy));
    }
  }


  // sort the results
  if(SortOutput == 1){
    qsort(InputStruct, nTotal, sizeof(struct InStruct), Compare11);
  }
  else if(SortOutput == 2){
    qsort(InputStruct, nTotal, sizeof(struct InStruct), Compare12);
  }
  else if(SortOutput == 3){
    qsort(InputStruct, nTotal, sizeof(struct InStruct), Compare13);
  }
  else if(SortOutput == 4){
    qsort(InputStruct, nTotal, sizeof(struct InStruct), Compare14);
  }


  // free allocated memory
  free(CompLookup);
}



/* Comparison function (in mandatory form) to send to qsort. See Prata, C
 * Primer Plus, 4th Ed. p. 654 for description.
 * Whichever has highest concentration (xj) returns -1.
 * If the concentrations are the same (usually only the case if they're both
 * zero), they are sorted by complex ID and then permutation ID
 */
int Compare11(const void *p1, const void *p2) {

  const struct PermSortStruct *ps1 = p1;
  const struct PermSortStruct *ps2 = p2;

  if (ps1->xj < ps2->xj) {
    return 1;
  }
  else if (ps1->xj > ps2->xj) {
    return -1;
  }
  else { // equal permutation concentration
    if (ps1->CompID < ps2->CompID) {
      return -1;
    }
    else if (ps1->CompID > ps2->CompID) {
      return 1;
    }
    else { // same complex ID
      if (ps1->PermID < ps2->PermID) {
        return -1;
      }
      else if (ps1->PermID > ps2->PermID) {
        return 1;
      }
      else { // shouldn't ever get here
        return 0;
      }
    }
  }
}



/* Comparison function (in mandatory form) to send to qsort. See Prata, C
 * Primer Plus, 4th Ed. p. 654 for description.
 * Used for sorting first by complex concentration, then by permutation
 * concentration
 */
int Compare12(const void *p1, const void *p2) {

  const struct PermSortStruct *ps1 = p1;
  const struct PermSortStruct *ps2 = p2;


  if (ps1->CompID == ps2->CompID) { // same complex, sort by perm conc
    if (ps1->xj < ps2->xj) {
      return 1;
    }
    else if (ps1->xj > ps2->xj) {
      return -1;
    }
    else { // same permutation concentration (sort by Perm ID)
      if (ps1->PermID < ps2->PermID) {
        return -1;
      }
      else if (ps1-> PermID > ps2->PermID) {
        return 1;
      }
      else { // same PermID
        return 0;
      }
    }
  }
  else { // different complexes, sort by complex conc
    if (ps1->xjc < ps2->xjc) {
      return 1;
    }
    else if (ps1->xjc > ps2->xjc) {
      return -1;
    }
    else { // same complex concentration (sort by complex ID)
      if (ps1->CompID < ps2->CompID) {
        return -1;
      }
      else if (ps1-> CompID > ps2->CompID) {
        return 1;
      }
      else { // shouldn't ever get here
        return 0;
      }
    }
  }
}



/* Comparison function (in mandatory form) to send to qsort. See Prata, C
 * Primer Plus, 4th Ed. p. 654 for description.
 * They are sorted by complex ID, and then permutation ID
 */
int Compare13(const void *p1, const void *p2) {

  const struct PermSortStruct *ps1 = p1;
  const struct PermSortStruct *ps2 = p2;


  if (ps1->CompID < ps2->CompID) {
    return -1;
  }
  else if (ps1->CompID > ps2->CompID) {
    return 1;
  }
  else { // same complex ID
    if (ps1->PermID < ps2->PermID) {
      return -1;
    }
    else if (ps1->PermID > ps2->PermID) {
      return 1;
    }
    else { // shouldn't ever get here
      return 0;
    }
  }
}



/* Comparison function (in mandatory form) to send to qsort. See Prata, C
 * Primer Plus, 4th Ed. p. 654 for description.
 * Used for sorting first by number of strands in complex, and then
 * alphabetically, i.e., an order of complexes might be: A, B, AA, AB, BB
 * later on sorted by permutation ID number
 */
int Compare14(const void *p1, const void *p2) {

  const struct PermSortStruct *ps1 = p1;
  const struct PermSortStruct *ps2 = p2;
  int s1;
  int s2;


  if (ps1->CompID == ps2->CompID) { // same complex
    if (ps1->PermID < ps2->PermID) {
      return -1;
    }
    if (ps1->PermID > ps2->PermID) {
      return 1;
    }
    else { // same Perm ID number
      return 0;
    }
  }
  else { // sifferent complexes
    s1 = sumint(ps1->Aj,ps1->numSS);
    s2 = sumint(ps2->Aj,ps2->numSS);
    if (s1 > s2) {
      return 1;
    }
    else if (s1 < s2) {
      return -1;
    }
    else { // have same number of strands
      for(int i=0 ; i<(ps1->numSS) ; ++i){
        if (ps1->Aj[i] > ps2->Aj[i]) {
          return -1;
        }
        else if (ps1->Aj[i] < ps2->Aj[i]) {
          return 1;
        }
      }
      return 0; // shouldn't ever get here
    }
  }
}


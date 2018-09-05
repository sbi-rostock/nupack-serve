/*
   complexes.c is part of the NUPACK software suite
   Copyright (c) 2007 Caltech. All rights reserved.
   Coded by: Robert Dirks, 3/2006 and Justin Bois 1/2007

   Multistranded partition function
   Calculates all distinct complexes and
   circular permutations for N strands.

   Computes and outputs all of the partition functions

   See NUPACK manual for usage notes
*/


#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>

#include <thermo/core.h>

#include "complexesStructs.h"
#include "complexesUtils.h"
#include "permBG.h"
#include "ReadCommandLine.h"

extern int nStrands;
globalArgs_t globalArgs;



int main( int argc, char **argv) {

  char **seqs; // list of all seqs
  int *seqlength; // list of all seqlengths

  multiset *allSets;

  int nSets;
  int totalSets; // including sets in the input list
  int setStart=0; // index of set to start with, for use with -listonly

  int maxLength;

  int maxListComplexSize = 0;
  int maxComplexSize = 0;
  int nNewComplexes = 0;
  int nNewPerms = 0;

  long double pf;

  double totalOrders;
  int nTotalOrders = 0;
  int totalOrders2 = 0;

  long double TEMP_K;

  // permutations
  permutation *currentPerm=NULL;
  int permId;
  char line[MAXLINE];
  char line2[MAXLINE];
  char *token;

  long double **permPr = NULL;

  int lastCxId = 1; // index complex id numbers
                    // (in case some are not used due to no possible secondary
                    // structures)

  char LIST_STARTS[] = "[";
  char LIST_ENDS[] = "]";
  char COMMA[] = ",";


  // provenance blocks
  int len_header = 1000;
  int len_parameters = 1000;
  int len_complexes = 1000;
  int len_provenance;


  // global argument defaults
  globalArgs.T = 37.0;
  globalArgs.dangles = 1;
  globalArgs.dopairs = 0;
  globalArgs.parameters = RNA;
  globalArgs.listonly = 0;
  globalArgs.cutoff = 0.001; // cutoff bp probability to report
  globalArgs.onlyOneMFE = 1;
  globalArgs.sodiumconc = 1.0;
  globalArgs.magnesiumconc = 0.0;
  globalArgs.uselongsalt = 0;
  strcpy(globalArgs.inputFilePrefix, "NoInputFile");


  TEMP_K = globalArgs.T + ZERO_C_IN_KELVIN;


  /* read number of sequences
   */
  char newline;
  printf("Enter number of different sequences: ");
  scanf("%d%c", &nStrands, &newline);

  // allocate function variables
  seqs = (char**) malloc(sizeof(char*) * nStrands);
  seqlength = (int*) malloc(sizeof(int) * nStrands);

  maxLength = 0;

  /* read sequences
   */
  for(int i=0 ; i<=(nStrands-1) ; ++i) {
    printf("Enter sequence %d:\n", i+1);
    scanf("%s", line);
    seqlength[i] = strlen(line);

    if(seqlength[i] > maxLength){
      maxLength = seqlength[i];
    }

    seqs[i] = (char*) malloc(sizeof(char) * (seqlength[i]+1));
    strcpy(seqs[i], line);
  }

  /* read max complex size
   */
  char *q, r[MAXLINE];
  while (fgets(r, MAXLINE, stdin)){
      maxComplexSize = strtol(r, &q, 10);
      if (q == r || *q != '\n') {
        printf("Enter max complex size to completely enumerate: ");
      } else break;
  }


  // read information from .list file
  maxListComplexSize = maxComplexSize;

  nNewPerms = nStrands;

  // determine total # of distinct strand orders (lovasz, 3.23b)
  totalOrders = nNewPerms;

  for(int i=1; i<=maxComplexSize ; ++i){
    for(int j=1 ; j<=i ; ++j){
      totalOrders += pow(nStrands, gcd(j, i))/i;
    }
  }
  nTotalOrders = totalOrders + 0.1;

  // generate all multisets
  nSets = binomial_coefficient(maxComplexSize + nStrands,maxComplexSize) - 1;
  if(nSets < 1) {
    fprintf(stderr,"Integer overflow occurred while counting permutations!\n");
    exit(1);
  }


  totalSets = nSets + nNewComplexes;


  /* generate all necklaces for each length with order nStrands starts
   */
  int cursize = 1;
  permutation * allPermutations = (permutation*)
                    malloc(nTotalOrders * sizeof(permutation));
  int added = 0;
  int offset = 0;
  for(cursize = 1; cursize <= maxComplexSize ; cursize++) {
    added = makePermutations(allPermutations + offset,cursize,nStrands);
    offset += added;
  }
  totalOrders2 = offset;

  /* read permutations
   */
  for(int x=1 ; x<=nNewPerms ; ++x){
    printf("Enter permutation %d: ", x);
    fgets(line, MAXLINE, stdin);

    int curStrand;
    int curStrandIndex;
    int curNumStrands;

    strncpy(line2,line,MAXLINE);
    token = strtok(line," ,\t\n");
    curNumStrands = 0;
    while(NULL != token) {
      curNumStrands++;
      token = strtok(NULL, " ,\t\n");
    }
    allPermutations[offset].nSeqs = curNumStrands;
    allPermutations[offset].code = (int *) malloc(curNumStrands * sizeof(int));
    allPermutations[offset].strand_sums = (int *) malloc(nStrands * sizeof(int));
    allPermutations[offset].symmetryFactor = 1;
    for(curStrand = 0; curStrand < nStrands ; curStrand++) {
      allPermutations[offset].strand_sums[curStrand] = 0;
    }
    token = strtok(line2," ,\t\n");

    curStrandIndex = 0;
    while(NULL != token) {
      sscanf(token, "%d", &curStrand);
      allPermutations[offset].code[curStrandIndex] = curStrand;
      allPermutations[offset].strand_sums[curStrand - 1] ++;
      curStrandIndex ++;
      token = strtok(NULL, " ,\t\n");
    }
    ++offset;
  }
  /*
   * generate all necklaces for each length with order nStrands ends */


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
  len_provenance = complexes_header(header, argc, argv);
  for(int y=0 ; y<len_provenance ; ++y){
    printf("%c", header[y]);
  }

  // free provenance block
  free(header);
  header = NULL;
  /*
   * echo provenance header ends */


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
  len_provenance = complexes_parameters(parameters, nStrands, seqs,
        nTotalOrders);
  for(int y=0 ; y<len_provenance ; ++y){
    printf("%c", parameters[y]);
  }

  // free provenance block
  free(parameters);
  parameters = NULL;
  /*
   * echo provenance parameters ends */


  /* complexes calculation starts
   */
  totalOrders2 = offset;
  qsort(allPermutations, totalOrders2, sizeof(permutation),
        &comparePermutations);

  totalSets = CountSets(allPermutations, totalOrders2, nStrands);

  allSets = (multiset*) malloc(sizeof(multiset) * totalSets);

  int maxSeqLength = FillSets(allSets, allPermutations,
                          totalSets, totalOrders2,
                          nStrands, seqlength) ;

  maxListComplexSize = GetMaxComplexSize(allSets,totalSets);

  int* nicks = (int*) malloc(sizeof(int) * maxListComplexSize);

  // allocate memory for pfSeq;
  char* pfSeq = (char*) malloc(sizeof(char) * (maxSeqLength + 1));

  permPr = (long double**) malloc(sizeof(long double*) * nStrands);

  for(int j=0; j<nStrands ; ++j) { // calloc initialize to zero
    permPr[j] = (long double*) calloc(seqlength[j], sizeof(long double));
  }

  for(int i=setStart ; i<=(totalSets-1) ; ++i){
    allSets[i].nMfePerms = 0;
    allSets[i].mfePerms = (int*) calloc(10, sizeof(int));
  }


  int status = setStart;
  for(int i=setStart ; i<=(totalSets-1) ; ++i){

    status += allSets[i].nPerms;
    allSets[i].pf = 0; // initialize pf


    currentPerm = allSets[i].perms;
    permId = 1;
    while( currentPerm != NULL){
      resetNicks(maxListComplexSize, nicks);

      int seqCode = (currentPerm->code)[0] - 1;
      strcpy(pfSeq, seqs[seqCode]);

      // set sequences and nicks
      if(allSets[i].nSeqs >= 2){
        nicks[0] = seqlength[seqCode] - 1;
      }

      for(int k=0 ; k<=(allSets[i].nSeqs-2) ; ++k){
        seqCode = (currentPerm->code)[k+1] - 1;

        strcat(pfSeq, "+");
        strcat(pfSeq, seqs[seqCode]);

        if(k != (allSets[i].nSeqs-2)){
          nicks[k+1] = nicks[k] + seqlength[seqCode];
        }

      }

      strncpy(currentPerm->seq, pfSeq,
            allSets[i].totalLength + allSets[i].nSeqs);

      // call library function to compute pseudoknot-free partition function
      int tmpLength = strlen(pfSeq); // store current sequence length
      int seqNum[MAXSEQLENGTH + 1];  // store current sequence as ints
      convertSeq(pfSeq, seqNum, tmpLength);
      pf = pfuncFullWithSym(seqNum, 3, globalArgs.parameters,
            globalArgs.dangles, globalArgs.T, globalArgs.dopairs,
            currentPerm->symmetryFactor, globalArgs.sodiumconc,
            globalArgs.magnesiumconc, globalArgs.uselongsalt);

      /* echo provenance complexes starts
       */
      // allocate provenance block
      char *complexes = malloc(sizeof(char) * len_complexes);
      if(!complexes){
        exit(1);
      }
      for(int y=0 ; y<len_complexes ; ++y){
        complexes[y] = 0;
      }

      // fill provenance block
      if(i == setStart){
        len_provenance = complexes_results(complexes, lastCxId, permId, nStrands, allSets, i, pf, TEMP_K, LIST_STARTS);
      } else if ((i > setStart) && (i < (totalSets-1))){
        len_provenance = complexes_results(complexes, lastCxId, permId, nStrands, allSets, i, pf, TEMP_K, COMMA);
      } else{
        len_provenance = complexes_results(complexes, lastCxId, permId, nStrands, allSets, i, pf, TEMP_K, LIST_ENDS);
      }
      for(int y=0 ; y<len_provenance ; ++y){
        printf("%c", complexes[y]);
      }

      // free provenance block
      free(complexes);
      complexes = NULL;
      /*
       * echo provenance complexes ends */

      permId++;

      currentPerm->pf = pf;
      allSets[i].pf += pf;

      currentPerm = currentPerm->next;
    }

    // keep complex Ids consecutive
    if(allSets[i].pf > 0.0){
      lastCxId++;
    }

  }

  for(int j=0 ; j<nStrands ; ++j){ // free
    free(permPr[j]);
    permPr[j] = NULL;
  }
  free(permPr);
  permPr = NULL;

  free( nicks); nicks = NULL;

  for(int i=0 ; i<=(nStrands-1) ; ++i){
    free(seqs[i]);
    seqs[i] = NULL;
  }

  free(seqlength);
  seqlength = NULL;

  free(seqs);
  seqs = NULL;

  for(int i=setStart ; i<=(totalSets-1) ; ++i){
    free(allSets[i].code);
    allSets[i].code = NULL;

    free(allSets[i].mfePerms);
    allSets[i].mfePerms = NULL;

    currentPerm = allSets[i].perms;
    while( currentPerm != NULL){
      free(currentPerm->code); currentPerm->code = NULL;
      free(currentPerm->baseCode); currentPerm->code = NULL;
      free(currentPerm->seq); currentPerm->seq = NULL;
      free(currentPerm->strand_sums); currentPerm->strand_sums = NULL;
      currentPerm = currentPerm->next;
    }
  }
  free(allPermutations); allPermutations = NULL;
  free(allSets); allSets = NULL;

  free(pfSeq);
  pfSeq = NULL;
  /*
   * complexes calculation ends */

  return 0;
}


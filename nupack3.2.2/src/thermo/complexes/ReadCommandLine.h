#include <stdio.h>

#include "complexesStructs.h"

int ReadCommandLine(int, char**);

int ReadInputFileComplexes(char *filePrefix, int *nStrands,
                           char ***seqs, int **seqlength,
                           int *maxLength, int *maxComplexSize);

void printHeader(int nStrands, char **seqs, int maxComplexSize,
                 int totalOrders, int nNewPerms, int nSets, int nNewComplexes,
                 FILE *F_cx, int nargs, char **argv, int isPairs);
int complexes_header(char*, int, char**);
int complexes_parameters(char*, int, char**, int);
int complexes_results(char*, int, int, int, multiset*, int, long double,
        long double, char*);
void print_deprecation_info(FILE *out);

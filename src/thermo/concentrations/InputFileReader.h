#ifndef NUPACK_THERMO_CONCENTRATIONS_INPUTFILEREADER_H__
#define NUPACK_THERMO_CONCENTRATIONS_INPUTFILEREADER_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void getSize(int *numSS, int *numTotal, int *nTotal, int *LargestCompID,
       int **numPermsArray);

double ReadInputFiles(int ***A, double **G, int **CompIDArray,
        int **PermIDArray, double **x0, double** concentrations, int *numSS,
        int *numSS0, int *numTotal, int *numPermsArray, double *kT,
        double* temperature, int Toverride, struct InStruct* InputStruct);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_THERMO_CONCENTRATIONS_INPUTFILEREADER_H__ */

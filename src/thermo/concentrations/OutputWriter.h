#ifndef NUPACK_THERMO_CONCENTRATIONS_OUTPUTWRITER_H__
#define NUPACK_THERMO_CONCENTRATIONS_OUTPUTWRITER_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void WriteOutput(double* X, double* G, int* CompIDArray, int LargestCompID,
        int numSS, int numTotal, int nTotal, double kT, int SortOutput,
        double MolesWaterPerLiter, int NUPACK_VALIDATE,
        struct InStruct* InputStruct);

// provenance functions
int concentrations_header(char*, int, char**);
int concentrations_parameters(char*, int, double*, double);
int concentrations_results(char*, int, int, double kT, double, int,
        struct InStruct*);

// comparison functions for sorting
int Compare11(const void*, const void*);
int Compare12(const void*, const void*);
int Compare13(const void*, const void*);
int Compare14(const void*, const void*);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_THERMO_CONCENTRATIONS_OUTPUTWRITER_H__ */

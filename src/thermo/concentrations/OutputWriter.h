#ifndef NUPACK_THERMO_CONCENTRATIONS_OUTPUTWRITER_H__
#define NUPACK_THERMO_CONCENTRATIONS_OUTPUTWRITER_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void WriteOutput(double*, double*, int*, int, int, int, int, double, int,
        double, int, int, struct InStruct*);

// provenance functions
int concentrations_header(char*, int, char**);
int concentrations_parameters(char*, int, double*, double);
int concentrations_results(char*, int, int, double kT, double, int, int,
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

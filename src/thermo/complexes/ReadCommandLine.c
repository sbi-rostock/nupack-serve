/*
  ReadCommandLine.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks, 3/2006 and Justin Bois 1/2007

  Reads command line input for multiWrapper.c.
  Uses the package getopt.h to retrieve option flags.
*/

#include "ReadCommandLine.h"
#include <string.h>
#include <stdlib.h>
#include <getopt.h> // Takes options from the command line
#include <time.h>
#include <ctype.h>

#include <shared.h>
#include <thermo/core.h>
#include "complexesStructs.h"

extern globalArgs_t globalArgs;
extern char PARAM_FILE[];


/* ********** SCRIPT FOR DISPLAYING HELP ************************** */
void DisplayHelpComplexes(void);

/* ****************************************************************** */
int ReadCommandLine(int nargs, char **args) {

  // Returns 1 if input file is specified and 0 otherwise.

  int options;  // Counters used in getting flags
  int ShowHelp=0; // ShowHelp = 1 if help option flag is selected
  int option_index;
  char line[MAXLINE];
  char param[MAXLINE];
  float temp;

  /* version 3 output */

  static struct option long_options [] =
    {
      {"T",             required_argument,  NULL, 'e'},
      {"dangles",       required_argument,  NULL, 'f'},
      {"material",      required_argument,  NULL, 'g'},
      {"help",          no_argument,        NULL, 'h'},
      {"cutoff",        required_argument,  NULL, 'l'},
      {"degenerate",    no_argument,        NULL, 'n'},
      {"sodium",        required_argument,  NULL, 'o'},
      {"magnesium",     required_argument,  NULL, 'p'},
      {"longhelixsalt", no_argument,        NULL, 'q'},
      {"validate",      no_argument,        NULL, 'r'},
      {0, 0, 0, 0}
    };

  SetExecutionPath(nargs, args);


  NUPACK_VALIDATE=0;
  // Get the option flags
  while (1)
    {

      /* getopt_long stores the option index here. */
      option_index = 0;
      options = getopt_long_only (nargs, args, 
				  "abcde:f:g:h:ijkl:mno:p:qrs", 
                                  long_options, 
				  &option_index);
		
      // Detect the end of the options.
      if (options == -1)
        break;

      switch (options) {

      case 'e':
        strcpy( line, optarg);
        if( sscanf(line, "%f", &(temp)) != 1) {
          printf("Invalid T value\n");
          exit(1);
        }
        globalArgs.T = (long double) temp;
        break;

      case 'f':
        strcpy( line, optarg);
        if( sscanf(line, "%s", param) != 1) {
          printf("Invalid parameters value\n");
          exit(1);
        }
        if (isdigit(param[0])) {
          globalArgs.dangles = atoi(param);
        }
        else {
          if (!strcmp(param,"none")) {
            globalArgs.dangles = 0;
          }
          else if (!strcmp(param,"some")) {
            globalArgs.dangles = 1;
          }
          else if (!strcmp(param,"all")) {
            globalArgs.dangles = 2;
          }
          else {
            printf("Invalid dangles value\n");
            exit(1);
          }
        }
        break;

      case 'g':
        strcpy( line, optarg);
        if( sscanf(line, "%s", param) != 1) {
          printf("Invalid parameters value\n");
          exit(1);
        }
        else if( !strcmp(param, "dna") || !strcmp(param, "dna1998")) {
          globalArgs.parameters = DNA;

        } else if( !strcmp(param, "rna") || !strcmp(param, "rna1995") ) {
          globalArgs.parameters = RNA;
        } else if( !strcmp(param, "rna37") || !strcmp(param, "rna1999") ) {
          if(strcmp(PARAM_FILE,"rna37") == 0) {
            printf("Parameter specification using rna37 has been deprecated. Please use rna1999 instead\n");
          }
          globalArgs.parameters = RNA37;
        } else {
          globalArgs.parameters = USE_SPECIFIED_PARAMETERS_FILE;
          strcpy( PARAM_FILE, param);
        }
        break;

      case 'h':
        ShowHelp = 1;
        break;

      case 'l':
        strcpy( line, optarg);
        if( sscanf(line, "%f", &(temp)) != 1) {
          printf("Invalid cutoff value\n");
          exit(1);
        }
        globalArgs.cutoff = (long double) temp;
        break;

      case 'n':
        globalArgs.onlyOneMFE = 0;
        break;

      case 'o':
	      strcpy( line, optarg);
	      globalArgs.sodiumconc = str2double(line);
	      break;
	
      case 'p':
	      strcpy( line, optarg);
	      globalArgs.magnesiumconc = str2double(line);
	      break;
	
      case 'q':
	      globalArgs.uselongsalt = 1;
	      break;

      case '?':
        // getopt_long already printed an error message.
        break;

      case 'r':
        NUPACK_VALIDATE=1;
        globalArgs.permsOn = 1;
        globalArgs.cutoff = 0.0;

        break;

      default:
        abort ();
      }
    }

  if (ShowHelp) {
    DisplayHelpComplexes();
    exit(1);
  }


  // Check salt inputs to make sure we're ok
  if ((globalArgs.sodiumconc != 1.0 || globalArgs.magnesiumconc != 0.0) && globalArgs.parameters != DNA) {
    printf("WARNING: No salt corrections availabe for RNA.  Using 1 M Na and 0 M Mg.\n");
    globalArgs.sodiumconc = 1.0;
    globalArgs.magnesiumconc = 0.0;
  }

  if (globalArgs.sodiumconc  <= 0.0) {
    printf("ERROR: In valid sodium concentration.  Must have [Na+] > 0.\n");
    exit(1);
  }

  if (globalArgs.magnesiumconc  < 0.0) {
    printf("ERROR: In valid magnesium concentration.  Must have [Mg2+] >= 0.\n");
    exit(1);
  }

  if (globalArgs.sodiumconc < 0.05 || globalArgs.sodiumconc > 1.1) {
    printf("WARNING: Salt correction only verified for 0.05 M < [Na+] < 1.1 M.\n");
    printf("         [Na+] = %g M may give erroneous results.\n",(float) globalArgs.sodiumconc);
  }

  if (globalArgs.magnesiumconc > 0.2) {
    printf("WARNING: Salt correction only verified for [Mg2+] <= 0.2 M.\n");
    printf("         [Mg2+] = %g M may give erroneous results.\n",(float) globalArgs.magnesiumconc);
  }
  // The range of validity of the magnesium correction is unknown

  if (globalArgs.uselongsalt && globalArgs.magnesiumconc > 0.0) {
    printf("WARNING: No magnesium correction parameters are available for the long\n");
    printf("         helix mode of salt correction.  Using [Mg2+] = 0.\n");
    globalArgs.magnesiumconc = 0.0;
  }
  

  // Get the the input file
  if (optind == nargs) { // There's no input from the user
    return 0;
  }
  else {
    strcpy(globalArgs.inputFilePrefix,args[optind]);
  }

  return 1;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void DisplayHelpComplexes() {
  printf("Please read the NUPACK User Guide for detailed instructions.\n");
  printf("Usage: complexes [OPTIONS] PREFIX\n");
  printf("Calculate equilibrium properties of all possible unpseudoknotted complexes\n");
  printf("of the input strands up to user-defined size L_max\n");
  PrintNupackThermoHelp();
  printf("Additional options:\n");
  printf(" -cutoff CUTOFF   set the minimum stored probability/expected value\n");
  printf("\n");
}
/* ******************************************************************************** */

int ReadInputFileComplexes( char *filePrefix, int *nStrands,
                            char ***seqs, int **seqlength,
                            int *maxLength, int *maxComplexSize) {

  // Returns 1 is input file is ok, 0 if it has an error, and 2 if it is not found.

  FILE *F_inp;
  char *token;
  char line[MAXLINE];
  char line2[MAXLINE];
  int i;
  char filename[MAXLINE];

  strcpy( filename, globalArgs.inputFilePrefix);
  strcat( filename, ".in");

  strcpy( filePrefix, globalArgs.inputFilePrefix);

  F_inp = fopen( filename, "r");
  if( !F_inp) {
    return 2;
  }

  fgets( line, MAXLINE, F_inp);
  while( line[0] == '%' || line[0] == '>') {
    fgets( line, MAXLINE, F_inp);
  }

  //Read in # of strands
  token = strtok( line, "%>");
  if( sscanf( token, "%d", nStrands ) != 1) {
    printf("Error in %s: %s\n", filename, token);
    fclose( F_inp);
    return 0;
  }

  //allocate function variables
  (*seqs) = (char **) malloc( *nStrands * sizeof( char*));
  (*seqlength) = (int *) malloc( *nStrands * sizeof( int));

  //read in sequences
  *maxLength = 0;
  for( i = 0; i < *nStrands; i++) {
    fgets( line, MAXLINE, F_inp);
    while( line[0] == '%' || line[0] == '>') {
      fgets( line, MAXLINE, F_inp);
    }
    token = strtok( line, "%>");

    if( sscanf( token, "%s", line2 ) != 1) {
      printf("Error in %s: %s\n", filename, token);
      fclose( F_inp);
      return 0;
    }
    (*seqlength)[i] = strlen( line2);
    if( (*seqlength)[i] > *maxLength) *maxLength = (*seqlength)[i];


    (*seqs)[i] = (char*) malloc( ((*seqlength)[i]+1)*sizeof( char));

    strcpy( (*seqs)[i], line2);
  }

  //read in Lmax
  fgets( line, MAXLINE, F_inp);
  while( line[0] == '%' || line[0] == '>') {
    fgets( line, MAXLINE, F_inp);
  }
  token = strtok( line, "%>");

  if( sscanf( token, "%d", maxComplexSize ) != 1) {
    printf("Error in %s: %s\n", filename, token);
    fclose( F_inp);
    return 0;
  }

  //check if complex list size is included
  char * check = fgets(line, MAXLINE, F_inp);
  if(!check) {
    fclose( F_inp);
    return 1;
  }
  while(check && (line[0] == '%' || line[0] == '>') ) {
    check = fgets( line, MAXLINE, F_inp);
  }
  if (!check) {
    fclose( F_inp);
    return 1;
  }


  fclose( F_inp);
  return 1;
}

/* ******** */
void printHeader( int nStrands, char **seqs, int maxComplexSize,
                  int nTotalOrders, int nNewPerms, int nSets, int nNewComplexes,
                  FILE *F_cx, int nargs, char **argv, int isPairs) {

  char comment[2] = "%";
  char timestr[25];
  time_t curtime;
  struct tm *loctime;
  int i;
  int initial_perms = nTotalOrders - nNewPerms;
  
  curtime = time(NULL); //current time
  loctime = localtime( &curtime);


  fprintf(F_cx,"%s NUPACK %s\n",comment, NUPACK_VERSION);
  fprintf(F_cx,"%s Program: complexes\n",comment);

  strncpy(timestr,asctime(loctime),24);
  timestr[24] = '\0';

  fprintf( F_cx, "%s Start time: %s PST\n%s\n", comment,
           timestr, comment);

  fprintf( F_cx, "%s Command: ", comment);
  for( i = 0; i < nargs; i++) {
    fprintf( F_cx, "%s ", argv[ i]);
  }

  fprintf( F_cx, "\n");

  fprintf( F_cx, "%s Maximum complex size to enumerate: %d\n",
           comment, maxComplexSize);

  // Show that we're using a cutoff if this is a pairs file
  if (isPairs && globalArgs.cutoff > 0.0) {
    fprintf(F_cx, "%s Minimum output pair probability: %Lg\n",
            comment, globalArgs.cutoff);
  }
  if (globalArgs.v3) {
    fprintf( F_cx, "%s Number of complexes from enumeration: %d\n",
             comment, nSets);
    fprintf( F_cx, "%s Additional complexes from .list file: %d\n",
             comment, nNewComplexes);
    fprintf( F_cx, "%s Total number of permutations to calculate: %d\n",
             comment, nTotalOrders );
  }
  else {
    fprintf( F_cx, "%s Number of complexes from enumeration: %d\n",
             comment, initial_perms);
    fprintf( F_cx, "%s Additional complexes from .list file: %d\n",
             comment, nNewPerms);
    fprintf( F_cx, "%s Total number of complexes: %d\n",
             comment, nTotalOrders );
  }

  fprintf( F_cx, "%s Parameters: ", comment);

  if( globalArgs.parameters == DNA) {
    fprintf( F_cx, "DNA, 1998\n");
  }
  else if(  globalArgs.parameters == RNA) {
    fprintf( F_cx, "RNA, 1995\n");
  }
  else if( globalArgs.parameters == RNA37) {
    fprintf( F_cx, "RNA, 1999\n");
  }

  fprintf( F_cx, "%s Dangles setting: %d\n", comment,
           globalArgs.dangles);

  fprintf( F_cx, "%s Temperature (C): %.1f\n", comment,
           (float) globalArgs.T);

  fprintf(F_cx,"%s Sodium concentration: %.4f M\n",comment, (float) globalArgs.sodiumconc);
  fprintf(F_cx,"%s Magnesium concentration: %.4f M\n",comment,(float) globalArgs.magnesiumconc);
  fprintf( F_cx, "%s\n", comment);
  fprintf( F_cx, "%s Do not change the comments below this line, as they may be read by other programs!",
           comment);
  fprintf( F_cx, "\n%s\n%s Number of strands: %d\n", comment, comment,
           nStrands);
  fprintf( F_cx, "%s id sequence\n", comment);
  for( i = 0; i < nStrands; i++) {
    fprintf( F_cx, "%s %2d %s\n", comment, i+1, seqs[i]);
  }

  fprintf( F_cx, "%s T = %.1f\n", comment, (float) globalArgs.T);
}

/* ******** */
void print_deprecation_info(FILE *out) {
  char *dep_mess = 
  "Relative to NUPACK 3.0, the following changes were introduced to\n"
  "the complexes executable:\n"
  "  -ordered is on by default\n"
  "  output files .cx and .cx-epairs are not written\n"
  "Use the -v3.0 option to revert to NUPACK 3.0 behavior.\n\n";

  fprintf(out, "%s", dep_mess);
  
  return;
}



/* create a string containing the provenance's header informations:
 * - version
 * - command invocation
 * return the length of the generated string
 */
int complexes_header(char *provenance, int argc, char **argv) {

  char PROVENANCE_STARTS[] = "{ ";
  char FIELD_VERSION[] = "\"version\": \"";
  char FIELD_COMMAND[] = "\"command\": \"";
  char QUOTE[] = "\"";
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
        + strlen(QUOTE) + strlen(FIELD_NEXT);
  len_provenance += len_entry_version;

  // store the version
  provenance = strcat(
                strcat(
                  strcat(
                    strcat(provenance, FIELD_VERSION),
                    NUPACK_VERSION),
                  QUOTE),
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
        + strlen(QUOTE) + strlen(FIELD_NEXT);
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
                  QUOTE),
                FIELD_NEXT);


  return len_provenance;
}



/* create a string containing all provenance's parameter informations:
 * - dangles
 * - concentration sodium
 * - concentration magnesium
 * - temperature
 * - strands
 * return the length of the generated string
 */
int complexes_parameters(char* provenance, int nStrands, char **seqs,
        int nTotalOrders) {

  char PROVENANCE_ENDS[] = " }\n";
  char FIELD_TEMPERATURE[] = "\"temperature (C)\": ";
  char FIELD_DANGLES[] = "\"dangles\": \"";
  char FIELD_CONC_NA[] = "\"concentration Na (M)\": ";
  char FIELD_CONC_MG[] = "\"concentration Mg (M)\": ";
  char FIELD_STRANDS[] = "\"strands\": ";
  char FIELD_NEXT[] = ", ";
  char LIST_STARTS[] = "[";
  char LIST_ENDS[]   = "]";
  char PAIR_STARTS[] = "(";
  char PAIR_ENDS[]   = ")";
  char COMMA[] = ",";
  char QUOTE[] = "\"";


  int len_entry_temperature;
  int len_entry_dangles;
  int len_entry_conc_na;
  int len_entry_conc_mg;
  int len_entry_strands;
  int len_provenance = 0;


  /* retrieve each entry's value, and store it as provenance
   */

  provenance = strcpy(provenance, "");


  /* dangles
   */

  // retrieve the dangle's length, and add it to the provenance
  int len_nupack_dangles = (snprintf(NULL, 0, "%d", globalArgs.dangles) + 1);
  len_entry_dangles = strlen(FIELD_DANGLES) + len_nupack_dangles
        + strlen(QUOTE) + strlen(FIELD_NEXT);
  len_provenance += len_entry_dangles;

  // retrieve the dangle's value
  char nupack_dangles[len_nupack_dangles];
  snprintf(nupack_dangles, len_nupack_dangles, "%d", globalArgs.dangles);

  // store the dangles
  provenance = strcat(
                strcat(
                  strcat(
                    strcat(provenance, FIELD_DANGLES),
                    nupack_dangles),
                  QUOTE),
                FIELD_NEXT);


  /* temperature
   */

  // retrieve the temperature's length, and add it to the provenance
  int len_nupack_temperature = (snprintf(NULL, 0, "%.1f",
        ((float)globalArgs.T)) + 1);
  len_entry_temperature = strlen(FIELD_TEMPERATURE) + len_nupack_temperature
        + strlen(FIELD_NEXT);
  len_provenance += len_entry_temperature;

  // retrieve the temperature's value
  char nupack_temperature[len_nupack_temperature];
  snprintf(nupack_temperature, len_nupack_temperature, "%.1f",
        ((float)globalArgs.T));

  // store the temperature
  provenance = strcat(
                  strcat(
                    strcat(provenance, FIELD_TEMPERATURE),
                    nupack_temperature),
                  FIELD_NEXT);


  /* sodium concentration
   */

  // retrieve the sodium's concentration length, and add it to the
  // provenance
  int len_nupack_conc_na = (snprintf(NULL, 0, "%.4f",
        ((float)globalArgs.sodiumconc)) + 1);
  len_entry_conc_na = strlen(FIELD_CONC_NA) + len_nupack_conc_na
        + strlen(FIELD_NEXT);
  len_provenance += len_entry_conc_na;

  // retrieve the sodium's concentration value
  char nupack_conc_na[len_nupack_conc_na];
  snprintf(nupack_conc_na, len_nupack_conc_na, "%.4f",
        ((float)globalArgs.sodiumconc));

  // store the sodium concentration
  provenance = strcat(
                  strcat(
                    strcat(provenance, FIELD_CONC_NA),
                    nupack_conc_na),
                  FIELD_NEXT);


  /* magnesium concentration
   */

  // retrieve the magnesium's concentration length, and add it to the
  // provenance
  int len_nupack_conc_mg = (snprintf(NULL, 0, "%.4f",
        ((float)globalArgs.magnesiumconc)) + 1);
  len_entry_conc_mg = strlen(FIELD_CONC_MG) + len_nupack_conc_mg
        + strlen(FIELD_NEXT);
  len_provenance += len_entry_conc_mg;

  // retrieve the magnesium's concentration value
  char nupack_conc_mg[len_nupack_conc_mg];
  snprintf(nupack_conc_mg, len_nupack_conc_mg, "%.4f",
        ((float)globalArgs.magnesiumconc));

  // store the magnesium concentration
  provenance = strcat(
                  strcat(
                    strcat(provenance, FIELD_CONC_MG),
                    nupack_conc_mg),
                  FIELD_NEXT);


  /* strands
   */

  // retrieve the strands' length, and add it to the provenance
  int len_nupack_strands = strlen(LIST_STARTS);
  for(int j=0 ; j<nStrands ; ++j) {
    int len_nupack_strand_a = strlen(QUOTE) + (snprintf(NULL, 0, "%d", j) + 1)
        + strlen(QUOTE);
    int len_nupack_strand_b = strlen(QUOTE) + strlen(seqs[j]) + strlen(QUOTE);
    len_nupack_strands += (strlen(PAIR_STARTS)
        + len_nupack_strand_a
        + strlen(COMMA) +
        + len_nupack_strand_b
        + strlen(PAIR_ENDS));
    if(j<nStrands-1){
      len_nupack_strands += strlen(COMMA);
    }
  }
  len_nupack_strands += strlen(LIST_ENDS);

  len_entry_strands = (strlen(FIELD_STRANDS) + len_nupack_strands);

  len_provenance += len_entry_strands;

  // retrieve the strands' value
  char *nupack_strands = malloc(sizeof(char) * len_nupack_strands);
  if(!nupack_strands){
    exit(1);
  }
  for(size_t x=0 ; x < len_nupack_strands ; ++x){
    nupack_strands[x] = 0;
  }
  nupack_strands = strcpy(nupack_strands, LIST_STARTS);

  for(int j=0 ; j<nStrands ; ++j) {

    // strand's value a
    int len_nupack_strand_a = strlen(QUOTE) + (snprintf(NULL, 0, "%d", j) + 1)
        + strlen(QUOTE);
    char nupack_strand_a[len_nupack_strand_a];
    snprintf(nupack_strand_a, len_nupack_strand_a, "%s%d%s", QUOTE, j, QUOTE);

    // strand's value b
    int len_nupack_strand_b = strlen(QUOTE) + strlen(seqs[j]) + strlen(QUOTE);
    char *nupack_strand_b = malloc(sizeof(char) * len_nupack_strand_b);
    if(!nupack_strand_b){
      exit(1);
    }
    for(size_t x=0 ; x < len_nupack_strand_b ; ++x){
      nupack_strand_b[x] = 0;
    }
    nupack_strand_b = strcat(
                        strcat(
                          strcpy(nupack_strand_b, QUOTE),
                          seqs[j]),
                        QUOTE);

    int len_nupack_strand = (strlen(PAIR_STARTS)
            + len_nupack_strand_a + strlen(COMMA) + len_nupack_strand_b
            + strlen(PAIR_ENDS));
    if(j<nStrands-1){
      len_nupack_strand += strlen(COMMA);
    }

    char *nupack_strand = malloc(sizeof(char) * len_nupack_strand);
    if(!nupack_strand){
      exit(1);
    }
    for(size_t x=0 ; x < len_nupack_strand ; ++x){
      nupack_strand[x] = 0;
    }

    nupack_strand = strcat(
                      strcat(
                        strcat(
                          strcat(
                            strcpy(nupack_strand, PAIR_STARTS),
                            nupack_strand_a),
                          COMMA),
                        nupack_strand_b),
                      PAIR_ENDS);
    if(j<nStrands-1){
      nupack_strand = strcat(nupack_strand, COMMA);
    }

    nupack_strands = strcat(nupack_strands, nupack_strand);

    // free each strand
    free(nupack_strand);
    nupack_strand = NULL;

    free(nupack_strand_b);
    nupack_strand_b = NULL;
  }
  nupack_strands = strcat(nupack_strands, LIST_ENDS);

  // store the strands
  provenance = strcat(
                strcat(
                  strcat(provenance, FIELD_STRANDS),
                  nupack_strands),
                PROVENANCE_ENDS);


  // free all memory objects
  free(nupack_strands);
  nupack_strands = NULL;

  return len_provenance;
}


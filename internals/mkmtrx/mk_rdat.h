//
//  mk_rdat.h
//  
//
//  Created by Eric Strobel on 2/9/23.
//

#ifndef mk_rdat_h
#define mk_rdat_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "./cotrans_mtrx.h"
#include "./mkmtrx_defs.h"
#include "./mkmtrx_structs.h"

#define VL_RDAT 0          //TECprobe-VL RDAT mode
#define SL_RDAT 1          //TECprobe-SL RDAT mode
#define LM_RDAT 2          //TECprobe-LM RDAT mode

#define MAX_MD_STRING 512  //max metadata string length
#define MAX_OTHER_CHEM 32  //max number of other chemical strings
#define MAX_COMMENTS 32    //max number of user-supplied comments

#define RAW 0              //reactivities are raw background-subtracted mutation rates
#define NORMALIZED 1       //reactivities are normalized

typedef struct rdat_metadata {   //stores metadata for rdat generation
    char fn[MAX_MD_STRING+1];    //output rdat file name
    char nm[MAX_MD_STRING+1];    //sample name
    char str[MAX_MD_STRING+1];   //structure
    int  offset;                 //sequence position offset
    int  rtype;                  //reactivity type, can be RAW or NORMALIZED
    int  smooth;                 //flag to indicate if neighboring transcript smoothing was applied
    char xtype[MAX_MD_STRING+1]; //experiment type
    char xmeth[MAX_MD_STRING+1]; //experiment method
    char probe[MAX_MD_STRING+1]; //chemical probe
    char tris[MAX_MD_STRING+1];  //tris concentration
    char KCl[MAX_MD_STRING+1];   //KCl concentration
    char EDTA[MAX_MD_STRING+1];  //EDTA concentration
    char DTT[MAX_MD_STRING+1];   //DTT concentration
    char MgCl2[MAX_MD_STRING+1]; //MgCl2 concentration
    char NTPs[MAX_MD_STRING+1];  //NTP concentration
    char BSA[MAX_MD_STRING+1];   //BSA concentration
    char other_chem[MAX_OTHER_CHEM][MAX_MD_STRING+1]; //array of other chemicals in reaction
    char temp[MAX_MD_STRING+1];  //temperature
    char comment[MAX_COMMENTS][MAX_LINE]; //array to store user-supplied comments
    int oc_cnt;                  //count of other chemical entries provided
    int cmnt_cnt;                //count of user-supplied comments provided
} rdat_metadata;

/* mk_rdat: manages rdata config parsing/checking and
 construction of output rdat file*/
int mk_rdat(FILE * fp_config, cotrans_matrix * mtrx, int mode, int reps);

/* parde_rdat_config: read rdat config file and store values in rdat_metadata struct */
int parse_rdat_config(FILE * fp_config, rdat_metadata * rdat_meta);

/* check_rdat_config: check that required rdat metadata has been supplied
 and is consistent with the shapemapper2 output */
int check_rdat_config(rdat_metadata * rdat_meta, cotrans_matrix * mtrx);

/* print_rdat: print rdat file */
int print_rdat(rdat_metadata * rdat_meta, cotrans_matrix * mtrx, int mode, int reps);

/* set_TF_value: set true or false value as 1 and 0, respectively */
void set_TF_value(char * tf, char * setting_name, int * config_val);

#endif /* mk_rdat_h */

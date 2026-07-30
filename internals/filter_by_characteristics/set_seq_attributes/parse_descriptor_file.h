//
//  parse_descriptor_file.h
//  
//
//  Created by Eric Strobel on 11/12/25.
//

#ifndef parse_descriptor_file_h
#define parse_descriptor_file_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"
#include "../../utils/gen_utils.h"
#include "../../seq_utils/isIUPACbase.h"

#define DSCRPTR_TYPES 5 //number of descriptor types (including initial value)
#define TYPE_INIT 0     //initial descriptor type value
#define NUC_ID 1        //nucleotide identity type
#define PRX_DG 2        //proximal deltaG type
#define DST_DG 3        //distal deltaG type
#define SS_LEN 4        //subsequence length type

#define SEG1 0          //segment 1 index
#define SEG2 1          //segment 2 index

#define ACT_INIT 0      //initial action value
#define INCLUDE  1      //include action
#define EXCLUDE  2      //exclude action
#define SPLIT    3      //split action
#define BIN      4      //bin action

#define PTYP_INIT   0   //prediction type initiator
#define MFE_STRUCT  1   //predict minimum free energy structure
#define ALL_STRUCTS 2   //predict all structures

#define OP_INIT    0    //initial operator value
#define IS_MATCH   1    //is match operator
#define NOT_MATCH  2    //not match operator
#define STOCH_FLD  3    //stochastic bin operator
#define MFNRG_FLD  4    //minimum free energy bin operator
#define IS_EQUAL2  5    //equal to operator
#define NOT_EQL_2  6    //not equal to operator
#define GRTRorEQL  7    //greater than or equal to operator
#define LESSorEQL  8    //less than or equal to operator
#define GRTR_THAN  9    //greater than operator
#define LESS_THAN 10    //less than operator
#define EQL_WIDTH 11    //equal width binning
#define EQL_FRQNC 12    //equal frequency binning
#define HLX_PTTRN 13    //helix pattern binning

#define LONG_INTEGER 0  //value is a long int
#define FLOATING_PT  1  //value is a double
#define MIN_FREE_NRG 2  //run MFE structure prediction for helix pattern binning
#define STOCHASTIC   3  //run stochastic structure prediction for helix pattern binning
 
/* nt_window: structure for storing sequence window bounds */
//NOTE: for nucleotide identity descriptors, bound 1 is the index at
//which the string starts and bound 2 is the length of the string.
//for all other descriptors, bounds 1 and 2 are both indices.
typedef struct nt_window {
    int b1;
    int b2;
} nt_window;

/* des_action: structure for storing descriptor action parameters*/
typedef struct des_action {
    int act;
    int ptyp;
    int op[2];
    int op_cnt;
    int val_typ;
    int val_cnt;
    void * val;
} des_action;

/* descriptor: structure for storing descriptor parameters*/
typedef struct descriptor {
    int typ;
    char * nm;
    char * typ_s;
    char * act_s;
    char * val_s;
    char * sq;
    des_action act;
    int wndw_cnt;
    nt_window wndw[2];
    int tot_sp;
} descriptor;

/* parse_descriptor_file: parse descriptor file and store contents in descriptor structure */
int parse_descriptor_file(char * des_path, descriptor ** des, int * typ_cnt);

/* count_descriptors: count number of lines in descriptor file */
void count_descriptors(char * des_path, int * xpctd, int * typ_cnt);

/* store_field: store field from descriptor line and return pointer to next field */
char * store_field(char ** field, char * str);

/* ret_type: return descriptor type value */
int ret_type(char * str);

/* parse_descriptor action: parse descriptor file action string and store information in action structure */
void parse_descriptor_action(int typ, char * str, des_action * act);

/* parse_nuc_iden: parse nucleotide identity string */
void parse_nuc_iden(char * str, char ** sq, nt_window * wndw, int * wndw_cnt);

/* parse_window: parse window field in descriptor file line */
int parse_window(char * str, nt_window * wndw);

/* print_descriptor: print descriptor parameters */
void print_descriptor(descriptor * des);

#endif /* parse_descriptor_file_h */

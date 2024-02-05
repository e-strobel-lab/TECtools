//
//  mkmtrx_defs.h
//  
//
//  Created by Eric Strobel on 10/12/22.
//

#ifndef mkmtrx_defs_h
#define mkmtrx_defs_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "./cotrans_mtrx.h"
#include "./mkmtrx_defs.h"

#define NUCLEOTIDE 0          //position of nucleotide column in shapemapper output file
#define SEQUENCE 1            //position of sequence column in shapemapper output file
#define MOD_EFFECTIVE 4       //position of modified effective read depth column in shapemapper output file
#define UNT_EFFECTIVE 11      //position of untreated effective read depth column in shapemapper output file
#define REACTIVITY_PROFILE 23 //position of reactivity_profile column in shapemapper output file

#define MIN_INIT MAX_ROW+1    //initialization value for min_transcript variable
#define MAX_INIT 0            //initialization value for max_transcript variable

#define PRFL_SRCH_DEPTH 3              //expected directory hierarchy depth
#define PRFL_SUBDIRS PRFL_SRCH_DEPTH-1 //search depth at which shapemapper output file is stored

#define MULTI 0               //flag that data is from TECprobe-ML experiment
#define SINGLE 1              //flag that data is from TECprobe-SL experiment

#define MAX_MD_STRING 512     //max metadata string length
#define MAX_OTHER_CHEM 32     //max number of other chemical strings
#define MAX_COMMENTS 32       //max number of user-supplied comments

#endif /* mkmtrx_defs_h */

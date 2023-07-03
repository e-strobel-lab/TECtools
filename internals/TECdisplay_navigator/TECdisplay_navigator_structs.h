//
//  TECdisplay_navigator_structs.h
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#ifndef TECdisplay_navigator_structs_h
#define TECdisplay_navigator_structs_h

#include <stdio.h>

#include "../global/global_defs.h"

#include "./TECdisplay_navigator_defs.h"

typedef struct values_input {     //storage for file names and sample names of values files
    char nm[MAX_NAME+1];          //array to store sample name
    char fn[MAX_LINE+1];          //array to store file name
    FILE * fp;                    //pointer to values file
} values_input;

typedef struct base_params {      //storage for base parameters
    char seq;                     //base identity
    int pos;                      //base position
    int ins;                      //number of base in insertion contig
    int typ;                      //base type
} base_params;

typedef struct pair_params {      //storage for base pair parameters
    base_params bs[2];            //parameters of bases in pair
    int att;                      //pair attribute (e.g. WC_PAIR/STRONG/WEAK, see definitions)
} pair_params;

typedef struct constraints {      //storage for variant constraints
    char nm[MAX_NAME+1];          //name of constraints
    base_params base[MAX_VBASES]; //constrained bases
    pair_params pair[MAX_PAIRS];  //constrainted pairs
    int bcnt;                     //number of constrained bases
    int pcnt;                     //number of constrained pairs
} constraints;

//NOTE: this was previously used to store data too. this functionality
//is not currently being used, but I've left the structure in case that
//changes in the future
typedef struct sample_values {   //storage for sample data
    char nm[MAX_NAME+1];         //sample name
    //int bnd;                   //number of reads in bound channel
    //int unbnd;                 //number of reads in unbound channel
    //double frac_bnd;           //fraction of reads in bound channel
} sample_values;

#endif /* TECdisplay_navigator_structs_h */

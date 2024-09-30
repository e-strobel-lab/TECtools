//
//  TECdisplay_navigator_defs.h
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#ifndef TECdisplay_navigator_defs_h
#define TECdisplay_navigator_defs_h

#include <stdio.h>

#include "../global/global_defs.h"

#define MAX_VALS 32         //maximum number of input values files
#define MAX_CONSTRAINTS 32  //maximum number of constraints in constraints file
#define MAX_NAME 256        //maximum name length
#define MAX_CODE 8
#define MAX_VBASES 256      //maximum number of base constraints; larger than SEQ2BIN_MAX for constant indels storage
#define MAX_PAIRS 256       //maximum number of base pair constraints
#define VALUE_FIELDS 3      //number of data columns in values file
#define XPCTD_FIELDS 4      //expected number of fields in values file
#define MAX_COL_NM 128      //maximum column name length

//pair attributes
#define PAIR_TYPE_INIT -1   //initial pair type value
#define ANY_PAIR 0          //any base pair allowed
#define WC_PAIR 1           //Watson-Crick pair
#define STRONG 2            //strong pair
#define WEAK 3              //weak pair
#define WEAK_AT 4           //weak pair, AT/AU only
#define WEAK_AU 5           //weak pair, AT/AU only
#define WEAK_GT 6           //weak pair, GT/GU only
#define WEAK_GU 7           //weak pair, GT/GU only
#define MISMATCH 8          //mismatch
#define NO_CONSTRAINT 9     //no pairing constraint

#define MAX_POS_ARRAY 8     //maximum position string length in vbase string
#define POS1 0              //index of pair mate 1
#define POS2 1              //index of pair mate 2

#define BASE_CONSTRAINT 0   //read_vbase code for base constraint
#define PAIR_CONSTRAINT 1   //read_vbase code for pair constraint
#define CNST_CONSTRAINT 2   //read_vbase code for constant constraint
#define NAME 3              //read_vbase code for name

#define MAX_SAMPLE 16       //maximum number of samples
#define BOUND_VAL 0         //position of bound value field in data columns
#define UNBOUND_VAL 1       //position of unbound value field in data columns
#define FRAC_BND_VAL 2      //position of fraction bound value field in data columns

#define NEG1 -1             //negative 1

#endif /* TECdisplay_navigator_defs_h */

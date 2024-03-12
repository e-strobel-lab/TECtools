//
//  assemble_TECprobeLMP_data_defs.h
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#ifndef assemble_TECprobeLMP_data_defs_h
#define assemble_TECprobeLMP_data_defs_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "../mkmtrx/cotrans_mtrx.h"

#define MODE_INIT -1            //initial value of mode
#define REACTIVITY 0            //value for REACTIVITY mode
#define LENGTH_DISTRIBUTION 1   //value for LENGTH_DISTRIUTION mode

#define TOT_SAMPLES 3           //total sample number
#define MAX_IPT 16              //maximum number of input files per sample
#define MAX_TRANSCRIPT MAX_ROW  //maximum transcript length
 
#define S1 0 //sample 1 index
#define S2 1 //sample 2 index
#define S3 2 //sample 3 index

#endif /* assemble_TECprobeLMP_data_defs_h */

//
//  process_TECprobeVL_profiles_defs.h
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#ifndef process_TECprobeVL_profiles_defs_h
#define process_TECprobeVL_profiles_defs_h

#include <stdio.h>

#include "../../global/global_defs.h"
#include "../../mkmtrx/mkmtrx_defs.h"
#include "../global/store_SM2_profile.h"

#define MIN 0           //min value index
#define MAX 1           //max value index

#define COPY 0          //op code to copy field value to merged profile
#define RSUM 1          //op code to sum field values in merged profile
#define RATE 2          //op code to calculate mutation rate in merged profile
#define BSUB 3          //op code to calculate background-subtracted reactivity in merged profile
#define MASK 4          //op code to mask value in merged profile

#define SEQ_COL 1       //index of Nucleotide column

#define MOD_MUT 2       //index of Modified_mutations column
#define MOD_EFF_DEP 4   //index of Modified_effective_depth column
#define MOD_RATE 5      //index of Modified_rate column

#define UNT_MUT 9       //index of Untreated_mutations column
#define UNT_EFF_DEP 11  //index of Untreated_effective_depth column
#define UNT_RATE 12     //index of Untreated_rate column

#define DEN_MUT 16      //index of Denatured_mutations column
#define DEN_READ_DEP 17 //index of Denatured_read_depth column
#define DEN_EFF_DEP 18  //index of Denatured_effective_depth column
#define DEN_RATE 19     //index of Denatured_rate column

#define RCT_PRFL 23     //index of Reactivity_profile column

#endif /* process_TECprobeVL_profiles_defs_h */

//
//  print_input_filenames.h
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#ifndef print_input_filenames_h
#define print_input_filenames_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "./assemble_TECprobeLM_data_defs.h"
#include "./assemble_TECprobeLM_data_structs.h"

void print_input_filenames(char * out_dir, char * out_nm, char * ipt_nm[TOT_SAMPLES][MAX_IPT], int ipt_cnt[TOT_SAMPLES]);

#endif /* print_input_filenames_h */

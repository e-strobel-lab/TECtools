//
//  print_reactivity_output.h
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#ifndef print_reactivity_output_h
#define print_reactivity_output_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "./assemble_TECprobeLMP_data_defs.h"
#include "./assemble_TECprobeLMP_data_structs.h"

void print_reactivity_output(char * out_dir, char * out_nm, mode_parameters * mode_params, int ipt_cnt[TOT_SAMPLES], int nrchd_len[TOT_SAMPLES], double vals[MAX_TRANSCRIPT][TOT_SAMPLES][MAX_IPT]);

#endif /* print_reactivity_output_h */

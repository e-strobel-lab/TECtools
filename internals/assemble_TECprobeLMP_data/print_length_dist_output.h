//
//  print_length_dist_output.h
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#ifndef print_length_dist_output_h
#define print_length_dist_output_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "./assemble_TECprobeLMP_data_defs.h"
#include "./assemble_TECprobeLMP_data_structs.h"

void print_length_dist_output(char * out_dir, char * out_nm, mode_parameters * mode_params, int ipt_cnt[TOT_SAMPLES], int max_index, double vals[MAX_TRANSCRIPT][TOT_SAMPLES][MAX_IPT]);

#endif /* print_length_dist_output_h */

//
//  print_merged_VL_profiles.c
//  
//
//  Created by Eric Strobel on 2/21/24.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"

#include "../global/store_SM2_profile.h"

#include "print_merged_VL_profiles.h"

void print_merged_VL_profiles(SM2_analysis_directory * mrg, output_files * outfiles)
{
    char hdrs[PRFL_CLMNS][MAX_NAME] = {
        "Nucleotide",
        "Sequence",
        "Modified_mutations",
        "Modified_read_depth",
        "Modified_effective_depth",
        "Modified_rate",
        "Modified_off_target_mapped_depth",
        "Modified_low_mapq_mapped_depth",
        "Modified_mapped_depth",
        "Untreated_mutations",
        "Untreated_read_depth",
        "Untreated_effective_depth",
        "Untreated_rate",
        "Untreated_off_target_mapped_depth",
        "Untreated_low_mapq_mapped_depth",
        "Untreated_mapped_depth",
        "Denatured_mutations",
        "Denatured_read_depth",
        "Denatured_effective_depth",
        "Denatured_rate",
        "Denatured_off_target_mapped_depth",
        "Denatured_low_mapq_mapped_depth",
        "Denatured_mapped_depth",
        "Reactivity_profile",
        "Std_err",
        "HQ_profile",
        "HQ_stderr",
        "Norm_profile",
        "Norm_stderr"
    };
    
    int tl = 0;
    int col = 0;
    int i = 0;
    
    for (tl = mrg->min_tl; tl < mrg->max_tl; tl++) {
        for (col = 0; col < PRFL_CLMNS; col++) {
            if (!col) {
                fprintf(outfiles->ofp[tl], "%s", hdrs[col]);
            } else {
                fprintf(outfiles->ofp[tl], "\t%s", hdrs[col]);
            }
        }
        fprintf(outfiles->ofp[tl], "\n");
        
        //printf("%d\n", mrg->data[tl].tot_nt_cnt);
        for (i = 0; i < mrg->data[tl].tot_nt_cnt; i++) {
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].nucleotide[i]);
            fprintf(outfiles->ofp[tl], "%c\t", mrg->data[tl].sequence[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].mod_mutations[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].mod_read_depth[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].mod_eff_depth[i]);
            fprintf(outfiles->ofp[tl], "%.6f\t", mrg->data[tl].mod_rate[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].mod_off_target_depth[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].mod_low_mapq_depth[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].mod_mapped_depth[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].unt_mutations[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].unt_read_depth[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].unt_eff_depth[i]);
            fprintf(outfiles->ofp[tl], "%.6f\t", mrg->data[tl].unt_rate[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].unt_off_target_depth[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].unt_low_mapq_depth[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].unt_mapped_depth[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].den_mutations[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].den_read_depth[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].den_eff_depth[i]);
            fprintf(outfiles->ofp[tl], "%.6f\t", mrg->data[tl].den_rate[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].den_off_target_depth[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].den_low_mapq_depth[i]);
            fprintf(outfiles->ofp[tl], "%d\t", mrg->data[tl].den_mapped_depth[i]);
            fprintf(outfiles->ofp[tl], "%.6f\t", mrg->data[tl].reactivity_profile[i]);
            fprintf(outfiles->ofp[tl], "%.6f\t", mrg->data[tl].std_err[i]);
            fprintf(outfiles->ofp[tl], "%.6f\t", mrg->data[tl].hq_profile[i]);
            fprintf(outfiles->ofp[tl], "%.6f\t", mrg->data[tl].hq_stderr[i]);
            fprintf(outfiles->ofp[tl], "%.6f\t", mrg->data[tl].dataset_norm_profile[i]);
            fprintf(outfiles->ofp[tl], "%.6f", mrg->data[tl].norm_stderr[i]);
            fprintf(outfiles->ofp[tl],"%c", (i+1 == mrg->data[tl].tot_nt_cnt) ?  '\0' : '\n');
        }
    }
}

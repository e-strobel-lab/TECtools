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

/* print_merged_VL_profiles: generated merged shapemapper2 output files */
void print_merged_VL_profiles(SM2_analysis_directory * mrg, output_files * outfiles)
{
    //list of headers used in shapemapper2 output files
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
    
    int col = 0; //column index
    int i = 0;   //general purpose index
    int j = 0;   //general purpose index
    
    int * ix = &mrg->indx[0]; //set pointer to target indices
    
    for (i = 0; ix[i] <= mrg->max_id; i++) {           //for each target
        //print header line
        for (col = 0; col < PRFL_CLMNS; col++) {                //for each column
            if (!col) {                                         //if printing the first column
                fprintf(outfiles->ofp[ix[i]], "%s", hdrs[col]);    //omit leading tab
            } else {                                            //otherwise
                fprintf(outfiles->ofp[ix[i]], "\t%s", hdrs[col]);  //include leading tab
            }
        }
        fprintf(outfiles->ofp[ix[i]], "\n"); //print terminal newline
        
        //print data line for each nucleotide
        for (j = 0; j < mrg->data[ix[i]].tot_nt_cnt; j++) {
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].nucleotide[j]);
            fprintf(outfiles->ofp[ix[i]], "%c\t",   mrg->data[ix[i]].sequence[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].mod_mutations[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].mod_read_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].mod_eff_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%.6f\t", mrg->data[ix[i]].mod_rate[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].mod_off_target_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].mod_low_mapq_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].mod_mapped_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].unt_mutations[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].unt_read_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].unt_eff_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%.6f\t", mrg->data[ix[i]].unt_rate[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].unt_off_target_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].unt_low_mapq_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].unt_mapped_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].den_mutations[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].den_read_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].den_eff_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%.6f\t", mrg->data[ix[i]].den_rate[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].den_off_target_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].den_low_mapq_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%d\t",   mrg->data[ix[i]].den_mapped_depth[j]);
            fprintf(outfiles->ofp[ix[i]], "%.6f\t", mrg->data[ix[i]].reactivity_profile[j]);
            fprintf(outfiles->ofp[ix[i]], "%.6f\t", mrg->data[ix[i]].std_err[j]);
            fprintf(outfiles->ofp[ix[i]], "%.6f\t", mrg->data[ix[i]].hq_profile[j]);
            fprintf(outfiles->ofp[ix[i]], "%.6f\t", mrg->data[ix[i]].hq_stderr[j]);
            fprintf(outfiles->ofp[ix[i]], "%.6f\t", mrg->data[ix[i]].dataset_norm_profile[j]);
            fprintf(outfiles->ofp[ix[i]], "%.6f",   mrg->data[ix[i]].dataset_norm_stderr[j]);
            fprintf(outfiles->ofp[ix[i]],"%c", (j+1 == mrg->data[ix[i]].tot_nt_cnt) ?  '\0' : '\n');
        }
    }
}

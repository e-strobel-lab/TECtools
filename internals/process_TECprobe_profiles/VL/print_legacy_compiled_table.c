//
//  print_legacy_compiled_table.c
//  
//
//  Created by Eric Strobel on 2/24/24.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"

#include "../global/store_SM2_profile.h"

#include "print_legacy_compiled_table.h"

/* print_legacy_compiled_table: print aggregate TECprobe-VL profiles in the format used by Courtney's scripts */
void print_legacy_compiled_table(SM2_analysis_directory * an_dir, output_files * outfiles)
{
    
    char hdrs[PRFL_CLMNS+1][MAX_NAME] = {
        "Denatured_effective_depth",
        "Denatured_low_mapq_mapped_depth",
        "Denatured_mapped_depth",
        "Denatured_mutations",
        "Denatured_off_target_mapped_depth",
        "Denatured_rate",
        "Denatured_read_depth",
        "HQ_profile",
        "HQ_stderr",
        "Indiv_Norm_profile",
        "Indiv_Norm_stderr",
        "Modified_effective_depth",
        "Modified_low_mapq_mapped_depth",
        "Modified_mapped_depth",
        "Modified_mutations",
        "Modified_off_target_mapped_depth",
        "Modified_rate",
        "Modified_read_depth",
        "Norm_calc_profile",
        "Nucleotide",
        "Reactivity_profile",
        "Sequence",
        "Std_err",
        "Untreated_effective_depth",
        "Untreated_low_mapq_mapped_depth",
        "Untreated_mapped_depth",
        "Untreated_mutations",
        "Untreated_off_target_mapped_depth",
        "Untreated_rate",
        "Untreated_read_depth"
    };
    
    int col = 0;  //general purpose index
    int tl = 0;  //general purpose index
    int nt = 0; //nucleotide index
    
    FILE * p_lgcy_tbl = NULL;         //legacy table file pointer
    char lgcy_tbl_nm[MAX_LINE] = {0}; //legacy table record file name
    
    //generate legacy table file name
    if (snprintf(lgcy_tbl_nm, MAX_LINE, "./%s/000_full_table.csv", outfiles->out_dir) >= MAX_LINE) {
        printf("print_processing_record: error - legacy compiled table file name exceeded buffer. aborting...\n");
        abort();
    }
    
    //open legacy table file
    if ((p_lgcy_tbl = fopen(lgcy_tbl_nm, "w")) == NULL) {
        printf("print_processing_record: error - could not open legacy compiled table output file. Aborting program...\n");
        abort();
    }
    
    fprintf(p_lgcy_tbl, "Nucleotide_Position");
    for (col = 0; col < PRFL_CLMNS+1; col++) {
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            fprintf(p_lgcy_tbl, ",%s_%03d", hdrs[col], tl);
        }
    }
    fprintf(p_lgcy_tbl, "\n");
        
    for (nt = an_dir->trgt_start; nt < (an_dir->max_tl + an_dir->trgt_start); nt++) {
        fprintf(p_lgcy_tbl, "%d", nt+1-an_dir->trgt_start);
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].den_eff_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].den_low_mapq_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].den_mapped_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].den_mutations[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].den_off_target_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", an_dir->data[tl].den_rate[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].den_read_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", an_dir->data[tl].hq_profile[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", an_dir->data[tl].hq_stderr[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", an_dir->data[tl].norm_profile[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", an_dir->data[tl].norm_stderr[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].mod_eff_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].mod_low_mapq_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].mod_mapped_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].mod_mutations[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].mod_off_target_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", an_dir->data[tl].mod_rate[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].mod_read_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", an_dir->data[tl].dataset_norm_profile[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].nucleotide[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", an_dir->data[tl].reactivity_profile[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%c", an_dir->data[tl].sequence[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", an_dir->data[tl].std_err[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].unt_eff_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].unt_low_mapq_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].unt_mapped_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].unt_mutations[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].unt_off_target_depth[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", an_dir->data[tl].unt_rate[nt]);
            }
        }
        
        for (tl = an_dir->min_tl; tl <= an_dir->max_tl; tl++) {
            if ((nt + 1 - an_dir->trgt_start) > tl) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[tl].unt_read_depth[nt]);
            }
        }

        fprintf(p_lgcy_tbl, "\n");
    }
    
    //close legacy table record file
    if (fclose(p_lgcy_tbl) == EOF) {
        printf("print_processing_record: error - failed to close legacy compiled table file. Aborting program...\n");
        abort();
    }
}

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

#include "../process_TECprobe_profiles_defs.h"
#include "../process_TECprobe_profiles_structs.h"

#include "../UNV/store_SM2_profile.h"

#include "print_legacy_compiled_table.h"

/* print_legacy_compiled_table: print aggregate TECprobe profiles-VL in the format used by Courtney's visualization tools */
void print_legacy_compiled_table(SM2_analysis_directory * an_dir, char * out_dir, sample_names * sn)
{
    //headers used in legacy compiled table
    char hdrs[PRFL_CLMNS+2][MAX_NAME] = {
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
        "Norm_stderr",
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
    
    int i = 0;    //general purpose index
    int col = 0;  //general purpose index
    int nt = 0;   //nucleotide index
    
    FILE * p_lgcy_tbl = NULL;         //legacy table file pointer
    char lgcy_tbl_nm[MAX_LINE+1] = {0}; //legacy table record file name
    int ret = 0;  //variable for storing snprintf return value
    
    int * ix = &an_dir->indx[0]; //set pointer to target indices
    
    //generate legacy table file name
    ret = snprintf(lgcy_tbl_nm, MAX_LINE, "./%s/%s_full_table.csv", out_dir, sn->sn2use);
    if (ret >= MAX_LINE || ret < 0) {
        printf("print_processing_record: error - error when constructing legacy compiled table file name. aborting...\n");
        abort();
    }
    
    //open legacy table file
    if ((p_lgcy_tbl = fopen(lgcy_tbl_nm, "w")) == NULL) {
        printf("print_processing_record: error - could not open legacy compiled table output file. Aborting program...\n");
        abort();
    }
    
    //print column headers
    fprintf(p_lgcy_tbl, "Transcript Length");
    
    for (col = 0; col < PRFL_CLMNS+2; col++) {
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            fprintf(p_lgcy_tbl, ",%s_%03d", hdrs[col], ix[i]);
        }
    }
    fprintf(p_lgcy_tbl, "\n");
        
    //for each nt in the target RNA, print all shapemapper2 output file
    //values plus whole dataset normalized reactivity values in matrix
    //format. as in the original compiled shapemapper2 file format used
    //in Courtney's visualization tools, the values are in alphabetical
    //order and NAN values are masked as 0
    for (nt = an_dir->trgt_start; nt < (an_dir->len[MAX] + an_dir->trgt_start); nt++) {
        
        //print the current nucleotide number
        fprintf(p_lgcy_tbl, "%d", nt+1-an_dir->trgt_start);
        
        //print denatured effective read depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].den_eff_depth[nt]);
            }
        }
        
        //print denatured low mapq depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].den_low_mapq_depth[nt]);
            }
        }
        
        //print denatured mapped depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].den_mapped_depth[nt]);
            }
        }
        
        //print denatured mutation count
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].den_mutations[nt]);
            }
        }
        
        //print denatured off target depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].den_off_target_depth[nt]);
            }
        }
        
        //print denatured mutation rate
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", (!isnan(an_dir->data[ix[i]].den_rate[nt])) ? an_dir->data[ix[i]].den_rate[nt] : 0);
            }
        }
        
        //print denatured read depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].den_read_depth[nt]);
            }
        }
        
        //print hq reactivity profile
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", (!isnan(an_dir->data[ix[i]].hq_profile[nt])) ? an_dir->data[ix[i]].hq_profile[nt] : 0);
            }
        }
        
        //print hq profile standard error
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", (!isnan(an_dir->data[ix[i]].hq_stderr[nt])) ? an_dir->data[ix[i]].hq_stderr[nt] : 0);
            }
        }
        
        //print individual target normalized reactivity profile
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", (!isnan(an_dir->data[ix[i]].norm_profile[nt])) ? an_dir->data[ix[i]].norm_profile[nt] : 0);
            }
        }
        
        //print individual transcrip length profile standard error
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", (!isnan(an_dir->data[ix[i]].norm_stderr[nt])) ? an_dir->data[ix[i]].norm_stderr[nt] : 0);
            }
        }
        
        //print modified effective depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].mod_eff_depth[nt]);
            }
        }
        
        //print modified low mapq depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].mod_low_mapq_depth[nt]);
            }
        }
        
        //print modified mapped depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].mod_mapped_depth[nt]);
            }
        }
        
        //print modified mutation count
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].mod_mutations[nt]);
            }
        }
        
        //print modified off target depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].mod_off_target_depth[nt]);
            }
        }
        
        //print modified mutation rate
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", (!isnan(an_dir->data[ix[i]].mod_rate[nt])) ? an_dir->data[ix[i]].mod_rate[nt] : 0);
            }
        }
        
        //print modified read depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].mod_read_depth[nt]);
            }
        }
        
        //print dataset normalized reactivity profile
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", (!isnan(an_dir->data[ix[i]].dataset_norm_profile[nt])) ? an_dir->data[ix[i]].dataset_norm_profile[nt] : 0);
            }
        }
        
        //print dataset normalized profile standard error
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", (!isnan(an_dir->data[ix[i]].dataset_norm_stderr[nt])) ? an_dir->data[ix[i]].dataset_norm_stderr[nt] : 0);
            }
        }
        
        //print nucleotide number
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].nucleotide[nt]);
            }
        }
        
        //print (non-normalized) reactivity profile
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", (!isnan(an_dir->data[ix[i]].reactivity_profile[nt])) ? an_dir->data[ix[i]].reactivity_profile[nt] : 0);
            }
        }
        
        //print sequence
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%c", an_dir->data[ix[i]].sequence[nt]);
            }
        }
        
        //print (non-normalized) profile standard error
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", (!isnan(an_dir->data[ix[i]].std_err[nt])) ? an_dir->data[ix[i]].std_err[nt] : 0);
            }
        }
        
        //print untreated effective depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].unt_eff_depth[nt]);
            }
        }
        
        //print untreated low mapq depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].unt_low_mapq_depth[nt]);
            }
        }
        
        //print untreated mapped depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].unt_mapped_depth[nt]);
            }
        }
        
        //print untreated mutation count
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].unt_mutations[nt]);
            }
        }
        
        //print untreated off target depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].unt_off_target_depth[nt]);
            }
        }
        
        //print untreated mutation rate
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%.6f", (!isnan(an_dir->data[ix[i]].unt_rate[nt])) ? an_dir->data[ix[i]].unt_rate[nt] : 0);
            }
        }
        
        //print untreated read depth
        for (i = 0; ix[i] <= an_dir->max_id; i++) {
            if ((nt + 1 - an_dir->trgt_start) > an_dir->data[ix[i]].trg_nt_cnt) {
                fprintf(p_lgcy_tbl, ",");
            } else {
                fprintf(p_lgcy_tbl, ",%d", an_dir->data[ix[i]].unt_read_depth[nt]);
            }
        }

        fprintf(p_lgcy_tbl, "\n"); //terminate line with newline char
    }
    
    //close legacy table record file
    if (fclose(p_lgcy_tbl) == EOF) {
        printf("print_processing_record: error - failed to close legacy compiled table file. Aborting program...\n");
        abort();
    }
}

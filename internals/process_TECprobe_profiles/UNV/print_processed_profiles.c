//
//  print_processed_profiles.c
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

#include "../../global/global_defs.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "../process_TECprobe_profiles_defs.h"
#include "../process_TECprobe_profiles_structs.h"

#include "../UNV/store_SM2_profile.h"

#include "print_processed_profiles.h"

/* print_processed_profiles: generate output directories and files containing processed data*/
void print_processed_profiles(SM2_analysis_directory * an_dir, char * out_dir, sample_names * sn)
{
    FILE * ofp = NULL; //output file pointer
    
    int i = 0;   //general purpose index
    int ret = 0; //variable for storing snprintf return value
    
    char trg_dir_nm[MAX_NAME+1] = {'\0'}; //array to store target directory name
    char prf_dir_nm[MAX_NAME+1] = {'\0'}; //array to store reactivity profile directory name
    char profile_nm[MAX_NAME+1] = {'\0'}; //array to store reactivity profile file name
    
    mk_out_dir(out_dir); //make parent output directory
    
    int * ix = &an_dir->indx[0]; //set pointer to target indices
    
    for (i = 0; ix[i] <= an_dir->max_id; i++) { //for every target
        
        //generate an analysis directory with the format:
        //<three-digit target id>_analysis
        ret = snprintf(trg_dir_nm, MAX_NAME, "%s/%03d_analysis", out_dir, ix[i]);
        if (ret >= MAX_NAME || ret < 0) {
            printf("print_processed_profiles: error - error when constructing target analysis directory name. aborting...\n");
            abort();
        }
        mk_out_dir(trg_dir_nm);
        
        //generate shapemapper 2 output directories with the format:
        //<sample_name>_<three-digit target id>_out
        ret = snprintf(prf_dir_nm, MAX_NAME, "%s/%s_%03d_out", trg_dir_nm, sn->sn2use, ix[i]);
        if (ret >= MAX_NAME || ret < 0) {
            printf("print_processed_profiles: error - error when constructing SM2 output directory name. aborting...\n");
            abort();
        }
        mk_out_dir(prf_dir_nm);
        
        //generate output profile files with the format:
        //<sample_name>_<three-digit target id>_nt_profile.txt
        ret = snprintf(profile_nm, MAX_NAME, "%s/%s_%03d_nt_profile.txt", prf_dir_nm, sn->sn2use, ix[i]);
        if (ret >= MAX_NAME || ret < 0) {
            printf("print_processed_profiles: error - error when constructing output file name. aborting...\n");
            abort();
        }
                
        if ((ofp = fopen(profile_nm, "w")) == NULL) {
            printf("print_processed_profiles: error - could not open output file %s. Aborting program...\n", profile_nm);
            abort();
        }
        
        print_profile(&an_dir->data[ix[i]], ofp);
        
        if (fclose(ofp) == EOF) {
            printf("print_processed_profiles: error - could not close output file %s. aborting program...\n", profile_nm);
            abort();
        }
        ofp = NULL;
    }
}


/* print_profile: print shapemapper2 profile to file */
void print_profile(SM2_profile * prf, FILE * ofp)
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
    
    //print header line
    for (col = 0; col < PRFL_CLMNS; col++) {  //for each column
        if (!col) {                           //if printing the first column
            fprintf(ofp, "%s", hdrs[col]);    //omit leading tab
        } else {                              //otherwise
            fprintf(ofp, "\t%s", hdrs[col]);  //include leading tab
        }
    }
    fprintf(ofp, "\n"); //print terminal newline
    
    //print data line for each nucleotide
    for (i = 0; i < prf->tot_nt_cnt; i++) {
        fprintf(ofp, "%d\t",   prf->nucleotide[i]);
        fprintf(ofp, "%c\t",   prf->sequence[i]);
        fprintf(ofp, "%d\t",   prf->mod_mutations[i]);
        fprintf(ofp, "%d\t",   prf->mod_read_depth[i]);
        fprintf(ofp, "%d\t",   prf->mod_eff_depth[i]);
        fprintf(ofp, "%.6f\t", prf->mod_rate[i]);
        fprintf(ofp, "%d\t",   prf->mod_off_target_depth[i]);
        fprintf(ofp, "%d\t",   prf->mod_low_mapq_depth[i]);
        fprintf(ofp, "%d\t",   prf->mod_mapped_depth[i]);
        fprintf(ofp, "%d\t",   prf->unt_mutations[i]);
        fprintf(ofp, "%d\t",   prf->unt_read_depth[i]);
        fprintf(ofp, "%d\t",   prf->unt_eff_depth[i]);
        fprintf(ofp, "%.6f\t", prf->unt_rate[i]);
        fprintf(ofp, "%d\t",   prf->unt_off_target_depth[i]);
        fprintf(ofp, "%d\t",   prf->unt_low_mapq_depth[i]);
        fprintf(ofp, "%d\t",   prf->unt_mapped_depth[i]);
        fprintf(ofp, "%d\t",   prf->den_mutations[i]);
        fprintf(ofp, "%d\t",   prf->den_read_depth[i]);
        fprintf(ofp, "%d\t",   prf->den_eff_depth[i]);
        fprintf(ofp, "%.6f\t", prf->den_rate[i]);
        fprintf(ofp, "%d\t",   prf->den_off_target_depth[i]);
        fprintf(ofp, "%d\t",   prf->den_low_mapq_depth[i]);
        fprintf(ofp, "%d\t",   prf->den_mapped_depth[i]);
        fprintf(ofp, "%.6f\t", prf->reactivity_profile[i]);
        fprintf(ofp, "%.6f\t", prf->std_err[i]);
        fprintf(ofp, "%.6f\t", prf->hq_profile[i]);
        fprintf(ofp, "%.6f\t", prf->hq_stderr[i]);
        fprintf(ofp, "%.6f\t", prf->dataset_norm_profile[i]);
        fprintf(ofp, "%.6f",   prf->dataset_norm_stderr[i]);
        fprintf(ofp, "%c", (i+1 == prf->tot_nt_cnt) ?  '\0' : '\n');
    }
}

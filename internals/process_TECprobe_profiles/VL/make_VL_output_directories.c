//
//  make_VL_output_directories.c
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

#include "../../global/global_defs.h"
#include "../../mkmtrx/mkmtrx_defs.h"
#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"

#include "make_VL_output_directories.h"

/* make_VL_output_directories: generate output directories and files */
void make_VL_output_directories(SM2_analysis_directory * an_dir, output_files * outfiles, sample_names * sn)
{
    int i = 0; //general purpose index
    
    char tl_dir_nm[MAX_NAME+1] = {'\0'};  //array to store transcrip length directory name
    char prf_dir_nm[MAX_NAME+1] = {'\0'}; //array to store reactivity profile directory name
    char profile_nm[MAX_NAME+1] = {'\0'}; //array to store reactivity profile file name
    
    mk_out_dir(outfiles->out_dir); //make parent output directory
    
    for (i = an_dir[0].min_tl; i <= an_dir[0].max_tl; i++) { //for eavery transcrip length
        
        //generate a transcript length analysis directory with the format:
        //<three-digit transcript length>_analysis
        if (snprintf(tl_dir_nm, MAX_NAME, "%s/%03d_analysis", outfiles->out_dir, i) >= MAX_NAME) {
            printf("merged SM2_profiles: error - transcript length analysis directory name exceeded buffer. aborting...\n");
            abort();
        }
        mk_out_dir(tl_dir_nm);
        
        //generate shapemapper 2 output directories with the format:
        //<sample_name>_<three-digit transcript length>_out
        if (snprintf(prf_dir_nm, MAX_NAME, "%s/%s_%03d_out", tl_dir_nm, sn->sn2use, i) >= MAX_NAME) {
            printf("merged SM2_profiles: error - SM2 output directory name exceeded buffer. aborting...\n");
            abort();
        }
        mk_out_dir(prf_dir_nm);
        
        //generate output profile files with the format:
        //<sample_name>_<three-digit transcript length>_nt_profile.txt
        if (snprintf(outfiles->ofn[i], MAX_NAME, "%s/%s_%03d_nt_profile.txt", prf_dir_nm, sn->sn2use, i) >= MAX_NAME) {
            printf("merged SM2_profiles: error - output file name exceeded buffer. aborting...\n");
            abort();
        }
                
        if ((outfiles->ofp[i] = fopen(outfiles->ofn[i], "w")) == NULL) {
            printf("read_SM2out_directory: error - could not open output file %s. Aborting program...\n", profile_nm);
            abort();
        }
    }
}

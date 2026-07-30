//
//  run_RNAStructure.c
//  
//
//  Created by Eric Strobel on 2/27/26.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../seq_utils/mk_fasta.h"

#include "./parse_descriptor_file.h"

#include "run_RNAStructure.h"

/* run_RNAStructure_Fold: predict RNA structures using the RNAStructure Fold algorithm */
int run_RNAStructure_Fold(char * path2RNAStructure, char * nm, char * in_dir, char * out_dir, char * fa_sffx, char * path2ct, int ct_path_maxlen, int ptyp)
{
    char command[MAX_LINE+1] = {0}; //RNAStructure run command
    char path2fa[MAX_LINE+1] = {0}; //path to input fasta file
    int ret = 0;                    //return value
    
    //generate input fasta file path
    ret = snprintf(path2fa, MAX_LINE, "%s%s%s", in_dir, nm, fa_sffx);
    if (ret >= MAX_LINE || ret < 0) {
        printf("run_RNAStructure_Fold: error - failed to generate fasta file path. aborting...\n");
        abort();
    }
    
    //generate output ct file path
    ret = snprintf(path2ct, ct_path_maxlen, "%s%s.ct", out_dir, nm);
    if (ret >= MAX_LINE || ret < 0) {
        printf("run_RNAStructure_Fold: error - failed to generate ct file path. aborting...\n");
        abort();
    }
    
    char * ptyp_str = NULL;     //pointer to string used to specify prediction type
    char all_str[2] = " ";      //when keeping non-mfe structures that approximate the mfe, no option is specified
    char mfe_str[7] = " -mfe "; //when keeping only the mfe structure, the '-mfe' option is specified
    
    if (ptyp == MFE_STRUCT) {         //if making mfe structure prediction
        ptyp_str = &mfe_str[0];       //set prediction type pointer to mfe option string
        
    } else if (ptyp == ALL_STRUCTS) { //if keeping all predicted structures that approximate the mfe
        ptyp_str = &all_str[0];       //set prediction type pointer to string that only contains a space char
    }
    
    //construct RNA structure run command
    ret = snprintf(command, MAX_LINE, "%s/exe/Fold %s %s%s> %sout.tmp", path2RNAStructure, path2fa, path2ct, ptyp_str, out_dir);
    if (ret >= MAX_LINE || ret < 0) {
        printf("run_RNAStructure_Fold: error - failed to generate command. aborting...\n");
        abort();
    }
    system(command); //run RNAStructure
    
    return 1;
}

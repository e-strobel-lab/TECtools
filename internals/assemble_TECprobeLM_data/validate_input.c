//
//  validate_input.c
//  
//
//  Created by Eric Strobel on 2/8/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "../global/global_defs.h"
#include "../mkmtrx/cotrans_mtrx.h"
#include "../utils/io_management.h"

#include "../process_TECprobe_profiles/VL/parse_VL_sample_name.h"

#include "./assemble_TECprobeLM_data_defs.h"
#include "./assemble_TECprobeLM_data_structs.h"

#include "validate_input.h"

/* validate_input: perform basic checks to assess input sample compatiblity */
void validate_input (input_data * ipt, int mode)
{
    //check that number of input files for each sample is the same
    if (ipt->cnt[S2] != ipt->cnt[S1] || ipt->cnt[S3] != ipt->cnt[S1]) {
        printf("validate_input: error - input file counts not equal. aborting...\n");
        abort();
    }
    
    int sm = 0;                 //sample index
    int i = 0;                  //general purpose index
    int ipt_cnt = ipt->cnt[S1]; //input count
    
    char * sffx2rmv = NULL; //pointer to sffx 2 remove, depending on mode
    char tbl_sffx[12]  = {"_full_table"};     //suffix to remove from sample name
    char algn_sffx[17] = {"_alignment_rates"};//suffix to remove from sample name
    
    if (mode == REACTIVITY) {      //if running REACTIVITY mode
        sffx2rmv = &tbl_sffx[0];   //remove full table suffix
    } else if (mode == LEN_DIST) { //if running LEN_DIST mode
        sffx2rmv = &algn_sffx[0];  //remove the alignment rate suffix
    } else {
        printf("validate_input: unexpected mode. aborting...\n");
        abort();
    }
    
    for (sm = 0; sm < TOT_SAMPLES; sm++) { //for each sample
        
        for (i = 0; i < ipt_cnt; i++) {    //for each sub-sample
            
            //allocate memory for sample name
            if ((ipt->sn[sm][i] = malloc((strlen(ipt->fn[sm][i])+1) * sizeof(*(ipt->sn[sm][i])))) == NULL) {
                printf("store_ipt_name: error - memory allocation for input name storage failed. aborting...\n");
                abort();
            }
            
            //get parsable sample name
            get_sample_name(ipt->fn[sm][i], ipt->sn[sm][i]);
            remove_simple_suffix(ipt->sn[sm][i], sffx2rmv);
            
            //parse sample name, then compare the attributes of current
            //sub-sample name to the attributes of the first sub-sample name
            parse_VL_sample_name(&ipt->sn[sm][i][0], &ipt->cfg[sm][i]);
            validate_sample_compatibility(&ipt->cfg[sm][i], &ipt->cfg[sm][0]);
        }
        
        //compare the attributes of the current sample name
        //to the attributes of the first sample name
        validate_sample_compatibility(&ipt->cfg[sm][0], &ipt->cfg[0][0]);
    }
}


/* validate_sample_compatibility: compare sample attributes to assess whether samples are compatible */
void validate_sample_compatibility(configuration_MLT * cfgT, configuration_MLT * cfgR)
{
    //check that RNA names are identical
    if (strcmp(cfgT->input_name, cfgR->input_name)) {
        printf("validate_sample_compatibility: error - samples contain discordant input names. aborting...\n");
        abort();
    }
    
    //check that folding type is identical
    if (cfgT->cotranscriptional != cfgR->cotranscriptional) {
        printf("validate_sample_compatibility: error - samples contain discordant folding type. aborting...\n");
        abort();
    }
    
    //check that chemical probe is identical
    if (strcmp(cfgT->chemical_probe, cfgR->chemical_probe)) {
        printf("validate_sample_compatibility: error - samples were not probed using the same chemical probe. aborting...\n");
        abort();
    }
    
    //check that the data are from one sequencing run
    if (cfgR->run_count != 1) {
        printf("validate_sample_compatibility: error - expected input data to be individual runs, but detected merged/concatenated data. aborting...\n");
        abort();
    }
    
    //check that the smoothing flag is identical
    if (cfgT->smoothing != cfgR->smoothing) {
        printf("validate_sample_compatibility: error - application of neighboring transcript smoothing is not uniform across all samples. aborting...\n");
        abort();
    }

    //if either sample contained a ligand name, check that
    //ligand names and concentrations are identical
    if (cfgR->ligand_name[0] || cfgT->ligand_name[0]) {
        if (strcmp(cfgT->ligand_name, cfgR->ligand_name)) {
            printf("validate_sample_compatibility: error - samples contain discordant ligand names. aborting...\n");
            abort();
        }
        
        if (strcmp(cfgT->ligand_conc, cfgR->ligand_conc)) {
            printf("validate_sample_compatibility: error - samples contain discordant ligand concentrations. aborting...\n");
            abort();
        }
    }
}

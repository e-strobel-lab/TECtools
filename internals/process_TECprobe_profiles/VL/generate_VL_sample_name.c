//
//  generate_VL_sample_name.c
//  
//
//  Created by Eric Strobel on 1/26/24.
//

#include <stdio.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../utils/io_management.h"

#include "../../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"
#include "../../cotrans_preprocessor/run_script_gen/MLT/mk_MLT_run_nm.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"

#include "./parse_VL_sample_name.h"

#include "generate_VL_sample_name.h"

/* generate_VL_sample_name: manages sample name parsing and generation */
void generate_VL_sample_name (sample_names * sn)
{
    int i = 0;              //general purpose index
    char * first_sn = NULL; //pointer to first sample name
    
    //parse each sample name
    for (i = 0; i < sn->cnt; i++) {     //for each sample name

        if (!i) {                       //if reading the first sample name
            first_sn = &sn->ipt[i][0];  //set pointer to the first sample name
            
        } else if (!strcmp(sn->ipt[i], first_sn)) { //otherwise, compare to the first sample name
            printf("generate_VL_sample_name: error - detected duplicate input for sample name %s. aborting...", first_sn);
            abort();
        }
        
        parse_VL_sample_name(sn->ipt[i], &sn->cfg[i]); //parse sample name
        
        if (sn->cfg[i].run_count != 1) {
            printf("generate_VL_sample_name: error - expected directory to contain data that was not concatenated prior to analysis. aborting...\n");
            abort();
        }
    }
    
    //merge the parsed sample names
    merge_sample_names(sn);
}

/* merge_sample_names: confirm that sample name attributes match and generate merged sample name */
void merge_sample_names(sample_names *sn)
{
    int i = 0; //general purpose index
    
    //set merged sample name config attributes. in most cases,
    //attributes are set using the config from the first input
    //directory. comparisons with configs from all other input
    //directories are then performed below
    
    set_cfg_string(&sn->mrgd_cfg.input_name, sn->cfg[0].input_name, 0);   //set RNA name
    sn->mrgd_cfg.cotranscriptional = sn->cfg[0].cotranscriptional;        //set folding type
    set_cfg_string(&sn->mrgd_cfg.chemical_probe, sn->cfg[0].chemical_probe, 0); //set chemical probe
    sn->mrgd_cfg.concatenated = (sn->cnt > 1) ? 1 : 0;                    //set concatenated flag
    sn->mrgd_cfg.run_count = sn->cfg[0].run_count;                        //set intial run count val
    set_cfg_string(&sn->mrgd_cfg.runID[0], sn->cfg[0].runID[0], 0);       //set first run ID
    sn->mrgd_cfg.smoothing = sn->cfg[0].smoothing;                        //set smoothing flag
    
    if (sn->cfg[0].ligand_name != NULL && sn->cfg[0].ligand_conc != NULL) {
        set_cfg_string(&sn->mrgd_cfg.ligand_name, sn->cfg[0].ligand_name, 0); //set ligand name
        set_cfg_string(&sn->mrgd_cfg.ligand_conc, sn->cfg[0].ligand_conc, 0); //set ligand concentration
    }
    
    //compare attributes of the first input directory (now in the merged config)
    //to attributes of all other input directories. during this process, the
    //run ID of each additional directory is added to the merged config and
    //the run_count variable with the config struct is incremented
    
    for (i = 1; i < sn->cnt; i++) { //for each input directory
        
        //check that RNA names are identical
        if (strcmp(sn->cfg[i].input_name, sn->mrgd_cfg.input_name)) {
            printf("merge_sample_names: error - samples contain discordant input names. aborting...\n");
            abort();
        }
        
        //check that folding type is identical
        if (sn->cfg[i].cotranscriptional != sn->mrgd_cfg.cotranscriptional) {
            printf("merge_sample_names: error - samples contain discordant folding type. aborting...\n");
            abort();
        }
        
        //check that chemical probe is identical
        if (strcmp(sn->cfg[i].chemical_probe, sn->mrgd_cfg.chemical_probe)) {
            printf("merge_sample_names: error - samples did not use the same chemical probe. aborting...\n");
            abort();
        }
        
        //if MAX_RUNS hasn't been reached, copy the run ID of the current input
        //directory to the merged config and increment the run_count variable
        if (sn->mrgd_cfg.run_count == MAX_RUNS) {
            printf("merge_sample_names: error - input directory count exceeds the maximum allowed number (%d). aborting...\n", MAX_RUNS);
            abort();
        } else {
            set_cfg_string(&sn->mrgd_cfg.runID[sn->mrgd_cfg.run_count++], sn->cfg[i].runID[0], 0);
        }
        
        //check that the smoothing flag is identical
        if (sn->cfg[i].smoothing != sn->mrgd_cfg.smoothing) {
            printf("merge_sample_names: error - application of neighboring transcript smoothing is not uniform across all samples. aborting...\n");
            abort();
        }

        //if either sample contained a ligand name, check that
        //ligand names and concentrations are identical
        if (sn->mrgd_cfg.ligand_name != NULL || sn->cfg[i].ligand_name != NULL) {
            if (strcmp(sn->cfg[i].ligand_name, sn->mrgd_cfg.ligand_name)) {
                printf("merge_sample_names: error - samples contain discordant ligand names. aborting...\n");
                abort();
            }
            
            if (strcmp(sn->cfg[i].ligand_conc, sn->mrgd_cfg.ligand_conc)) {
                printf("merge_sample_names: error - samples contain discordant ligand concentrations. aborting...\n");
                abort();
            }
        }
    }
    
    //TODO: this temporary code preserves variable fields in the merged sample name if only one input was provided. need to decide how to handle variable fields when multiple inputs are provided.
    if (sn->cnt == 1) {
        sn->mrgd_cfg.field_count = sn->cfg[0].field_count;
        for (i = 0; i < sn->cfg[0].field_count; i++) {
            set_cfg_string(&sn->mrgd_cfg.field[i], sn->cfg[0].field[i], 0);
        }
    }
    
    mk_MLT_run_nm(&sn->mrg[0], &sn->mrgd_cfg); //generate a run name using the merged config
}

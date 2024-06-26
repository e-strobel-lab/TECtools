//
//  mk_MLT_run_nm.c
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "config_MLT_struct.h"

#include "mk_MLT_run_nm.h"

/* mk_MLT_run_nm: generate sample name from config file settings */
int mk_MLT_run_nm(char *sample_name, configuration_MLT * config_MLT)
{
    int i = 0;
    
    //the order of settings in sample_name is:
    //1. input name				(required)
    //2. cotranscriptional flag	(required)
    //3. chemical probe			(required)
    //4. ligand name			(optional)
    //5. ligand concentration	(optional)
    //6. concatenated flag 		(required, but no flag is added to sample_name if concatenated=FALSE)
    //7. run identifier			(required, there may be more than one identifier)
    //8. custom fields			(optional)
    //9. smoothing flag			(required, but no flag is added to sample_name if smoothing=FALSE)
    
    //1.  copy input_name to sample array
    strcpy(sample_name, config_MLT->input_name);
    strcat(sample_name, "_");
    
    //2.  append 'CoTxn' for cotranscriptionally probed RNA or
    //    append 'Equil' for RNA that was probed at equilibrium
    if (config_MLT->cotranscriptional) {
        strcat(sample_name, "CoTxn_");
    } else {
        strcat(sample_name, "Equil_");
    }
    
    //3.  append chemical probe identifier
    strcat(sample_name, config_MLT->chemical_probe);
    strcat(sample_name, "_");
    
    //if ligand name and conc were provided
    //4.  append ligand name
    //5.  append ligand concentration
    if (config_MLT->ligand_name != NULL && config_MLT->ligand_conc != NULL) {
        strcat(sample_name, config_MLT->ligand_name);
        strcat(sample_name, "_");
        strcat(sample_name, config_MLT->ligand_conc);
        strcat(sample_name, "_");
    }
    
    //6.  append 'CAT<run count>_' if reads were concatenated from multiple sequencing runs
    //    max run count is 8, so <run count> will always be a single digit >=2
    char cat_str[64] = {0};
    
    if (config_MLT->concatenated) {
        sprintf(cat_str, "CAT%d_", config_MLT->run_count);
        strcat(sample_name, cat_str);
    }
    
    //7a. append runID; there will always be at least 1 runID
    strcat(sample_name, config_MLT->runID[0]);
    
    //7b. append any additional runIDs that were provided
    for (i = 1; i < config_MLT->run_count; i++) {
        strcat(sample_name, "_");
        strcat(sample_name, config_MLT->runID[i]);
    }
    
    //8.  append custom fields
    for (i = 0; i < config_MLT->field_count; i++) {
        strcat(sample_name, "_");
        strcat(sample_name, config_MLT->field[i]);
    }
    
    //9.  append 'SM' if analysis will be performed on smoothed fastqs
    if (config_MLT->smoothing) {
        strcat(sample_name, "_SM");
    }
    
    return 1;
}

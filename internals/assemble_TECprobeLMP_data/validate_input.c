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

#include "./assemble_TECprobeLMP_data_defs.h"
#include "./assemble_TECprobeLMP_data_structs.h"

#include "validate_input.h"

void validate_input (input_data * ipt)
{
    //check that number of input files for each sample is the same
    if (ipt->cnt[S2] != ipt->cnt[S1] || ipt->cnt[S3] != ipt->cnt[S1]) {
        printf("validate_input: error - input file counts not equal. aborting...\n");
        abort();
    }
    
    int sm = 0;                 //sample index
    int i = 0;                  //general purpose index
    int ipt_cnt = ipt->cnt[S1]; //input count
    
    //char tmp_sn[MAX_NAME] = {0}; //temp sample name //TODO: commented out to temporarily silence warning while working on other things
    
    char tbl_sffx[22] = {"_Normalized_Full_Table"};
    char * tbl_sffx_ptr = NULL;
    
    for (sm = 0; sm < TOT_SAMPLES; sm++) {
        for (i = 0; i < ipt_cnt; i++) {
            //allocate memory for sample name
            if ((ipt->sn[sm][i] = malloc((strlen(ipt->fn[sm][i])+1) * sizeof(*(ipt->sn[sm][i])))) == NULL) {
                printf("store_ipt_name: error - memory allocation for input name storage failed. aborting...\n");
                abort();
            }
            
            get_sample_name(ipt->fn[sm][i], ipt->sn[sm][i]);
            //printf("%s\n", ipt->sn[sm][i]);
            
            tbl_sffx_ptr = strstr(ipt->sn[sm][i], tbl_sffx);
            //printf("%d %d %s\n", sm, i, tbl_sffx_ptr);
            
            //printf("%llu\n", (uint64_t)(tbl_sffx_ptr) - (uint64_t)(ipt->sn[sm][i]));
            ipt->sn[sm][i][(uint64_t)(tbl_sffx_ptr) - (uint64_t)(ipt->sn[sm][i])] = '\0';
            printf("%s\n", ipt->sn[sm][i]);
            
            
            parse_VL_sample_name(&ipt->sn[sm][i][0], &ipt->cfg[sm][i]);
        }
        printf("\n");
    }
}

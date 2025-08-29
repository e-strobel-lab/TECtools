//
//  mk_smooth_script.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

#include "mk_smooth_script.h"

/* mk_smooth_script: generate shell script for concatenating 3' end split files
 the minimum target is concatenated as: n, n+1.
 the maximum target is concatenated as: n-1, n.
 intermediate targets are concatenated as: n-1, n, n+1
 */
void mk_smooth_script(TPROBE_names * nm, target3p_params trg_prms)
{
    FILE * out_fp = NULL;	//output file pointer
    
    int i = 0;
    int j = 0;
    int k = 0;
    
    static const char chnl_code[2][4] = {"UNT", "MOD"};	//array for channel codes, no need for ERR code
    static const char read_code[READ_MAX][3] = {"R1", "R2"}; //array for read codes
    
    //generate smooth_transition.sh shell script file
    if ((out_fp = fopen("./split/smooth_transition.sh", "w")) == NULL) {
        printf("mk_smooth_script: ERROR - could not generate shell script. Aborting program...\n");
        abort();
    }
    fprintf(out_fp, "#!/bin/bash\n"); //print she-bang
    
    //print command to make ./split_smooth directory,
    //which will store smoothed fastq files
    fprintf(out_fp, "mkdir ./split_smooth\n");
    
    //generate commands for concatenating minimum and maximum targets
    for (i = 0; i < 2; i++) { 			 //i loop iterates channel codes (UNT then MOD)
        for (j = 0; j < READ_MAX; j++) { //j loop iterates read codes    (R1 then R2)
            
            //the minimum target is concatenated as: min, min+1 > min_SM
            fprintf(out_fp, "cat ./split/%s_%s_%03d_%s.fq.gz ./split/%s_%s_%03d_%s.fq.gz > ./split_smooth/%s_%s_%03d_%s_SM.fq.gz\n",
                    nm->smpl[j], chnl_code[i], trg_prms.min,   read_code[j],  /* input  length min    */
                    nm->smpl[j], chnl_code[i], trg_prms.min+1, read_code[j],  /* input  length min+1  */
                    nm->smpl[j], chnl_code[i], trg_prms.min,   read_code[j]); /* output length min_SM */
            
            //the maximum target is concatenated as: max-1, max > max_SM
            fprintf(out_fp, "cat ./split/%s_%s_%03d_%s.fq.gz ./split/%s_%s_%03d_%s.fq.gz > ./split_smooth/%s_%s_%03d_%s_SM.fq.gz\n",
                    nm->smpl[j], chnl_code[i], trg_prms.max-1, read_code[j],  /* input  length max-1  */
                    nm->smpl[j], chnl_code[i], trg_prms.max,   read_code[j],  /* input  length max    */
                    nm->smpl[j], chnl_code[i], trg_prms.max,   read_code[j]); /* output length max_SM */
        }
    }
    
    //generate commands for concatenating intermediate targets
    for (i = trg_prms.min+1; i <= (trg_prms.max-1); i++) {	//i loop iterates 3' end length (min+1 to max-1)
        for (j = 0; j < 2; j++) {							//j loop iterates channel codes (UNT then MOD)
            for (k = 0; k < READ_MAX; k++) {				//k loop iterates read codes    (R1 then R2)
                
                //internal targets are concatenated as: n-1, n, n+1 > n_SM
                fprintf(out_fp, "cat ./split/%s_%s_%03d_%s.fq.gz ./split/%s_%s_%03d_%s.fq.gz ./split/%s_%s_%03d_%s.fq.gz > ./split_smooth/%s_%s_%03d_%s_SM.fq.gz\n",
                        nm->smpl[k], chnl_code[j], i-1, read_code[k],  /* input  length i-1  */
                        nm->smpl[k], chnl_code[j], i,   read_code[k],  /* input  length i    */
                        nm->smpl[k], chnl_code[j], i+1, read_code[k],  /* input  length i+1  */
                        nm->smpl[k], chnl_code[j], i,   read_code[k]); /* output length i_SM */
            }
        }
    }
    
    //close smooth_transition.sh file
    if ((fclose(out_fp)) == EOF) {
        printf("mk_smooth_script: ERROR - error occurred when closing smoothing shell script. Aborting program...\n");
        abort();
    }
}

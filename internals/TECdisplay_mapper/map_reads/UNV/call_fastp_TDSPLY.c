//
//  call_fastp_TDSPLY.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"

#include "call_fastp_TDSPLY.h"

/* call_fastp_TDSPLY: performs a system call to initiate fastp preprocessing, opens resulting files */
int call_fastp_TDSPLY(TDSPLY_names * nm, fastp_params prms)
{
    static const char smRNA_adpt1[41] = "--adapter_sequence=TGGAATTCTCGGGTGCCAAGG";		//small RNA adapter 1
    static const char smRNA_adpt2[44] = "--adapter_sequence_r2=GATCGTCGGACTGTAGAACTC";	//small RNA adapter 2
    static const char ovrlp_len_req[25] = "--overlap_len_require=15";                   //overlap length requirement
    static const char crrct[13] = "--correction";										//flag for error correction
    static const char merge[8] = "--merge";                                             //flag to merge reads
    static const char umi_STD[35] = "--umi --umi_loc=read2 --umi_len=16";               //standard UMI settings
    static const char umi_BRCD[38] = "--umi --umi_loc=per_read --umi_len=16";           //barcoded UMI settings
    static const char ipt1[3] = "-i";													//option for read 1 input
    static const char ipt2[3] = "-I";													//option for read 2 input
    
    const char * umi = NULL; //pointer to umi command to use
    
    //set umi command
    if (prms.mode == STD_TDSPLY) {
        umi = umi_STD;
    } else if (prms.mode == BRCD_TDSPLY) {
        umi = umi_BRCD;
    } else {
        printf("call_fastp_TDSPLY: error - unrecognized processing mode. aborting...\n");
        abort();
    }
    
    char merged_out[256] = {0};    //array to store merged output option
    char unmerged_out1[256] = {0}; //array to store read 1 output option
    char unmerged_out2[256] = {0}; //array to store read 2 output option
    
    int ret = 0; //variable for storing snprintf return value
    
    //construct output file options
    ret = snprintf(merged_out, 256, "--merged_out ./processed/%s.fq.gz", nm->mrg);
    if (ret >= 256 || ret < 0) {
        printf("call_fastp_TDSPLY: error - error when constructing merged_out name. aborting...\n");
        abort();
    }
    
    ret = snprintf(unmerged_out1, 256, "--out1 ./processed/%s_unmerged.fq.gz", nm->smpl[READ1]);
    if (ret >= 256 || ret < 0) {
        printf("call_fastp_TDSPLY: error - error when constructing unmerged_out1 name. aborting...\n");
        abort();
    }
    
    ret = snprintf(unmerged_out2, 256, "--out2 ./processed/%s_unmerged.fq.gz", nm->smpl[READ2]);
    if (ret >= 256 || ret < 0) {
        printf("call_fastp_TDSPLY: error - error when constructing unmerged_out2 name. aborting...\n");
        abort();
    }
    
    char command[(MAX_LINE*4)]; //array for fastp command. array size exceeds max possible string length
    char lmt[64] = {0};         //read processing limit for debugging purposes
    
    //construct fastp command
    sprintf(command, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s",
            prms.path,			//path for fastp
            smRNA_adpt1,		//adapter 1 trimming sequence
            smRNA_adpt2,		//adapter 2 trimming sequence
            ovrlp_len_req,		//overlap length
            crrct,				//error correction flag
            merge,              //merge read option
            umi, 				// umi settings
            ipt1,				//input 1 flag
            nm->file[READ1],	//filename for read 1
            ipt2,				//input 2 flag
            nm->file[READ2],	//filename for read 2
            merged_out,			//merged output
            unmerged_out1,		//unmerged output read1
            unmerged_out2);		//unmerged output read2
    
    //append optional read processing limit argument
    if (prms.limit) {
        sprintf(lmt, " --reads_to_process=%d", prms.limit);
        strcat(command, lmt);
    }
    
    //record fastp command
    FILE * out_fp = NULL;
    if ((out_fp = fopen("./fastp_command.txt", "w")) == NULL) {
        printf("call_fastp_TDSPLY: error - could not open fastp command file. Aborting program...\n");
        abort();
    }
    
    printf("%s\n", command);			//print fastp command to screen
    fprintf(out_fp, "%s\n", command);	//print fastp command to file
    
    if ((fclose(out_fp)) == EOF ) {
        printf("call_fastp_TDSPLY: error - error occurred when closing fastp command file. Aborting program...\n");
        abort();
    }
    
    system(command); //call fastp using command string
    printf("\n\n");  //print newlines after fastp output message
    return 1;
}

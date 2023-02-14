//
//  call_fastp.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

#include "call_fastp.h"

/* call_fastp: performs a system call to initiate fastp preprocessing, opens resulting files */
int call_fastp(names * nm, FILE ** ifp, fastp_params prms)
{
    static const char umi_MLT[38] = " --umi --umi_loc=per_read --umi_len=9";			//multi-length UMI settings
    static const char umi_SGL[35] = " --umi --umi_loc=read2 --umi_len=9";				//single length UMI settings
    static const char smRNA_adpt1[41] = "--adapter_sequence=TGGAATTCTCGGGTGCCAAGG";		//small RNA adapter 1
    static const char smRNA_adpt2[44] = "--adapter_sequence_r2=GATCGTCGGACTGTAGAACTC";	//small RNA adapter 2
    static const char ovrlp_len_req[25] = "--overlap_len_require=15";
    static const char crrct[13] = "--correction";										//flag for error correction
    static const char ipt1[3] = "-i";													//option for read 1 input
    static const char ipt2[3] = "-I";													//option for read 2 input
    static const char outFiles[46] = "-o ./split/R1out.fq.gz -O ./split/R2out.fq.gz";	//option for output files
    
    char command[(MAX_LINE*4)];        //array for fastp command. array size exceeds max possible string length
    char lmt[64] = {0};        //read processing limit for debugging purposes
    
    //construct fastp command
    sprintf(command, "%s %s %s %s %s %s %s %s %s %s",
            prms.path,			//path for fastp
            smRNA_adpt1,		//adapter 1 trimming sequence
            smRNA_adpt2,		//adapter 2 trimming sequence
            ovrlp_len_req,		//overlap length
            crrct,				//error correction flag
            ipt1,				//input 1 flag
            nm->file[READ1],	//filename for read 1
            ipt2,				//input 2 flag
            nm->file[READ2],	//filename for read 2
            outFiles);			//option for output files
    
    switch (prms.mode) {
        case MULTI:  strcat(command, umi_MLT); break;
        case SINGLE: strcat(command, umi_SGL); break;
        default:
            printf("call_fastp: error - unexpected processing mode. aborting...\n");
            abort();
            break;
    }
    
    //append optional read processing limit argument
    if (prms.limit) {
        sprintf(lmt, " --reads_to_process=%d", prms.limit);
        strcat(command, lmt);
    }
    
    //record fastp command
    FILE * out_fp = NULL;
    if ((out_fp = fopen("./fastp_command.txt", "w")) == NULL) {
        printf("call_fastp: error - could not open fastp command file. Aborting program...\n");
        abort();
    }
    
    printf("%s\n", command);			//print fastp command to screen
    fprintf(out_fp, "%s\n", command);	//print fastp command to file
    
    if ((fclose(out_fp)) == EOF ) {
        printf("call_fastp: error - error occurred when closing fastp command file. Aborting program...\n");
        abort();
    }
    
    system(command);						//call fastp using command string
    system("gunzip ./split/R*out.fq.gz");	//unzip fastp output files
    
    //open fastp output files
    if ((ifp[READ1] = fopen("./split/R1out.fq", "r")) == NULL) {
        printf("error: could not open R1out.fq as read one file. check that fastp path is correct. Aborting program...\n");
        abort();
    }
    
    if ((ifp[READ2] = fopen("./split/R2out.fq", "r")) == NULL) {
        printf("error: could not open R2out.fq as read two file. check that fastp path is correct. Aborting program...\n");
        abort();
    }
    
    return 1;
}

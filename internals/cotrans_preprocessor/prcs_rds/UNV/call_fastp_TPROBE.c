//
//  call_fastp_TPROBE.c
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

#include "call_fastp_TPROBE.h"

/* call_fastp_TPROBE: performs a system call to initiate fastp preprocessing, opens resulting files */
int call_fastp_TPROBE(char * fq1, char * fq2, FILE ** ifp, fastp_params prms, char * fp_out_dir)
{
    char umi_MLT[38]   = " --umi --umi_loc=per_read --umi_len=9";			//multi-length UMI settings
    char umi_SGL[35]   = " --umi --umi_loc=read2 --umi_len=9";				//single length UMI settings
    char umi_MUX[39]   = " --umi --umi_loc=per_read --umi_len=16";         //multiplex UMI settings
    char umi_DNAqc[39] = " --umi --umi_loc=per_read --umi_len=16";         //DNA QC UMI settings
    char smRNA_adpt1[41] = "--adapter_sequence=TGGAATTCTCGGGTGCCAAGG";		//small RNA adapter 1
    char smRNA_adpt2[44] = "--adapter_sequence_r2=GATCGTCGGACTGTAGAACTC";	//small RNA adapter 2
    char ovrlp_len_req[25] = "--overlap_len_require=15";
    char crrct[13] = "--correction";										//flag for error correction
    char ipt1[3] = "-i";													//option for read 1 input
    char ipt2[3] = "-I";													//option for read 2 input
    char outFile1[12] = "R1out.fq.gz";	                                    //read 1 output name
    char outFile2[12] = "R2out.fq.gz";                                     //read 2 output name
    
    char * umi2use = NULL;
    
    char fastp_command[MAX_LINE+1] = {0}; //array for fastp command. array size exceeds max possible string length
    char lmt[64] = {0};                   //read processing limit for debugging purposes
    
    char unzip_command[MAX_LINE+1] = {0};    //array for storing unzip command
    char open_rd1_command[MAX_LINE+1] = {0}; //array for storing open read one file command
    char open_rd2_command[MAX_LINE+1] = {0}; //array for storing open read two file command
    
    int ret = 0; //snprintf return value
    
    //construct fastp command
    ret = snprintf(fastp_command, MAX_LINE, "%s %s %s %s %s %s %s %s %s -o ./%s/%s -O ./%s/%s",
            prms.path,			//path for fastp
            smRNA_adpt1,		//adapter 1 trimming sequence
            smRNA_adpt2,		//adapter 2 trimming sequence
            ovrlp_len_req,		//overlap length
            crrct,				//error correction flag
            ipt1,				//input 1 flag
            fq1,	            //filename for read 1
            ipt2,				//input 2 flag
            fq2,	            //filename for read 2
            fp_out_dir,         //output file directory
            outFile1,           //output file 1 name
            fp_out_dir,         //output file directory
            outFile2);			//output file 2 name

    //test snprintf success
    if (ret >= MAX_LINE || ret < 0) {
        printf("call_fastp_TPROBE: error - failed to construct fastp run command. aborting...\n");
        abort();
    }
    
    //set pointer to UMI processing arguments for current mode
    switch (prms.mode) {
        case MULTI:           umi2use = umi_MLT;   break;
        case SINGLE:          umi2use = umi_SGL;   break;
        case MULTIPLEX:       umi2use = umi_MUX;   break;
        case DNA_PREP_QC_LIB: umi2use = umi_DNAqc; break;
        default:
            printf("call_fastp_TPROBE: error - unexpected processing mode. aborting...\n");
            abort();
            break;
    }
    
    //append mode-specific UMI processing arguments to fastp command
    if ((strlen(fastp_command) + strlen(umi2use)) > MAX_LINE) {
        printf("call_fastp_TPROBE: error - fastp command is too long. aborting...");
        abort();
    } else {
        strcat(fastp_command, umi2use);
    }
    
    //append optional read processing limit argument
    if (prms.limit) {
        sprintf(lmt, " --reads_to_process=%d", prms.limit); //NOTE: string is guaranteed to fit in array
        
        if ((strlen(fastp_command) + strlen(lmt)) > MAX_LINE) {
            printf("call_fastp_TPROBE: error - fastp command is too long. aborting...");
            abort();
        } else {
            strcat(fastp_command, lmt);
        }
    }
    
    //record fastp command
    FILE * out_fp = NULL;
    if ((out_fp = fopen("./fastp_command.txt", "w")) == NULL) {
        printf("call_fastp_TPROBE: error - could not open fastp command file. Aborting program...\n");
        abort();
    }
    
    printf("%s\n", fastp_command);			//print fastp command to screen
    fprintf(out_fp, "%s\n", fastp_command);	//print fastp command to file
    
    if ((fclose(out_fp)) == EOF ) {
        printf("call_fastp_TPROBE: error - error occurred when closing fastp command file. Aborting program...\n");
        abort();
    }
    
    system(fastp_command); //call fastp using command string
    
    //unzip processed fastq files
    ret = snprintf(unzip_command, MAX_LINE, "gunzip ./%s/R*out.fq.gz", fp_out_dir);
    if (ret >= MAX_LINE || ret < 0) {
        printf("call_fastp_TPROBE: error - unable to generate unzipping command. aborting...\n");
        abort();
    } else {
        system(unzip_command);  //unzip fastp output files
    }
    
    int ret_r1 = snprintf(open_rd1_command, MAX_LINE, "./%s/R1out.fq", fp_out_dir);
    int ret_r2 = snprintf(open_rd2_command, MAX_LINE, "./%s/R2out.fq", fp_out_dir);
    
    //test snprintf success
    if (ret_r1 >= MAX_LINE || ret_r1 < 0 || ret_r2 >= MAX_LINE || ret_r2 < 0) {
        printf("call_fastp_TPROBE: error - failed to generate commands for opening fastp output files. aborting...\n");
        abort();
        
    } else {
        
        //open fastp output files
        if ((ifp[READ1] = fopen(open_rd1_command, "r")) == NULL) {
            printf("error: could not open R1out.fq as read one file. Aborting program...\n");
            abort();
        }
        
        if ((ifp[READ2] = fopen(open_rd2_command, "r")) == NULL) {
            printf("error: could not open R2out.fq as read two file. Aborting program...\n");
            abort();
        }
    }
    return 1;
}

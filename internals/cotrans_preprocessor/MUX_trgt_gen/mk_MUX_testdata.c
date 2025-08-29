//
//  mk_MUX_testdata.c
//  
//
//  Created by Eric Strobel on 6/27/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"
#include "../../variant_maker/constant_seqs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../MLT_trgt_gen/mk_3pEnd_testdata.h"

#include "./mk_MUX_trgts.h"

#include "mk_MUX_testdata.h"

/* mk_MUX_testdata: generate TECprobe-MUX test data */
void mk_MUX_testdata(TPROBE_names * nm, compact_target * ctrg, int ctrg_cnt)
{
    srand(time(NULL)); //seed pseudorandom number generator
    
    char out_dir[MAX_LINE+1] = {0}; //array to store output directory
    
    int ret = 0;             //variable for storing snprintf return
    int read_pair_count = 0; //number of read pairs generated
    
    //construct output directory name
    ret = snprintf(out_dir, MAX_LINE, "%s_test_data", nm->trgts_prfx);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_MUX_testdata: error - error when constructing output directory name. aborting...\n");
        abort();
    }
    mk_out_dir(out_dir); //make output directory
    
    FILE * out_rd1 = NULL; //read 1 file pointer
    FILE * out_rd2 = NULL; //read 2 file pointer

    //generate read 1 output file
    ret = snprintf(nm->file[READ1], MAX_LINE+1, "./%s/%s_testdata_R1.fq", out_dir, nm->trgts_prfx);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_MUX_testdata: error - error when constructing output file name. aborting...\n");
        abort();
    }
    if ((out_rd1 = fopen(nm->file[READ1], "w")) == NULL) {
        printf("mk_MUX_testdata: error - could not generate test data read one file. Aborting program...\n");
        abort();
    }
    
    //generate read 2 output file
    ret = snprintf(nm->file[READ2], MAX_LINE+1, "./%s/%s_testdata_R2.fq", out_dir, nm->trgts_prfx);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_MUX_testdata: error - error when constructing output file name. aborting...\n");
        abort();
    }
    if ((out_rd2 = fopen(nm->file[READ2], "w")) == NULL) {
        printf("mk_MUX_testdata: error - could not generate test data read two file. Aborting program...\n");
        abort();
    }
    
    int i = 0;             //general purpose index
    int j = 0;             //general purpose index
    char chnl_bc[6] = {0}; //array for storing channel barcode, append to head of read 1

    //arrays that direct the generation of five types of channel barcodes for test data reads
    static const int mk_chnl[CHNL_CLASSES] = {  UNT,  MOD,  UNT,  MOD,  ERR};
    static const int mk_mtch[CHNL_CLASSES] = { FULL, FULL, PART, PART, FULL};
    
    //generate test data reads
    for (i = 0; i < ctrg_cnt; i++) { //for every compact target
        if (!ctrg[i].mul) {          //if the target was not identified as redundant
            
            //generate one test data read for each of the five types of channel barcodes
            for (j = 0; j < CHNL_CLASSES; j++) {
                mk_rndmzd_bc(chnl_bc, mk_chnl[j], mk_mtch[j]);
                read_pair_count = print_MUX_fq(out_rd1, out_rd2, chnl_bc, &ctrg[i]);
            }
        }
    }
    
    //close read 1 output file
    if ((fclose(out_rd1)) == EOF) {
        printf("mk_MUX_testdata: error - error occurred when closing test data read 1 file. Aborting program...\n");
        abort();
    }
    
    //close read 2 output file
    if ((fclose(out_rd2)) == EOF) {
        printf("mk_MUX_testdata: error - error occurred when closing test data read 2 file. Aborting program...\n");
        abort();
    }
    
    printf("made %d test data reads\n", read_pair_count);
    
}

/* print_MUX_fq: print TECprobe-MUX test data read to output fastq file */
int print_MUX_fq(FILE * out_rd1, FILE * out_rd2, char * chnl_bc, compact_target * ctrg)
{
    extern char sc1[17];             //structure cassette 1 sequence
    extern char RLA29synch_3p11[34]; //TECprobe-MX linker, v2
    
    static int cnt = 0; //number of read pairs generated
    
    char brcd[MAX_LINE+1] = {0}; //array to store TECprobe-MX barcode
    
    int i = 0;   //general purpose index
    int ret = 0; //variable for storing snprintf return
    
    char rd1[MAX_LINE+1] = {0};    //array for generating read 1 sequence
    char rd2[MAX_LINE+1] = {0};    //array for generating read 2 sequence
        
    opt_BC * p_opt_BC_crnt = (opt_BC *)ctrg->opt;               //set current target optional barcode values pointer
    opt_BC * p_opt_BC_ref  = (opt_BC *)p_opt_BC_crnt->ref->opt; //set reference target optional barcode values pointer
    
    bin2seq(brcd, &ctrg->bsq, MAX_LINE+1); //convert binary-encoded barcode sequence to character string
    
    //generate test data read 2, which contains:
    //1. channel barcode
    //2. structure cassette 1
    //3. target sequence of the reference target from which the current barcode was derived
    //4. RLA29synch_3p11 linker
    //5. current barcode
    ret = snprintf(rd2, MAX_LINE+1, "%s%s%s%s%s", chnl_bc, sc1, p_opt_BC_ref->tsq, RLA29synch_3p11, brcd);
    if (ret >= MAX_LINE || ret < 0) {
        printf("print_MUX_fq: error - error when constructing test data read. aborting...\n");
        abort();
    }
    
    for (i = 0; rd2[i]; i++) {    //convert test data read to uppercase
        rd2[i] = toupper(rd2[i]);
    }
    reverse_complement(rd1, rd2, REVCOMP); //revcomp read 2 to obtain read 1 sequence
        
    //print fastq line 1 (read id)
    //read id contains:
    //1. prefix 'testdata_barcode' to indicate that the subsequent number is a barcode id
    //2. id of the barcode target used to generate the read
    //3. generalized Illumina read suffix
    fprintf(out_rd1, "@testdata_barcode=%020llu_R1 1:N:0:INDEX\n", (long long unsigned int)ctrg->bid);
    fprintf(out_rd2, "@testdata_barcode=%020llu_R2 2:N:0:INDEX\n", (long long unsigned int)ctrg->bid);
    
    //print fastq line 2 (read sequence)
    fprintf(out_rd1, "%s\n", rd1);
    fprintf(out_rd2, "%s\n", rd2);
    
    //print fastq line 3
    fprintf(out_rd1, "+\n");
    fprintf(out_rd2, "+\n");
    
    //print fastq line 4 (qscore, set as string of 'I' that matches sequence length)
    int len = strlen(rd1);
    for (i = 0; i < (len); i++) {
        fprintf(out_rd1, "I");
        fprintf(out_rd2, "I");
    }
    fprintf(out_rd1, "\n");
    fprintf(out_rd2, "\n");
    
    return ++cnt; //increment read output counter
}

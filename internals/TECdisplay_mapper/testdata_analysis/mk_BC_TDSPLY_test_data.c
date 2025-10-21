//
//  mk_BC_TDSPLY_test_data.c
//  
//
//  Created by Eric Strobel on 10/16/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"
#include "../../seq_utils/revcomp.h"

#include "../../variant_maker/variant_maker_defs.h"

#include "../TECdisplay_mapper_structs.h"
#include "../TECdisplay_mapper_defs.h"

#include "../../utils/debug.h"

#include "./mk_TDSPLY_test_data.h"

#include "mk_BC_TDSPLY_test_data.h"

extern int debug;

/* mk_BC_TDSPLY_test_data: coordinates test data generation */
int mk_BC_TDSPLY_test_data(TDSPLY_names * nm, compact_target *ctrg, target_params *trg_prms, int ctrg_cnt)
{
    srand(time(NULL)); //seed pseudorandom number generator
    
    char out_dir[MAX_LINE+1] = {0}; //array to store output directory name
    
    int ret = 0; //variable for storing snprintf return value

    //construct output directory name
    ret = snprintf(out_dir, MAX_LINE, "%s_test_data", nm->trgs_prefix);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_BC_TDSPLY_test_data: error - error when constructing output directory name. aborting...\n");
        abort();
    }
    mk_out_dir(out_dir); //make output directory
    
    FILE * out_rd1 = NULL; //read 1 file pointer
    FILE * out_rd2 = NULL; //read 2 file pointer

    //generate read 1 output file
    ret = snprintf(nm->file[READ1], MAX_LINE, "./%s/%s_testdata_R1.fq", out_dir, nm->trgs_prefix);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_BC_TDSPLY_test_data: error - error when constructing read one output file name. aborting...\n");
        abort();
    }
    if ((out_rd1 = fopen(nm->file[READ1], "w")) == NULL) {
        printf("mk_BC_TDSPLY_test_data: error - could not generate test data read one file. Aborting program...\n");
        abort();
    }
    
    //generate read 2 output file
    ret = snprintf(nm->file[READ2], MAX_LINE, "./%s/%s_testdata_R2.fq", out_dir, nm->trgs_prefix);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_BC_TDSPLY_test_data: error - error when constructing read two output file name. aborting...\n");
        abort();
    }
    if ((out_rd2 = fopen(nm->file[READ2], "w")) == NULL) {
        printf("mk_BC_TDSPLY_test_data: error - could not generate test data read two file. Aborting program...\n");
        abort();
    }
    
    //channel barcode generation variables
    static const int mk_chnl[CHNL_CLASSES] = {  UNB,  BND,  UNB,  BND,  ERR};  //channel type
    static const int mk_mtch[CHNL_CLASSES] = { FULL, FULL, PART, PART, FULL};  //channel match type
        
    char chnl_bc[5] = {0}; //array for storing channel barcode, append to head of read 1

    int i = 0; //general purpose index
    int j = 0; //general purpose index
        
    //generate test data reads. test data is only made using non-redundant targets.
    //five mappable reads and are made for each non-redundant target. the five
    //mappable reads contain the following randomized channel barcodes:
    
    //  1. bound, full match
    //  2. bound, partial match
    //  3. unbound, full match
    //  4. unbound, partial match
    //  5. unmappable
     
    //generate test data reads
    for (i = 0; i < ctrg_cnt; i++) { //for every compact target
        if (!ctrg[i].mul) {          //if the target was not identified as redundant
            
            //generate one test data read for each of the five types of channel barcodes
            for (j = 0; j < CHNL_CLASSES; j++) {
                mk_rndmzd_TDSPLY_bc(chnl_bc, mk_chnl[j], mk_mtch[j]);
                
                //make reads from native sequence
                print_BC_TDSPLY_fq(out_rd1, out_rd2, &ctrg[i], chnl_bc, NAT);
            }
        }
    }
    
    //close read 1 output file
    if ((fclose(out_rd1)) == EOF) {
        printf("mk_BC_TDSPLY_test_data: error - error occurred when closing test data read 1 file. Aborting program...\n");
        abort();
    }
    
    //close read 2 output file
    if ((fclose(out_rd2)) == EOF) {
        printf("mk_BC_TDSPLY_test_data: error - error occurred when closing test data read 2 file. Aborting program...\n");
        abort();
    }
    
    return 1;
}


/* print_BC_TDSPLY_fq: construct read sequences and print to fastq file */
void print_BC_TDSPLY_fq(FILE * out_rd1, FILE * out_rd2, compact_target * ctrg, char * chnl_bc, int end_rnd_typ) {
    
    static int bc_cnt = 0; //number of read pairs generated
    
    int i = 0; //general purpose index
    
    char rd1[MAX_LINE+1] = {0}; //array for generating read 1 sequence
    char rd2[MAX_LINE+1] = {0}; //array for generating read 2 sequence
    
    const char UMI[13] = "CATCATCATCAT"; //12 nt UMI placeholder, append after channel barcode
    const char prm[51] = "TTATCAAAAAGAGTATTGACTCTTTTACCTCTGGCGGTGATAATGGTTGC"; //PRA1 promoter
    const char ldr[34] = "ATTCGGTGCTCTTCTCTTCGGCCTTCGGGCCAA"; //C3SC1 leader sequence, append 5' end of insert
    
    opt_BC * p_opt_BC_crnt = (opt_BC *)ctrg->opt;               //set current target optional barcode values pointer
    opt_BC * p_opt_BC_ref  = (opt_BC *)p_opt_BC_crnt->ref->opt; //set reference target optional barcode values pointer
    
    char brcd[MAX_LINE+1] = {0};           //array to store barcode
    bin2seq(brcd, &ctrg->bsq, MAX_LINE+1); //convert binary-encoded barcode sequence to character string
    
    int ret = 0; //variable for storing snprintf return
    
    //construct read 2
    ret = snprintf(rd2, MAX_LINE, "%s%s%s%s%s%s", chnl_bc, UMI, prm, ldr, p_opt_BC_ref->tsq, brcd); //construct read 2
    if (ret >= MAX_LINE || ret < 0) {
        printf("print_BC_TDSPLY_fq: error - error when constructing test data read. aborting...\n");
        abort();
    }
    
    reverse_complement(rd1, rd2, REVCOMP); //revcomp read 2 to obtain read 1 sequence
    
    //print fastq line 1 (read id)
    //read id contains:
    //1. prefix 'testdata_barcode' to indicate read is for use as barcoded TECdisplay test data
    //2. id of the barcode variant used to generate the read
    //5. generalized Illumina read suffix
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
    for (i = 0; i < len; i++) {
        fprintf(out_rd1, "I");
        fprintf(out_rd2, "I");
    }
    fprintf(out_rd1, "\n");
    fprintf(out_rd2, "\n");
    
    bc_cnt++; //increment read output counter
}

//
//  mk_DNA_qc_testdata.c
//  
//
//  Created by Eric Strobel on 9/10/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"
#include "../../seq_utils/revcomp.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

#include "mk_DNA_qc_testdata.h"

/* mk_DNA_qc_testdata: generate TECprobe-MUX test data */
void mk_DNA_qc_testdata(TPROBE_names * nm, compact_target * ctrg, int ctrg_cnt)
{
    srand(time(NULL)); //seed pseudorandom number generator
    
    char out_dir[MAX_LINE+1] = {0}; //array to store output directory
    
    int ret = 0;             //variable for storing snprintf return
    int read_pair_count = 0; //number of read pairs generated
    
    //construct output directory name
    ret = snprintf(out_dir, MAX_LINE, "%s_test_data", nm->trgts_prfx);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_DNA_qc_testdata: error - error when constructing output directory name. aborting...\n");
        abort();
    }
    mk_out_dir(out_dir); //make output directory
    
    FILE * out_rd1 = NULL; //read 1 file pointer
    FILE * out_rd2 = NULL; //read 2 file pointer

    //generate read 1 output file
    ret = snprintf(nm->file[READ1], MAX_LINE+1, "./%s/%s_testdata_R1.fq", out_dir, nm->trgts_prfx);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_DNA_qc_testdata: error - error when constructing output file name. aborting...\n");
        abort();
    }
    if ((out_rd1 = fopen(nm->file[READ1], "w")) == NULL) {
        printf("mk_DNA_qc_testdata: error - could not generate test data read one file. Aborting program...\n");
        abort();
    }
    
    //generate read 2 output file
    ret = snprintf(nm->file[READ2], MAX_LINE+1, "./%s/%s_testdata_R2.fq", out_dir, nm->trgts_prfx);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_DNA_qc_testdata: error - error when constructing output file name. aborting...\n");
        abort();
    }
    if ((out_rd2 = fopen(nm->file[READ2], "w")) == NULL) {
        printf("mk_DNA_qc_testdata: error - could not generate test data read two file. Aborting program...\n");
        abort();
    }
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    //generate test data reads
    for (i = 0; i < ctrg_cnt; i++) { //for every compact target
        if (!ctrg[i].mul) {          //if the target was not identified as redundant
            read_pair_count = print_DNA_qc_fq(out_rd1, out_rd2, &ctrg[i], &ctrg[0], ctrg_cnt, 1);
            read_pair_count = print_DNA_qc_fq(out_rd1, out_rd2, &ctrg[i], &ctrg[0], ctrg_cnt, 0);
        }
    }
    
    //close read 1 output file
    if ((fclose(out_rd1)) == EOF) {
        printf("mk_DNA_qc_testdata: error - error occurred when closing test data read 1 file. Aborting program...\n");
        abort();
    }
    
    //close read 2 output file
    if ((fclose(out_rd2)) == EOF) {
        printf("mk_DNA_qc_testdata: error - error occurred when closing test data read 2 file. Aborting program...\n");
        abort();
    }
    
    printf("made %d test data reads\n", read_pair_count);

    return;
}

/* print_DNA_qc_fq: print TECprobe-MUX test data read to output fastq file */
int print_DNA_qc_fq(FILE * out_rd1, FILE * out_rd2, compact_target * ctrg, compact_target * root_ctrg, int ctrg_cnt, int concordant)
{
    static int cnt = 0; //number of read pairs generated
    
    char brcd[MAX_LINE+1] = {0}; //array to store barcode
    
    int i = 0;   //general purpose index
    int ret = 0; //variable for storing snprintf return
    
    char rd1[MAX_LINE+1] = {0};    //array for generating read 1 sequence
    char rd2[MAX_LINE+1] = {0};    //array for generating read 2 sequence
        
    opt_BC * p_opt_BC_crnt = (opt_BC *)ctrg->opt;               //set current target optional barcode values pointer
    opt_BC * p_opt_BC_ref  = (opt_BC *)p_opt_BC_crnt->ref->opt; //set reference target optional barcode values pointer
    
    //use ref target association struct to track number of concordant/discordant reads generated
    compact_target * ref_ctrg = p_opt_BC_crnt->ref;       //set pointer to reference ctrg
    association * assoc = (association *)(ref_ctrg->utl); //set pointer to association struct
    if (concordant) {                                     //record whether current testa read is concordant or discordant
        assoc->src_cncrdnt++;
    } else {
        assoc->src_dscrdnt++;
    }
        
    compact_target * lnkd_ctrg  = (compact_target *)p_opt_BC_crnt->lnk; //set pointer to linked target
    compact_target * dscrdnt_ctrg = NULL;                               //pointer for setting discordant target
    
    bin2seq(brcd, &ctrg->bsq, MAX_LINE); //convert binary-encoded barcode sequence to character string
    
    char * p_bc1 = NULL;       //pointer to barcode 1
    char * p_bc2 = NULL;       //pointer to barcode 2
    uint64_t * p_bid1 = NULL;  //pointer to barcode 1 id
    uint64_t * p_bid2 = NULL;  //pointer to barcode 2 id
    
    //set barcode 1 and 2 pointers based on:
    //  1. whether the input barcode is barcode 1 or 2
    //  2. whether the the read will contain concordant or discordant barcodes
    
    if (p_opt_BC_crnt->num == 1) {     //input barcode is barcode 1
        p_bc1 = ctrg->csq;             //set input barcode sequence as barcode 1
        p_bid1 = &ctrg->bid;           //set input barcode id as barcode 1 id
        
        if (concordant) {              //read barcodes are to be concordant
            p_bc2 = lnkd_ctrg->csq;    //set linked barcode as barcode 2
            p_bid2 = &lnkd_ctrg->bid;  //set linked barcode id as barcode 2 id
            
        } else {
            //setting discordant barcode
            dscrdnt_ctrg = find_discordant_ctrg(root_ctrg, ctrg, ctrg_cnt); //select random discordant barcode
            p_bc2 = dscrdnt_ctrg->csq;                                      //set discordant barcode as barcode 2
            p_bid2 = &dscrdnt_ctrg->bid;                                    //set discordant barcode id as barcode 2 id
        }
        
    } else if (p_opt_BC_crnt->num == 2) { //input barcode is barcode 2
        p_bc2 = ctrg->csq;                //set input barcode sequence as barcode 2
        p_bid2 = &ctrg->bid;              //set input barcode id as barcode 2 id
        
        if (concordant) {                 //reads are to be concordant
            p_bc1 = lnkd_ctrg->csq;       //set linked barcode as barcode 1
            p_bid1 = &lnkd_ctrg->bid;     //set linked barcode id as barcode 1 id
            
        } else {
            //setting discordant barcode
            dscrdnt_ctrg = find_discordant_ctrg(root_ctrg, ctrg, ctrg_cnt); //select random discordant barcode
            p_bc1 = dscrdnt_ctrg->csq;                                      //set discordant barcode as barcode 1
            p_bid1 = &dscrdnt_ctrg->bid;                                    //set discordant barcode id as barcode 1 id
        }
        
    } else { //throw error if input barcode number is not 1 or 2
        printf("print_DNA_qc_fq: error - expected barcode number to be 1 or 2. aborting...\n");
        abort();
    }
    
    //generate test data read 2, which contains:
    //1. barcode 1
    //2. insert sequence
    //3. barcode 2
    ret = snprintf(rd2, MAX_LINE, "%s%s%s", p_bc1, p_opt_BC_ref->tsq, p_bc2);
    if (ret >= MAX_LINE || ret < 0) {
        printf("print_DNA_qc_fq: error - error when constructing test data read. aborting...\n");
        abort();
    }
    
    for (i = 0; rd2[i]; i++) {    //convert test data read to uppercase
        rd2[i] = toupper(rd2[i]);
    }
    reverse_complement(rd1, rd2, REVCOMP); //revcomp read 2 to obtain read 1 sequence
        
    //print fastq line 1 (read id)
    //read id contains:
    //1. prefix 'testdata_barcode' to indicate that the subsequent number are barcode ids
    //2. ids of the barcode targets used to generate the read
    //3. generalized Illumina read suffix
    fprintf(out_rd1, "@testdata_barcode=%020llu_barcode2=%020llu_R1 1:N:0:INDEX\n", (long long unsigned int)p_bid1, (long long unsigned int)p_bid2);
    fprintf(out_rd2, "@testdata_barcode=%020llu_barcode2=%020llu_R2 2:N:0:INDEX\n", (long long unsigned int)p_bid1, (long long unsigned int)p_bid2);
    
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
    
    return ++cnt; //increment read output counter
}

/* find_discordant_ctrg: find barcode ctrg that is discordant relative to the input ctrg*/
compact_target * find_discordant_ctrg(compact_target * root_ctrg, compact_target * ipt_ctrg, int ctrg_cnt)
{
    opt_BC * ipt_BC_val = (opt_BC *)(ipt_ctrg->opt); //pointer to input ctrg optional values
        
    compact_target * dscrdnt_ctrg = NULL;     //pointer to discordant compact_target
    compact_target * dscrdnt_ref_ctrg = NULL; //pointer to discordant compact_target reference
    opt_BC * dscrdnt_BC_val = NULL;           //pointer to discordant compact_target optional values
    association * dscrdnt_assoc = NULL;       //pointer to discordant compact_target association structure
    
    while (1) { //TODO: set a maximum number of attempts and use it to control the loop?
        
        dscrdnt_ctrg = &root_ctrg[rand() % ctrg_cnt];   //point dscrdnt_ctrg to random ctrg
        dscrdnt_BC_val = (opt_BC *)(dscrdnt_ctrg->opt); //set pointer to optional values struct
        
        //test whether the candidate discordant ctrg barcode number is complementary to
        //the input ctrg barcode number and whether the discordant ctrg was marked as
        //blacklisted or redundant. if these tests are passed, test whether the barcode
        //is discordant relative to the input barcode
        if (((ipt_BC_val->num == 1 && dscrdnt_BC_val->num == 2)  ||
             (ipt_BC_val->num == 2 && dscrdnt_BC_val->num == 1)) &&
              !dscrdnt_ctrg->bl && !dscrdnt_ctrg->mul) {
            
            //test for discordance by comparing the ctrg linked of the discordant ctrg
            //candidate to the native ctrg of the input ctrg. if these are distinct,
            //the barcodes are discordant
            if ((uint64_t)(dscrdnt_BC_val->lnk) != (uint64_t)(ipt_BC_val->ref) ) {
                dscrdnt_ref_ctrg = dscrdnt_BC_val->ref;                 //set pointer to discordant ctrg reference ctrg
                dscrdnt_assoc = (association *)(dscrdnt_ref_ctrg->utl); //set pointer to association structure of ref
                dscrdnt_assoc->rnd_dscrdnt++;                           //recorde usage of ctrg as discordant ctrg
                return dscrdnt_ctrg;                                    //return pointer to the discordant ctrg
            }
        }
    }
}

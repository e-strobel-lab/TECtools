//
//  mk_MUX_trgts.c
//  
//
//  Created by Eric Strobel on 6/25/25.
//

#include <stdio.h>
#include <string.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/io_management.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../../variant_maker/make_barcodes.h"
#include "../../variant_maker/constant_seqs.h"
#include "../../variant_maker/vmt_suffix.h"

#include "../../TECdisplay_mapper/map_reads/map_expected/parse_vmt_trgts.h"

#include "./parse_fa_trgts.h"
#include "./mk_barcoded_target_fastas.h"

#include "mk_MUX_trgts.h"

//global variables
uint64_t trgt_ID_mask = 0xFFFFFFFFFFFF0000; //barcode id mask that isolates the reference barcode id
uint64_t mutcode_mask = 0x000000000000FFFF; //barcode id mask that isolates mutation code

/* mk_MUX_trgts: manages TECprobe-MUX target generation */
int mk_MUX_trgts(TPROBE_names * nm, target * refs, opt_ref * ref_val, compact_target * ctrg, opt_BC * BC_val, FILE * fp_MUXtrgs, int trgt_ftype, target_params * trg_prms, int clcd_ctrg_cnt, TDSPLY_fasta * wt)
{
    extern char RLA29synch_3p11[34]; //linker sequence for TECprobe-MUX experiments
    
    int brcd_len = 0;  //barcode length
    int ctrg_cnt = 0;  //number of barcode targets
    
    
    //parse barcoded targets file and store as compact targets
    if (trgt_ftype == VMT_FILE) {
        parse_vmt_trgts(fp_MUXtrgs, trgt_ftype, refs, ref_val, ctrg, BC_val, trg_prms, wt, TPROBE_TRGS);
        
    } else if (trgt_ftype == FASTA_FILE) {
        parse_fa_trgts(fp_MUXtrgs, trgt_ftype, ctrg, BC_val, trg_prms);
        
    } else {
        printf("mk_MUX_trgts: error - unrecognized barcoded targets file type. aborting...\n");
    }

    mk_barcoded_target_fastas(nm, ctrg, trg_prms); //generate individual barcoded target fasta files
    ctrg_cnt = trg_prms->t_cnt;                    //set ctrg_cnt, which tracks total ctrgs below
    
    printf("barcode count = %d\ncalc'd trg count = %d\n", trg_prms->t_cnt, clcd_ctrg_cnt);
    
    uint64_t i = 0;        //general purpose index
    uint64_t mutCode = 0;  //variable for storing mutation code. each barcode id is an unsigned 64 bit integer.
                           //bits 0-15 are used to store a mutation code. bits >= 16 are used to store the id
                           //of the input barcode from which the mutated barcode is derived
    
    int tot_SUB_trgts = 0; //number of substitution targets generated
    int tot_INS_trgts = 0; //number of insertion targets generated
    int tot_DEL_trgts = 0; //number of deletion targets generated
    
    for (i = 0; i < trg_prms->t_cnt; i++) { //for every input barcode
        
        mutCode = 1; //set mutCode to 1 (input barcodes have mutCode 0)
        
        //generate sub and indel targets for the current barcode
        tot_SUB_trgts += mk_SUB_trgts(ctrg, BC_val, &ctrg_cnt, &ctrg[i], trg_prms->BClen, &mutCode);
        tot_INS_trgts += mk_INS_trgts(ctrg, BC_val, &ctrg_cnt, &ctrg[i], trg_prms->BClen, &mutCode);
        tot_DEL_trgts += mk_DEL_trgts(ctrg, BC_val, &ctrg_cnt, &ctrg[i], trg_prms->BClen, &mutCode, toupper(RLA29synch_3p11[strlen(RLA29synch_3p11)-1]));
    }
    
    //print the number of targets that were generated
    printf("total targets = %d\n", trg_prms->t_cnt + tot_SUB_trgts + tot_INS_trgts + tot_DEL_trgts);
    printf("%d compact targets\n", ctrg_cnt);
    
    return ctrg_cnt; //return the number of compact targets that were generated
}

/* mk_SUB_trgts: generates single substitution targets for input barcode */
int mk_SUB_trgts(compact_target * ctrg, opt_BC * BC_val, int * ctrg_cnt, compact_target * src_ctrg, int brcd_len, uint64_t * mutCode)
{
    int sub_trgts_made = 0; //number of substitution targets that were made
    
    int pos2var = 0; //position to vary
    int i = 0;       //general purpose index
    
    char var[3][MAX_LINE+1] = {{0}};    //array for constructing substitution targets
    
    if (*mutCode != MIN_SUB_CODE) {
        printf("mk_SUB_trgts: error - incorrect minimum substitution code (%llu). expected %d. aborting...\n", (long long unsigned int)(*mutCode), MIN_SUB_CODE);
        abort();
    }
    
    for (pos2var = 0; pos2var < brcd_len && src_ctrg->csq[pos2var]; pos2var++) { //for each position of the input barcode
        for (i = 0; i < brcd_len; i++) {     //iterate through barcode sequence
            if (i == pos2var) {              //if index i is being varied in this interation
                switch (src_ctrg->csq[i]) {  //substitute src_ctrg->csq[i] with the other 3 bases in var array
                    case 'A':
                        var[0][i] = 'T';
                        var[1][i] = 'G';
                        var[2][i] = 'C';
                        break;
                    case 'T':
                        var[0][i] = 'A';
                        var[1][i] = 'G';
                        var[2][i] = 'C';
                        break;
                    case 'G':
                        var[0][i] = 'A';
                        var[1][i] = 'T';
                        var[2][i] = 'C';
                        break;
                    case 'C':
                        var[0][i] = 'A';
                        var[1][i] = 'T';
                        var[2][i] = 'G';
                        break;
                    default:
                        printf("printSubs: error - unexpected character %c in barcode sequence. aborting...\n", src_ctrg->csq[i]);
                        abort();
                        break;
                };
            } else { //not varying index i, copy native base to var array
                var[0][i] = src_ctrg->csq[i];
                var[1][i] = src_ctrg->csq[i];
                var[2][i] = src_ctrg->csq[i];
            }
        }
        
        //add terminating null characters to each var array
        var[0][i] = '\0';
        var[1][i] = '\0';
        var[2][i] = '\0';
        
        //check that barcode was read completely
        if (i != brcd_len) {
            printf("mk_SUB_trgts: error - unexpected barcode length. aborting...\n");
            abort();
        }
        
        //store substitution variants as binary encoded sequences in compact targets
        for (i = 0; i < 3; i++) {
            set_BC_val(&ctrg[*ctrg_cnt], &BC_val[*ctrg_cnt], NULL, src_ctrg, MUTANT_BRCD); //set ptr to BC_val structure
            ctrg[*ctrg_cnt].bid = src_ctrg->bid | (*mutCode)++; //generate and store barcode id
            seq2bin_long(var[i], &ctrg[*ctrg_cnt].bsq, 1);      //store binary-encoded barcode sequence
            BC_val[*ctrg_cnt].typ = SUB;                        //set barcode type to substitution
            
            //allocate memory for storing character-encoded barcode
            if ((ctrg[*ctrg_cnt].csq = malloc((strlen(var[i])+1) * sizeof(*ctrg[*ctrg_cnt].csq))) == NULL) {
                printf("mk_SUB_trgts: error - memory allocation for character-encoded target sequence failed. aborting...\n");
                abort();
            }
            strcpy(ctrg[*ctrg_cnt].csq, var[i]); //store character-encoded barcode
            
            (*ctrg_cnt)++; //increment compact target count
        }
        sub_trgts_made += 3; //increment number of substitution targets by three
    }
    
    if (((*mutCode)-1) != MAX_SUB_CODE) {
        printf("mk_SUB_trgts: error - incorrect maximum substitution code (%llu). expected %d. aborting...\n", (long long unsigned int)((*mutCode)-1), MAX_SUB_CODE);
        abort();
    }
    
    return sub_trgts_made; //return number of substitution targets that were made
}

/* mk_INS_trgts: generate single insertion targets for input barcode */
int mk_INS_trgts(compact_target * ctrg, opt_BC * BC_val, int * ctrg_cnt, compact_target * src_ctrg, int brcd_len, uint64_t * mutCode)
{
    int ins_trgts_made = 0;  //number of insertion targets that were made
    
    int pos2var = 0; //current position in target to add insertion
    int i = 0;       //general purpose index
    int j = 0;       //general purpose index
    
    char var[4][MAX_LINE+1] = {{0}}; //array for constructing insertion targets
    
    if (*mutCode != MIN_INS_CODE) {
        printf("mk_INS_trgts: error - incorrect minimum insertion code (%llu). expected %d. aborting...\n", (long long unsigned int)(*mutCode), MIN_INS_CODE);
        abort();
    }
    
    for (pos2var = 0; pos2var < brcd_len && src_ctrg->csq[pos2var]; pos2var++) { //for each position of the input barcode
        
        //iterate through barcode sequence
        //j starts at 1 because first input barcode character is always removed
        for (i = 0, j = 1; i < brcd_len && j <= brcd_len; i++) {
            if (i == pos2var) {  //if index i is being varied in this iteration, insert nucleotide at this position
                var[0][i] = 'A';
                var[1][i] = 'T';
                var[2][i] = 'G';
                var[3][i] = 'C';
            } else {             //otherwise, copy input sequence
                var[0][i] = src_ctrg->csq[j];
                var[1][i] = src_ctrg->csq[j];
                var[2][i] = src_ctrg->csq[j];
                var[3][i] = src_ctrg->csq[j];
                j++;
            }
        }
        
        //add terminating null character
        var[0][i] = '\0';
        var[1][i] = '\0';
        var[2][i] = '\0';
        var[3][i] = '\0';
                
        //check that input and variant barcode indices have reached full length
        if (i != brcd_len || j != brcd_len) {
            printf("mk_INS_trgts: error - unexpected barcode length. aborting...\n");
            abort();
        }
        
        //store insertion variants as binary encoded sequences in compact targets
        for (i = 0; i < 4; i++) {
            set_BC_val(&ctrg[*ctrg_cnt], &BC_val[*ctrg_cnt], NULL, src_ctrg, MUTANT_BRCD); //set ptr to BC_val structure
            ctrg[*ctrg_cnt].bid = src_ctrg->bid | (*mutCode)++;  //generate and store barcode id
            seq2bin_long(var[i], &ctrg[*ctrg_cnt].bsq, 1);       //store binary-encoded barcode sequence
            BC_val[*ctrg_cnt].typ = INS;                         //set barcode type to insertion
            
            //allocate memory for storing character-encoded barcode
            if ((ctrg[*ctrg_cnt].csq = malloc((strlen(var[i])+1) * sizeof(*ctrg[*ctrg_cnt].csq))) == NULL) {
                printf("mk_INS_trgts: error - memory allocation for character-encoded target sequence failed. aborting...\n");
                abort();
            }
            strcpy(ctrg[*ctrg_cnt].csq, var[i]); //store character-encoded barcode
            
            (*ctrg_cnt)++; //increment compact target count
        }
        ins_trgts_made += 4; //incrmement number of insertion targets by four
    }
    
    if (((*mutCode)-1) != MAX_INS_CODE) {
        printf("mk_INS_trgts: error - incorrect maximum insertion code (%llu). expected %d. aborting...\n", (long long unsigned int)((*mutCode)-1), MAX_INS_CODE);
        abort();
    }
    
    return ins_trgts_made; //return number of insertion targets made
}

/* mk_DEL_trgts: generate single deletion targets for input barcode */
int mk_DEL_trgts(compact_target * ctrg, opt_BC * BC_val, int * ctrg_cnt, compact_target * src_ctrg, int brcd_len, uint64_t * mutCode, char upstrm_nt)
{
    int del_trgts_made = 0;  //number of deletion targets made
    
    int pos2var = 0; //current position in target to delete
    int i = 0;       //general purpose index
    int j = 0;       //general purpose index
    
    char var[MAX_LINE+1]; //array for constructing deletion targets
    
    if (*mutCode != MIN_DEL_CODE) {
        printf("mk_DEL_trgts: error - incorrect minimum deletion code (%llu). expected %d. aborting...\n", (long long unsigned int)(*mutCode), MIN_DEL_CODE);
        abort();
    }
    
    for (pos2var = 0; pos2var < brcd_len && src_ctrg->csq[pos2var]; pos2var++) { //for each position of the input barcode
        
        i = 0;                //set i to zero
        var[i++] = upstrm_nt; //set first variant base to the nt immediately upstream, since one nt will be deleted
        
        for (j = 0; i <= brcd_len && j < brcd_len; j++) { //iterate through barcode sequence
            if (j == pos2var) {                           //if pos2var is reached, do nothing so the nt is deleted
                ;
            } else {                       //othewise,
                var[i] = src_ctrg->csq[j]; //copy input barcode nt
                i++;
            }
        }
        var[i] = '\0'; //add terminating null
        
        //check that input and variant barcode indices have reached full length
        if (i != brcd_len || j != brcd_len) {
            printf("mk_DEL_trgts: error - unexpected barcode length. aborting...\n");
            abort();
        }
        
        //store deletion variant as binary encoded sequences in compact target
        set_BC_val(&ctrg[*ctrg_cnt], &BC_val[*ctrg_cnt], NULL, src_ctrg, MUTANT_BRCD); //set ptr to BC_val structure
        ctrg[*ctrg_cnt].bid = src_ctrg->bid | (*mutCode)++;  //generate and store barcode id
        seq2bin_long(var, &ctrg[*ctrg_cnt].bsq, 1);          //store binary-encoded barcode sequence
        BC_val[*ctrg_cnt].typ = DEL;                         //set barcode type to deletion
        
        //allocate memory for storing character-encoded barcode
        if ((ctrg[*ctrg_cnt].csq = malloc((strlen(var)+1) * sizeof(*ctrg[*ctrg_cnt].csq))) == NULL) {
            printf("mk_DEL_trgts: error - memory allocation for character-encoded target sequence failed. aborting...\n");
            abort();
        }
        strcpy(ctrg[*ctrg_cnt].csq, var); //store character-encoded barcode
        
        (*ctrg_cnt)++;     //increment compact target count
        del_trgts_made++;  //increment number of deletion targets made by one
    }
    
    if (((*mutCode)-1) != MAX_DEL_CODE) {
        printf("mk_DEL_trgts: error - incorrect maximum deletion code (%llu). expected %d. aborting...\n", (long long unsigned int)((*mutCode)-1), MAX_DEL_CODE);
        abort();
    }
    
    return del_trgts_made; //return number of deletion targets made
}


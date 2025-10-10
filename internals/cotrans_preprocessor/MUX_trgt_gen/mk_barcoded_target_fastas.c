//
//  mk_barcoded_target_fastas.c
//  
//
//  Created by Eric Strobel on 10/6/25.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../variant_maker/constant_seqs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/io_management.h"

#include "./mk_MUX_trgts.h"

#include "mk_barcoded_target_fastas.h"

/* mk_barcoded_target_fastas: generate directory that contains a fasta file for each target and a file that contains a list of barcode ids */
void mk_barcoded_target_fastas(TPROBE_names * nm, compact_target * ctrg, target_params * trg_prms)
{
    extern uint64_t trgt_ID_mask; //barcode id mask that isolates the reference barcode id
    extern char sc1[41];             //SC1 leader
    extern char RLA29synch_3p11[34]; //RLA29synch_3p11 linker
    
    uint64_t i = 0;   //general purpose index
    uint64_t nid = 0; //numerical barcode id
    
    FILE * fp_bcid = NULL; //pointer to barcode ids output file
    FILE * fp_trgt = NULL; //pointer for handling output target files
    
    opt_BC * p_opt_BC = NULL; //pointer for handling barcode target optional values
    
    int ret = 0;                    //return value for error checking
    char out_dir[MAX_LINE+1] = {0}; //output directory string
    char out_str[MAX_LINE+1] = {0}; //general purpose output string
    
    //make barcode targets output directory
    ret = snprintf(out_dir, MAX_LINE+1, "%s_barcode_targets", nm->trgts_prfx);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_barcoded_target_fastas: error - error when constructing barcode targets output directory name. aborting...\n");
        abort();
    }
    mk_out_dir(out_dir);
    
    //make barcode ids output file
    ret = snprintf(out_str, MAX_LINE+1, "./%s/%s_barcode_ids.txt", out_dir, nm->trgts_prfx);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_barcoded_target_fastas: error - error when constructing barcode ids output file name. aborting...\n");
        abort();
    }
    if ((fp_bcid = fopen(out_str, "w")) == NULL) {
        printf("mk_barcoded_target_fastas: failed to make barcode_ids.txt file. aborting...\n");
        abort();
    }
    
    for (i = 0; i < trg_prms->t_cnt; i++) {  //for every barcode in the input barcode file
        
        p_opt_BC = (opt_BC *)ctrg[i].opt;                   //set pointer for handling BC_val below
        nid = (ctrg[i].bid & trgt_ID_mask) >> MUTCODE_BITS; //set numerical id
        
        fprintf(fp_bcid, "%05llu\n", (long long unsigned int)nid); //print barcode id to barcode ids file
        
        //make fasta file for current target
        ret = snprintf(out_str, MAX_LINE+1, "./%s/%s_%05llu.fa", out_dir, nm->trgts_prfx, (long long unsigned int)nid);
        if (ret >= MAX_LINE || ret < 0) {
            printf("mk_barcoded_target_fastas: error - error when constructing barcode target fasta file name. aborting...\n");
            abort();
        }
        if ((fp_trgt = fopen(out_str, "w")) == NULL) {
            printf("mk_barcoded_target_fastas: failed to make fasta file for target %05llu. aborting...\n", (long long unsigned int)nid);
            abort();
        }
        fprintf(fp_trgt, ">%s\n", ctrg[i].cid);  //print target name
        fprintf(fp_trgt, "%s\n", p_opt_BC->tsq); //print target sequence
        
        //close current target fasta file
        if (fclose(fp_trgt) == EOF) {
            printf("mk_barcoded_target_fastas: failed to close fasta file for target %05llu. aborting...\n", (long long unsigned int)nid);
            abort();
        }
    }
    
    //close barcode ids file
    if ((fclose(fp_bcid)) == EOF) {
        printf("mk_barcoded_target_fastas: failed to close barcode ids file. aborting...\n");
        abort();
    }
    
    return;
}

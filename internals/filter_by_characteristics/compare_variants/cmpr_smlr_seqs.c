//
//  cmpr_smlr_seqs.c
//  
//
//  Created by Eric Strobel on 6/29/26.
//

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"
#include "../../seq_utils/basemap.h"

#include "../../TECdisplay_mapper/TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper/TECdisplay_mapper_structs.h"

#include "../../TECdisplay_mapper/map_reads/map_standard_TDSPLY_reads.h"
#include "../../TECdisplay_mapper/map_reads/map_expected/get_key.h"
#include "../../TECdisplay_mapper/map_reads/map_expected/set_trgt.h"
#include "../../TECdisplay_mapper/map_reads/map_expected/parse_vmt_trgts.h"
#include "../../TECdisplay_mapper/map_reads/map_expected/print_navigator_template.h"
#include "../../TECdisplay_navigator/parse_reference.h"
#include "../../cotrans_preprocessor/MUX_trgt_gen/set_barcoded_compact_target.h"
#include "../../variant_maker/constant_seqs.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"
#include "../set_seq_attributes/parse_descriptor_file.h"

#include "cmpr_smlr_seqs.h"


/* cmpr_smlr_seqs: perform a systematic comparison of closely-related sequences to identify sequence variants that exhibit large effects on function */
void cmpr_smlr_seqs(char * msa_path, int msa_typ, sequence_attributes * sq_att, int seq_cnt)
{
    if (msa_typ != VMT_FILE) { //if msa is not a VMT file, throw error and abort.
        printf("cmpr_smlr_seqs: error - msa type must be VMT. aborting...\n");
        abort();
    }
    
    FILE * fp_vmt = NULL;           //pointer to vmt file
    
    target_params trg_prms = {0};   //target params for storing target information
    TDSPLY_fasta wt = {0};          //fasta struct for storing wt sequence
    
    target * refs = {NULL};         //pointer for array of reference targets
    opt_ref * ref_val = {NULL};     //pointer for array of optional reference target structures
    
    compact_target * ctrg = NULL;   //pointer for allocating compact_target structures
    opt_mx_trg * trg_val = NULL;    //pointer for allocating target_value structures
    opt_sq_att * sq_att_val = NULL; //pointer for allocating sequence attribute structures
    
    wt_source wt_src = {0};         //wt_source struct for association with basemaps
    basemap * bmap = NULL;          //pointer for allocating basemap structs
    char ** cnstnt_indels = NULL;   //pointer->pointers for storing constant indel strings
    
    int i = 0; //general purpose index
    
    char sub_dir_msa_path[MAX_LINE+1] = {"../"}; //array for generating msa_path from current sub-directory
    
    //generate path to msa file from current subdirectory
    if ((strlen(sub_dir_msa_path) + strlen(msa_path)) <= MAX_LINE) {
        strcat(sub_dir_msa_path, msa_path);
    } else {
        printf("cmpr_smlr_seqs: error - msa path from output subdirectory is too long. aborting...\n");
        abort();
    }
    
    //open msa file
    if ((fp_vmt = fopen(sub_dir_msa_path, "r")) == NULL) {
        printf("cmpr_smlr_seqs: error - failed to open input vmt file. aborting...\n");
        abort();
    }
    
    //allocate memory for reference targets
    if ((refs = calloc(MAXREF, sizeof(*refs))) == NULL) {
        printf("cmpr_smlr_seqs: error - reference target memory allocation failed. aborting...\n");
        abort();
    }
    
    //allocate memory for reference target values
    if ((ref_val = calloc(MAXREF, sizeof(*ref_val))) == NULL) {
        printf("cmpr_smlr_seqs: error - reference target value memory allocation failed. aborting...\n");
        abort();
    }
    
    //allocate memory for compact targets
    if ((ctrg = calloc(seq_cnt, sizeof(*ctrg))) == NULL) {
        printf("cmpr_smlr_seqs: error - target memory allocation failed\n");
        abort(); //TODO: abort or return error code?
    }
    
    //allocate memory for optional values
    if ((trg_val = calloc(seq_cnt, sizeof(*trg_val))) == NULL) {
        printf("cmpr_smlr_seqs: error - opt values memory allocation failed\n");
        abort(); //TODO: abort or return error code?
    }
    
    //allocate memory for utility values
    if ((sq_att_val = calloc(seq_cnt, sizeof(*sq_att_val))) == NULL) {
        printf("cmpr_smlr_seqs: error - utl values memory allocation failed\n");
        abort(); //TODO: abort or return error code?
    }
    
    //allocate memory for basemaps
    if ((bmap = calloc(MAXREF, sizeof(*bmap))) == NULL) {
        printf("cmpr_smlr_seqs: error - basemap memory allocation failed\n");
        abort(); //TODO: abort or return error code?
    }
    
    //allocate memory for constant indel strings
    if ((cnstnt_indels = calloc(MAXREF, sizeof(*cnstnt_indels))) == NULL) {
        printf("cmpr_smlr_seqs: error - constant indel pointer memory allocation failed\n");
        abort(); //TODO: abort or return error code?
    }
    
    //TODO: need to incorporate pairing constraints when setting basemap
    //parse vmt file header lines and targets, which are stored in compact_target structs
    //then use parsed vmt information to generate a basemap for each reference that was
    //present in the vmt file. then associate sequence_attributes and basemaps with their
    //respective targets
    parse_header_lines(fp_vmt, &trg_prms, &wt);
    parse_vmt_trgts(fp_vmt, msa_typ, refs, ref_val, ctrg, trg_val, &trg_prms, &wt, SQ_ATT, CMPCT_TRGT_STRUCT);
    set_bmap_from_prsd_vmt(bmap, &wt_src, cnstnt_indels, refs, &wt, &trg_prms);
    set_sq_att_utl(ctrg, sq_att_val, sq_att, bmap, &trg_prms);
    
    //close input file
    if (fclose(fp_vmt) == EOF) { //close file
        printf("cmpr_smlr_seqs: error - failed to close input vmt file. aborting...\n");
        abort();
    }
    
    //make hash table
    compact_h_node **chtbl = NULL;                  //compact target hash table root
    compact_h_node_bank chn_bank = {NULL, NULL, 0}; //bank for compact target hash table nodes
    chn_bank.count = 0;                             //initialize compact hash table node count to zero
    compact_h_node_bank *crrnt_hn_bank = &chn_bank; //pointer for handling hash node bank

    //allocate TABLE_SIZE hash table node pointers
    if ((chtbl = calloc(TABLE_SIZE, sizeof(*chtbl))) == NULL) {
        printf("cmpr_smlr_seqs: error - hash table memory allocation failed. aborting...\n");
        abort();
    }
    //allocate BLOCK_SIZE hash table nodes
    if ((chn_bank.chn = calloc(BLOCK_SIZE, sizeof(*(chn_bank.chn)))) == NULL) {
        printf("cmpr_smlr_seqs: error - hash table node memory allocation failed. aborting...\n");
        abort();
    }
        
    mk_chtbl_TDSPLY(chtbl, &chn_bank, ctrg, refs, &trg_prms); //make hash table
}

/*set_bmap_from_prsd_vmt: use information from parsed vmt file to generate basemap */
//TODO: need to add ability to include user-specified pairing constraints
void set_bmap_from_prsd_vmt(basemap * bmap, wt_source * wt_src, char ** cnstnt_indels, target * refs, TDSPLY_fasta * wt, target_params * trg_prms)
{
    FILE * fp_ntHDR = NULL; //pointer for generating navigator template header file
    
    int i = 0;   //general purpose index
    int ret = 0; //return value for snprintf
    
    char fn[MAX_LINE+1] = {0}; //array for storing filename of navigator template header file
        
    for (i = 0; i < trg_prms->r_cnt; i++) { //for every reference target
        
        //generate a filename for the navigator template header file
        snprintf(fn, MAX_LINE, "%s_navtmp_ref.txt", refs[i].id);
        if (ret >= MAX_LINE || ret < 0) {
            printf("set_bmap_from_prsd_vmt: error - error when navtmp reference file name. aborting...\n");
            abort();
        }
        
        //open the navigator template header file
        if ((fp_ntHDR = fopen(fn, "w")) == NULL) {
            printf("set_bmap_from_prsd_vmt: failed to make navigator template header file for reference %d. aborting...\n", i);
            abort();
        }
        
        //print navigatore template header to file
        print_nav_tmp_hdr(fp_ntHDR, (opt_ref *)refs[i].opt, wt);
        
        //close navigator template output file
        if ((fclose(fp_ntHDR)) == EOF) {
            printf("set_bmap_from_prsd_vmt: error - error occurred when closing the navigator template file. Aborting program...\n");
            abort();
        }
        
        //open navigator template file in read mode
        if ((fp_ntHDR = fopen(fn, "r")) == NULL) {
            printf("set_bmap_from_prsd_vmt: failed to open navigator template header file for reference %d. aborting...\n", i);
            abort();
        }
        
        //generate basemap from navigatore template header
        set_wt_seq(wt_src, wt->nm, wt->sq);
        parse_reference(fp_ntHDR, bmap, wt_src, &cnstnt_indels[i], 1);
    
        //close navigator template file
        if ((fclose(fp_ntHDR)) == EOF) {
            printf("print_navigator template: error - error occurred when closing the navigator template file. Aborting program...\n");
            abort();
        }
    }
    
    return;
}

/* set_sq_att_utl: associate opt_sq_att struct with compact_target struct via the utility member and set values within opt_sq_att struct */
void set_sq_att_utl(compact_target * ctrg, opt_sq_att * utl, sequence_attributes * sq_att, basemap * bmap,  target_params * trg_prms)
{
    int r = 0; //reference index
    int v = 0; //absolute variant index
    int i = 0; //targets per reference index
    
    opt_sq_att * p_utl = NULL; //pointer for handling opt_sq_att structs
    
    for (r = 0, v = 0; r < trg_prms->r_cnt && r < MAXREF; r++) { //for every reference
        
        for (i = 0; i < trg_prms->t_per_r[r]; i++, v++) { //for every target associated with the current reference
                    
            if (!strcmp(sq_att[v].nm, ctrg[v].cid) && !strcmp(sq_att[v].sq1, ctrg[v].csq)) { //if id/seq in ctrg match those in sq_att
                                
                ctrg[v].utl = &utl[v];               //associate current opt_sq_att struct with ctrg utility member
                
                p_utl = (opt_sq_att *)(ctrg[v].utl); //dereference utility member as opt_sq_att struct
                p_utl->sq_att = &sq_att[v];          //associate current sequence_attributes struct with current opt_sq_att struct
                
                p_utl->bmap = &bmap[r];              //associate basemap for current reference with current opt_sq_att struct

            } else {
                //throw error if sequence attribute name/id does match that of the compact target struct
                printf("assign_sq_att_2_utl: error - sequence attribute name and/or sequence does not match that of the target to which it is assigned:\n\ntrg: %s %s\natt: %s %s\n\naborting...\n", ctrg[v].cid, ctrg[v].csq, sq_att[v].nm, sq_att[v].sq1);
                abort();
            }
        }
            
        return;
    }
}

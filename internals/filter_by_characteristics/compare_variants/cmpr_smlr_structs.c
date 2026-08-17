//
//  cmpr_smlr_structs.c
//  
//
//  Created by Eric Strobel on 6/29/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/data_dist.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"
#include "../set_seq_attributes/parse_descriptor_file.h"

#include "../set_seq_attributes/parse_ct_file.h"

#include "./mk_dot_bracket_htbl.h"
#include "./parse_comparison_file.h"

#include "cmpr_smlr_structs.h"

/* cmpr_smlr_structs: perform a systematic comparison of variants with closely-related structures to identify variants that exhibit large effects on function */
void cmpr_smlr_structs(sequence_attributes * sq_att, descriptor * des, int seq_cnt, int des_cnt, char * des_nm, comparison_values * cmp)
{
    /*
     a hash_table for looking up variants that exhibit a particular secondary structure is generated as follows:
     
     1. in cmpr_smlr_structs:
        - allocate compact_target, opt_db, and sp_list pointer memory for every
          structProps associated with the descriptor being used for comparisons

     2. in set_db_compact_target:
        - generate a compact_target for every structProps structure
        - associate each structProps structure with its corresponding compact_target
        - add each structProps structure to sp_list

     3. in mk_dot_bracket_htbl:
        - map every compact_target (generated from structProps) to the hash table.
          -the number of subsequence compact_targets that match the first instance of a structure are using mul
           (in the first instance ctrg) and blacklisted so that only the first instance of a structure is
           considered during downstream processing
     
        - using init_opt_db, structProps memory associated with an opt_db struct is allocated for every
          structProps structure that matches the first instance of a structure.
     
        - using map_structProps_2_htbl, each structProps is associated with its corresponding compact_target
          such that all variants that exhibit a given structures are associated by the compact target
     */
    
    
    compact_target *ctrg = NULL;    //pointer for allocating compact_target structs
    opt_db * db_val = NULL;         //pointer for allocating opt_db structs
    
    structProps ** sp_list = NULL;  //pointer for generating structProps structure list
    
    int ctrg_cnt = 0;         //number of compact targets
    
    int i = 0;                //general purpose index
    int fnd_des = 0;          //flag that descriptor to be assessed was found
    int des_ix = 0;           //index of descriptor to be assessed
    int des_typ = TYPE_INIT;  //type of descriptor to be assessed
    
    //search array of descriptors to identify match to des_nm, which is the descriptor to be assessed
    for (i = 0, fnd_des = 0, des_ix = 0; i < des_cnt && !fnd_des; i++) {
        if (!strcmp(des[i].nm, des_nm)) { //if a match to des_nm is found
            ctrg_cnt = des[i].tot_sp;     //set compact target count to the total number of structProps associated with the descriptor
            des_ix = i;                   //set des_ix to index of the current descriptor match
            des_typ = des[i].typ;         //set des_typ to the type of the current descriptor
            fnd_des = 1;                  //set flag that match to current descriptor was found
        }
    }
    
    if (!fnd_des) { //if no match to des_nm was found, throw error and abort
        printf("cmpr_smlr_structs: error - could not find descriptor with name '%s'. aborting...\n", des_nm);
        abort();
    }
    
    printf("cmpr_smlr_structs: tot_sp=%d\n", ctrg_cnt);
    
    //allocate memory for compact targets
    if ((ctrg = calloc(ctrg_cnt, sizeof(*ctrg))) == NULL) {
        printf("cmpr_smlr_structs: error - target memory allocation failed\n");
        abort(); //TODO: abort or return error code?
    }
    
    //allocate memory for optional values
    if ((db_val = calloc(ctrg_cnt, sizeof(*db_val))) == NULL) {
        printf("cmpr_smlr_structs: error - opt values memory allocation failed\n");
        abort(); //TODO: abort or return error code?
    }
    
    //allocate memory for structProps list
    if ((sp_list = calloc(ctrg_cnt, sizeof(*sp_list))) == NULL) {
        printf("cmpr_smlr_structs: error - structProps list memory allocation failed\n");
        abort(); //TODO: abort or return error code?
    }
    
    //generate compact_targets from structProps structures
    set_db_compact_target(ctrg, db_val, sp_list, ctrg_cnt, sq_att, seq_cnt, des_ix, des_typ);
    
    /********* hash table initialization and construction **********/
    compact_h_node **htbl = NULL;                  //hash table root
    compact_h_node_bank hn_bank = {NULL, NULL, 0}; //bank for hash table nodes
    
    hn_bank.count = 0; //initialize hash table node count to zero
    
    //allocate TABLE_SIZE hash table node pointers
    if ((htbl = calloc(TABLE_SIZE, sizeof(*htbl))) == NULL) {
        printf("cmpr_smlr_structs: error - hash table memory allocation failed. aborting...\n");
        abort();
    }
    
    //allocate BLOCK_SIZE hash table nodes
    if ((hn_bank.chn = calloc(BLOCK_SIZE, sizeof(*(hn_bank.chn)))) == NULL) {
        printf("cmpr_smlr_structs: error - hash table node memory allocation failed. aborting...\n");
        abort();
    }
    
    //generate htbl using dot-bracket structures
    compact_h_node_bank *crrnt_hn_bank = &hn_bank;    //pointer for handling hash node bank
    mk_dot_bracket_htbl(htbl, crrnt_hn_bank, ctrg, ctrg_cnt, sp_list/*, target_params * trg_prms*/);
    set_sp_z_scores(ctrg, ctrg_cnt, cmp);
    find_ntbl_var_prs(ctrg, ctrg_cnt, cmp);
    
    return;
}

/* set_db_compact_target: generate compact_targets from structProps structures and associate the structProps struct with the compact_target */
void set_db_compact_target(compact_target * ctrg, opt_db * db_val, structProps ** sp_list, int ctrg_cnt, sequence_attributes * sq_att, int seq_cnt, int des_ix, int des_typ)
{
    int i = 0;      //general purpose index
    int j = 0;      //general purpose index
    int c = 0;      //compact_target index
    int sp_cnt = 0; //number of structProps associated with the current attribute
    
    proximal_deltaG * p_prx_dG = NULL; //pointer for dereferencing proximal_deltaG attributes
    distal_deltaG * p_dst_dG = NULL;   //pointer for dereferencing distal_deltaG attributes
    structProps * p_sp = NULL;         //pointer to current structProps structure
        
    for (i = 0, c = 0; i < seq_cnt; i++) { //for every variant sequence
        
        //set values for structProps of current variant sequence
        switch (des_typ) {
            case PRX_DG:                                             //descriptor is PRX_DG
                p_prx_dG = (proximal_deltaG *)sq_att[i].att[des_ix]; //dereference attribute as PRX_DG
                p_sp = &p_prx_dG->sp;                                //set pointer to root structProps
                sp_cnt = p_prx_dG->sp_cnt;                           //set number of structProps associated with current attribute
                break;
                
            case DST_DG:                                             //descriptor is DST_DG
                p_dst_dG = (distal_deltaG *)sq_att[i].att[des_ix];   //dereference attribute as DST_DG
                p_sp = &p_dst_dG->sp;                                //set pointer to root structProps
                sp_cnt = p_dst_dG->sp_cnt;                           //set number of structProps associated with current attribute
                break;
                
            default: //unrecognized descriptor type, throw error and abort
                printf("set_db_compact_target: error - unrecognized descriptor type. aborting...\n");
                abort();
                break;
        }
    
        //generate compact_target from current structProps and
        //and associate current structProps with the compact_target
        for (j = 0; j < sp_cnt; j++) {    //for every structProps associated with the current attribute
                        
            if (j > 0) {                  //if not processing first structProps
                if (p_sp->nxt != NULL) {  //if next pointer to next structProps in linked list is not null
                    p_sp = p_sp->nxt;     //set p_sp to next structProps struct
                    
                } else { //fewer than sp_cnt structProps in linked list, throw error and abort
                    printf("set_db_compact_target: error - fewer than expected structProps structures in the linked list. aborting...\n");
                    abort();
                }
            }
            
            //generate compact_target using dot-bracket structure from current structProps
            if (c < ctrg_cnt) { //check that current ctrg is not exceed allocated ctrgs
                
                //allocate memory for storing dot-bracket sequence in compact_target structure
                if ((ctrg[c].csq = malloc((strlen(p_sp->db)+1) * sizeof(*ctrg[c].csq))) == NULL) {
                    printf("set_db_compact_target: error - failed to allocate compact target char seq memory. aborting...\n");
                    abort();
                }
                
                strcpy(ctrg[c].csq, p_sp->db);                                  //store dot-bracket structure string in compact_target
                dot_bracket_to_bin_hash(ctrg[c].csq, &ctrg[c].bsq, MAX_BLOCKS); //convert dot-bracket structure to binary
                ctrg[c].opt = &db_val[c];                                       //point opt to current opt_db struct
                sp_list[c] = p_sp;                                              //add current structProps to sp_list
                c++;                                                            //increment ctrg index
                
            } else { //not enough compact_targets, throw error and abort
                printf("set_db_compact_target: error - not enough compact targets to store all structure properties entries. aborting...\n");
                abort();
            }
        }
    }
    
    return;
}

/* set_z_scores: calculate z-scores for TECdisplay comparison values */
void set_sp_z_scores(compact_target * ctrg, int ctrg_cnt, comparison_values * cmp)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    double mn = 0; //mean
    double sd = 0; //standard deviation
    
    opt_db * p_db_val = NULL; //pointer for dereferencing opt_db values
    
    double * vals = NULL; //array for storing values for mean/stddev calculations
    
    for (i = 0; i < ctrg_cnt; i++) { //for every target
                
        mn = 0; //reset mean to zero
        sd = 0; //reset stddev to zero
        
        if (!ctrg[i].bl) { //if the target was not blacklisted
            
            //allocate memory for storing values to be used when calculating mean/stddev
            if ((vals = malloc(ctrg[i].cnt * sizeof(*vals))) == NULL) {
                printf("set_sp_z_scores: error - failed to allocate memory for storing comparison values. aborting...\n");
                abort();
            }
            
            p_db_val = ctrg[i].opt; //dereference opt_db vals to simplify code below
            
            //copy values to vals array
            for (j = 0; j < ctrg[i].cnt; j++) {
                vals[j] = p_db_val->sp[j]->sq_att->td_vals[cmp->ix];
            }
            
            mn = mean(vals, ctrg[i].cnt);               //calculate mean
            sd = stddev(vals, ctrg[i].cnt, POPULATION); //calculate stddev
            
            //calculate a z-score for every value
            for (j = 0; j < ctrg[i].cnt; j++) {
                p_db_val->sp[j]->z = calc_Zscore(vals[j], mn, sd);
            }
            
            free(vals);  //free memory
            vals = NULL; //reset vals pointer to NULL
        }
    }
    
    return;
}

/* find_ntbl_var_prs: find variant pairs that exhibit closely related structures but substantial functional differences */
void find_ntbl_var_prs(compact_target * ctrg, int ctrg_cnt, comparison_values * cmp)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
    int l = 0; //general purpose index
        
    structProps ** sp_srtd = NULL; //ptr-to-ptr for generating sorted structProps array
    opt_db * p_db_val = NULL;      //pointer for dereferencing opt_db values
    
    sequence_attributes * sq_att1 = NULL; //sequence attributes comparison pointer 1
    sequence_attributes * sq_att2 = NULL; //sequence attributes comparison pointer 2
    
    double min_z_dif = 5.75;    //minimum z-score difference TODO: make function argument
    double max_dG_dif = 0.5; //maximum  deltaG difference TODO: make function argument
    int fnd_bnd = 0;         //flag that max_dG_dif bound was found
        
    mct_diffs difs = {0}; //structure for tracking min_con_table differences
    
    FILE * ofp = NULL; //output file pointer
    
    //open output file for storing variant pair comparisons
    if ((ofp = fopen("difs.txt", "w")) == NULL) {
        printf("find_ntbl_var_prs: error - failed to open output file. aborting...\n");
        abort();
    }
    
    //print related structures
    for (i = 0; i < ctrg_cnt; i++) { //for every target
        if (!ctrg[i].bl) {           //if the target was not blacklisted
            
            //allocate structProps pointers for sorting
            if ((sp_srtd = calloc(ctrg[i].cnt, sizeof(*sp_srtd))) == NULL) {
                printf("find_ntbl_var_prs: error - failed to allocate memory for sorted structProps array\n");
                abort();
            }
            
            p_db_val = ctrg[i].opt; //dereference opt_db vals to simplify code below
            
            //copy structProps pointers from current target to sp_srtd
            for (j = 0; j < ctrg[i].cnt; j++) {
                sp_srtd[j] = p_db_val->sp[j];
            }
            
            //sort sp_srtd structProps pointers by deltaG
            qsort(sp_srtd, ctrg[i].cnt, sizeof(*sp_srtd), sp_cmpfnc);

            //print header to output file TODO: right now this prints for each variant
            fprintf(ofp, "id\tseq\tdb_an\tdG\tTECdisplay\tz-score\tdelta z-score\tnp_tot\tn_sub\tp_tot\tp_sub\tp_swp\n");
            
            //identify variant pairs in which the delta G difference is less than or equal to max_dG_dif and
            //the z-score difference is greater than or equal to min_z_diff
            for (j = 1; j < ctrg[i].cnt; j++) {
                
                //search upstream of the current variant until the delta G difference is greater than max_dG_dif
                //which causes fnd_bnd to be set to 1
                for (k = 1, fnd_bnd = 0; !fnd_bnd && j-k >= 0; k++) {
                    
                    //test whether the delta G difference of the current two variants is greater than max_dG_dif
                    //note: array is sorted so fabs is not needed
                    if (sp_srtd[j]->dG - sp_srtd[j-k]->dG <= max_dG_dif) {
                        
                        //delta G difference is less than max_dg_dif
                        //test z-score difference
                        if (fabs(sp_srtd[j]->z - sp_srtd[j-k]->z) >= min_z_dif) {
                            
                            //z-score difference is greater than or equal to min_z_dif
                            
                            //zero mct_difs structure
                            difs.np_tot = difs.n_sub = difs.p_tot = difs.p_sub = difs.p_swp = 0;
                            
                            //set sequence attributes pointers
                            sq_att1 = sp_srtd[j]->sq_att;
                            sq_att2 = sp_srtd[j-k]->sq_att;
                                                        
                            //determine differences between the structures of the current two variants
                            cnt_dif_bps(&difs, &sp_srtd[j]->mct, &sp_srtd[j-k]->mct);
                            
                            //print data to output file
                            fprintf(ofp, "%s\t%s\t%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\n", sq_att1->nm, sq_att1->sq1, sp_srtd[j]->db_an, sp_srtd[j]->dG, sq_att1->td_vals[cmp->ix], sp_srtd[j]->z, fabs(sp_srtd[j]->z - sp_srtd[j-k]->z), difs.np_tot, difs.n_sub, difs.p_tot, difs.p_sub, difs.p_swp);
                            fprintf(ofp, "%s\t%s\t%s\t%f\t%f\t%f\n\n", sq_att2->nm, sq_att2->sq1, sp_srtd[j-k]->db_an, sp_srtd[j-k]->dG, sq_att2->td_vals[cmp->ix], sp_srtd[j-k]->z);
                        }
                    } else {         //deltaG difference is greater than max_dG_dif
                        fnd_bnd = 1; //set fnd_bnd to 1 to end loop
                    }
                }
            }
            free(sp_srtd);  //free memory
            sp_srtd = NULL; //reset sp_srtd to NULL
        }
    }
    
    fclose(ofp); //close output file
    
    return;
}

/* cmpfunc: simple float comparison function used to sort structProps array. */
int sp_cmpfnc(const void * a, const void * b)
{
    structProps * sp_a = *(structProps **)(a); //dereference a as structProps pointer
    structProps * sp_b = *(structProps **)(b); //dereference b as structProps pointer
    
    if (sp_a->dG - sp_b->dG > 0) {
        return 1;
        
    } else if (sp_a->dG - sp_b->dG < 0) {
        return -1;
        
    } else {
        return 0;
    }
}

/* cnt_dif_bps: count differences between two related structures */
void cnt_dif_bps(mct_diffs * difs, min_con_table * mct1, min_con_table * mct2)
{
    int i = 0; //general purpose index
        
    //check that connectivity tables are not discordant
    for (i = 0; i < mct1->len; i++) {
        if (mct1->p2ix[i] != mct2->p2ix[i]) {
            printf("cnt_dif_bs: error - cannot count differences in discordant connectivity tables. aborting...\n");
            abort();
        }
    }
    
    char b1 = 0; //base 1
    char p1 = 0; //base paired to base 1
    
    char b2 = 0; //base 2
    char p2 = 0; //base paired to base 2
    
    int * tstd = NULL; //pointer for allocating array to track which nucleotides have been tested
    
    //allocate tstd array memory
    if ((tstd = calloc(mct1->len, sizeof(*tstd))) == NULL) {
        printf("ct_dif_bps: error - failed to allocate array for tracking tested nts. aborting...\n");
        abort();
    }
    
    for (i = 0; i < mct1->len; i++) { //for every position
        
        if (!tstd[i]) {                            //if the current position has not yet been tested
            if (mct1->p2ix[i] == -1) {             //if the current position is unpaired
                if (mct1->sq[i] != mct2->sq[i]) {  //track whether the base at the unpaired position differs
                    difs->n_sub++;
                }
                tstd[i] = 1;                       //mark current position as tested
                
            } else {
                
                //set b1/p1 and b2/p2 to simplify the code below
                b1 = mct1->sq[i];
                p1 = mct1->sq[mct1->p2ix[i]];
                         
                b2 = mct2->sq[i];
                p2 = mct2->sq[mct2->p2ix[i]];
                         
                         
                if (b1 != b2 || p1 != p2) {      //if the base pairs are different
                    if (b1 == p2 && p1 == b2) {  //check whether the orientation of the pair was swapped
                        difs->p_swp++;           //if so, track as swapped pair
                    } else {
                        difs->p_sub++;           //if not, track as substitution pair
                    }
                    difs->p_tot++; //track totaly differences
                }
                
                tstd[i] = tstd[mct1->p2ix[i]] = 1; //mark current position/paired to position as tested
            }
        }
    }
    
    for (i = 0; i < mct1->len; i++) {
        if (!tstd[i]) {
            printf("cnt_dif_bps: error - detected untested positions. aborting...\n");
            abort();
        }
    }
    
    difs->np_tot = difs->n_sub + difs->p_tot; //calculate total differences including nucleotides and pairs
    
    return;
}

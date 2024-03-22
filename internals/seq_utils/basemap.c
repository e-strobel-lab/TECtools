//
//  basemap.c
//  
//
//  Created by Eric Strobel on 5/30/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/test_possible_pairs.h"
#include "../seq_utils/seq2bin_hash.h"

#include "basemap.h"

//TODO: add check that number of variable bases does not exceed seq2bin max key

/* set_wt_seq: initialize wt sequence structure */
void set_wt_seq(wt_source * wt, char * nm, char * sq)
{
    //allocate memory for wt sequence name
    if ((wt->nm = malloc((strlen(nm)+1) * sizeof(*(wt->nm)))) == NULL) {
        printf("set_wt_seq: error - memory allocation for sequence name failed. aborting...\n");
        abort();
    }
    
    //allocate memory for wt sequence
    if ((wt->sq = malloc((strlen(sq)+1) * sizeof(*(wt->sq)))) == NULL) {
        printf("set_wt_seq: error - memory allocation for source sequence failed. aborting...\n");
        abort();
    }
    
    //allocate memory for wt reference position
    if ((wt->pos = calloc((strlen(sq)+1), sizeof(*(wt->pos)))) == NULL) {
        printf("set_wt_seq: error - memory allocation for reference position failed. aborting...\n");
        abort();
    }
    
    //allocate memory for wt lookup array
    if ((wt->lkp = calloc((strlen(sq)+1), sizeof(*(wt->lkp)))) == NULL) {
        printf("set_wt_seq: error - memory allocation for wt lookup array failed. aborting...\n");
        abort();
    }
    
    strcpy(wt->nm, nm); //copy sequence name to wt sequence struct
    
    //determine reference sequence position numbers, check that the
    //wt sequence contains only DNA bases or spacer (.) characters
    //and store the wt sequence
    int i = 0;       //general purpose index
    int j = 0;       //general purpose index
    int ref_pos = 0; //position within wt reference sequence
    
    for (i = 0, ref_pos = 0; sq[i] && i < MAXLEN; i++) { //for every base of the wt sequence
        if (isDNAbase(sq[i])) {         //if the nucleotide is DNA base
            wt->sq[i] = toupper(sq[i]); //set wt sequence as uppercase
            ref_pos++;                  //increment the reference position
            wt->pos[i] = ref_pos;       //store the reference position for the current index
            
        } else if (sq[i] == '.') {      //if the nucleotide is a spacer character
            wt->sq[i] = toupper(sq[i]); //set wt sequence
            wt->pos[i] = 0;             //set the reference position to zero
            if (!ref_pos) {             //if the first nucleotide of the wt sequence has not been found yet
                wt->start_pos--;        //decrement the starting position of the wt sequence
            }                           //^^this determines the start of negative numbered leader positions
            
        } else { //unrecognized DNA base error
            printf("set_wt_seq: error - source sequence contains unrecognized DNA base %c. aborting...\n", wt->sq[i]);
            abort();
        }
    }
    wt->sq[i] = '\0'; //terminate wild-type sequence string
    
    //check that wt sequence does not exceed maximum sequence length
    if (i == MAXLEN && sq[i]) {
        printf("set_wt_seq: error - source sequence length (%d) exceeds maximum sequence length (%d). aborting...\n", i, MAXLEN);
    }
    
    //set negative numbered leader positions in the
    //wt reference sequence position array
    for (i = 0, j = wt->start_pos; j < 0; i++, j++) { //starting at the first negative numbered leader position
        wt->pos[i] = j;                               //set the reference position, and increment by one until
    }                                                 //j = 0, which indicates the end of negative positions

    //set wt sequence lookup pointers
    set_lookup(&wt->lkp0, &wt->lkp_start, wt->lkp, wt->pos, wt->sq);
    
    return;
}

/* init_vtmp_mem: initialize memory in basemap struct for storing variant template */
int init_vtmp_mem(basemap * bmap, int nm_len, int sq_len)
{
    int i = 0; //general purpose index
    
    //for now, allocate everything as MAXLEN. may switch to non-constant
    //dynamic allocation in the future, but this is fine for now.
    
    if ((bmap->nm = calloc(MAXLEN+1, sizeof(*(bmap->nm)))) == NULL) {     //allocate memory for variant name
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    if ((bmap->rS = calloc(MAXLEN+1, sizeof(*(bmap->rS)))) == NULL) {     //allocate memory for variant sequence
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    if ((bmap->nts = calloc(MAXLEN+1, sizeof(*(bmap->nts)))) == NULL) {   //allocate memory for variable base array
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }

    if ((bmap->lkp = calloc(MAXLEN+1, sizeof(*(bmap->lkp)))) == NULL) {   //allocate memory for variable base array
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    if ((bmap->pos = calloc(MAXLEN+1, sizeof(*(bmap->pos)))) == NULL) {   //allocate memory for variable base array
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    if ((bmap->typ = calloc(MAXLEN+1, sizeof(*(bmap->typ)))) == NULL) {   //allocate memory for variable base array
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    //initializing all pairing memory here even if it's not used later. it's
    //not that much memory and is easier to keep initialization all in one place
    for (i = 0; i < MAX_SEC_STRUCT; i++) {
        if ((bmap->rP[i] = calloc(MAXLEN+1, sizeof(*(bmap->rP)))) == NULL) {     //allocate memory for pairing constraints
            printf("init_vtmp_mem: error - variant template memory allocation failed\n");
            abort();
        }
        
        if ((bmap->prs[i] = calloc(MAXLEN+1, sizeof(*(bmap->prs)))) == NULL) {   //allocate memory for pair map array
            printf("init_vtmp_mem: error - variant template memory allocation failed\n");
            abort();
        }
    }
    
    if ((bmap->p_vb = calloc(MAXLEN+1, sizeof(*(bmap->p_vb)))) == NULL) { //allocate memory for variable base table pointers
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    if ((bmap->v_ins = calloc(MAXLEN+1, sizeof(*(bmap->v_ins)))) == NULL) { //allocate memory for constant insertions
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    if ((bmap->c_ins = calloc(MAXLEN+1, sizeof(*(bmap->c_ins)))) == NULL) { //allocate memory for constant insertions
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    if ((bmap->dels = calloc(MAXLEN+1, sizeof(*(bmap->dels)))) == NULL) { //allocate memory for deletions
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    if ((bmap->vb_ix = calloc(MAXLEN+1, sizeof(*(bmap->vb_ix)))) == NULL) { //allocate memory for variable base indices
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    if ((bmap->vi_ix = calloc(MAXLEN+1, sizeof(*(bmap->vi_ix)))) == NULL) { //allocate memory for variable insertion indices
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    if ((bmap->ci_ix = calloc(MAXLEN+1, sizeof(*(bmap->ci_ix)))) == NULL) { //allocate memory for constant insertion indices
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    if ((bmap->d_ix = calloc(MAXLEN+1, sizeof(*(bmap->d_ix)))) == NULL) { //allocate memory for deletion indices
        printf("init_vtmp_mem: error - variant template memory allocation failed\n");
        abort();
    }
    
    return 1;
}

/* set_lookup: set wt sequence positional lookup pointers */
void set_lookup(char *** lkp0, char *** lkp_start, char ** lkp, int * wt_pos, char * sq)
{
    int i = 0; //general purpose index
        
    if (wt_pos[0] < 0) {              //if wt seq has spacers for a leader sequence...
        *lkp0 = &lkp[abs(wt_pos[0])]; //lkp0 index in lkp is the abs value of the first leader position...
        *lkp_start = &lkp[0];         //and lkp_start points to lkp[0], which is the start of the sequence string
    } else if (wt_pos[0] == 1){       //if wt sequence does not have a leader sequence...
        *lkp0 = &lkp[0];              //lkp0 points to lkp[0] and...
        *lkp_start = &lkp[1];         //lkp_start points to lkp[1], which is the start of the sequence string
    } else {                          //else throw error
        printf("set_lookup: unexpected starting wt sequence position (%d). aborting...\n", wt_pos[0]);
        abort();
    }
    
    for (i = 0; sq[i] && i < MAXLEN; i++) { //for every sequence character
        if (wt_pos[i]) {                    //if the character has a wt sequence position value
            (*lkp0)[wt_pos[i]] = &sq[i];    //point the lookup pointer indexed at that position to the character
        }
    }
    
    return;
}


/* get_vbases: parses input variant template string and points variable bases to vbase table entry */
void get_vbases(struct basemap *bmap, int mode)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int vb_cnt = 0;
    
    int set_pos1 = 0;         //flag that wt seq position 1 has been reached
    int rdng_ldr = 0;         //flag that leader nts are currently being read
    int last_ref_pos = 0;     //last observed reference position
    int ins_contig = 0;       //current number of contiguous insertions
    
    insertion * crrnt_ins = NULL; //pointer for setting insertion structure values
    
    for (i = 0, vb_cnt = 0; bmap->rS[i] && i < MAXLEN; i++) {
        
        //track number of variable bases in reference sequences
        if (isIUPACbase(bmap->rS[i]) && !isDNAbase(bmap->rS[i])) {
            vb_cnt++;
        }
        
        /* set basemap reference sequence and variable nucleotide array. */
        switch (toupper(bmap->rS[i])) {
                
            /* constant bases are set as constants in the basemap reference
             sequence and variable nucleotide arrays and have type CONSTANT */
                
            case 'A': bmap->nts[i] = 'A'; bmap->typ[i] = CONSTANT; break; //constant A
            case 'T': bmap->nts[i] = 'T'; bmap->typ[i] = CONSTANT; break; //constant T
            case 'U': bmap->nts[i] = 'T'; bmap->typ[i] = CONSTANT; break; //constant U is masked as T
            case 'G': bmap->nts[i] = 'G'; bmap->typ[i] = CONSTANT; break; //constant G
            case 'C': bmap->nts[i] = 'C'; bmap->typ[i] = CONSTANT; break; //constant C
            
            case '.': //constant alignment gap
                if (bmap->wt->sq[i] != '.') { //check that spacer marker corresponds to spacer in wt sequence
                    printf("get_vbases: error - spacer marker in variant template %s corresponds to non-spacer character in wild type sequence. aborting...\n", bmap->nm);
                    abort();
                } else if (rdng_ldr) { //throw error if leader is interrupted by a spacer character
                    printf("get_vbases: error - all leader nucleotides must be contiguous with position 1 of the sequence. there cannot be spacers between leader nucleotides. aborting...");
                    abort();
                } else {
                    bmap->nts[i] = '.';    //set nts char to '.' spacer
                    bmap->typ[i] = SPACER; //set nucleotide type to SPACER
                }
                break;
                
            case '-': //constant deletion
                //there is no limit on the number of deletions
                //the deletions array is large enough to delete
                //the entire wild type sequence
                if (!isDNAbase(bmap->wt->sq[i])) { //check that deletion marker corresponds to deleted wt base
                    printf("get_vbases: error - deletion marker in variant template %s corresponds to non-DNA base in wild type sequence. aborting...\n", bmap->nm);
                    abort();
                } else if (rdng_ldr) { //throw error if leader is interrupted by a deletion character
                    printf("get_vbases: error - all leader nucleotides must be contiguous with position 1 of the sequence. there cannot be deletion markers between leader nucleotides. aborting...");
                    abort();
                } else {
                    bmap->nts[i] = '-';                   //set nts char to '-' deletion
                    bmap->typ[i] = DELETION;              //set nucleotide type to DELETION
                    bmap->dels[i].pos = bmap->wt->pos[i]; //track position of deletion
                    bmap->dels[i].seq = bmap->wt->sq[i];  //track identity of deletion
                    bmap->d_ix[bmap->d_cnt++] = i;        //store index of deletion, increment d_cnt
                }
                break;
                
            /* the identity of each variable base is set in the basemap
             reference sequence and is indicated in the variable nucleotide
             array with a '*' character. the corresponding pointer, p_vb[i],
             is set to point to the entry for the variable base in the
             varbase table, vb_tbl
             */
            case 'N': bmap->nts[i] = '*'; bmap->p_vb[i] = &(vb_tbl[N]); break;
            case 'R': bmap->nts[i] = '*'; bmap->p_vb[i] = &(vb_tbl[R]); break;
            case 'Y': bmap->nts[i] = '*'; bmap->p_vb[i] = &(vb_tbl[Y]); break;
            case 'S': bmap->nts[i] = '*'; bmap->p_vb[i] = &(vb_tbl[S]); break;
            case 'W': bmap->nts[i] = '*'; bmap->p_vb[i] = &(vb_tbl[W]); break;
            case 'M': bmap->nts[i] = '*'; bmap->p_vb[i] = &(vb_tbl[M]); break;
            case 'K': bmap->nts[i] = '*'; bmap->p_vb[i] = &(vb_tbl[K]); break;
            case 'B': bmap->nts[i] = '*'; bmap->p_vb[i] = &(vb_tbl[B]); break;
            case 'D': bmap->nts[i] = '*'; bmap->p_vb[i] = &(vb_tbl[D]); break;
            case 'H': bmap->nts[i] = '*'; bmap->p_vb[i] = &(vb_tbl[H]); break;
            case 'V': bmap->nts[i] = '*'; bmap->p_vb[i] = &(vb_tbl[V]); break;
            default:
                printf("get_vbases: error - unrecognized character %c in variable base reference sequence. aborting...\n", bmap->rS[i]);
                abort();
                break;
        }
        
        if (bmap->nts[i] == '*') {           //if the current nucleotide is a variable base
            bmap->typ[i] = VARIABLE;         //set nucleotide type to VARIABLE
            bmap->vb_ix[bmap->vb_cnt++] = i; //store the index of the variable base
        }
        
        /* check that non-variable positions in the variant template match the wt sequence
         and determine positional information for variable and constant insertions
         */
        if (isIUPACbase(bmap->wt->sq[i])) {  //wt ref seq nt is an IUPAC base (can only be DNA base)
            
            bmap->pos[i] = bmap->wt->pos[i]; //set basemap reference position to that of the current wt base
            
            if (!set_pos1) {                 //if position 1 of wt sequence was not yet set
                if (bmap->pos[i] == 1) {     //check that current position is position 1
                    set_pos1 = 1;            //set flag that position 1 was set
                    rdng_ldr = 0;            //turn off reading leader flag
                } else {
                    printf("get_vbases: error - set position %d before setting position 1. aborting...\n", bmap->pos[i]);
                    abort();
                }
            }

            last_ref_pos = bmap->wt->pos[i]; //set last observed reference position to that of the current wt base
            ins_contig = 0;                  //set insertion contig length to 0
                        
            /* check that variant template is a match to wt seq at non-randomized, non-gap positions.
             all non-randomized, non-gap nucleotides in the variant template must match the wt reference
             sequence at the same position.
             */
            if (isIUPACbase(bmap->nts[i]) &&       //if variant template nt is an IUPAC base and the variant
                bmap->nts[i] != bmap->wt->sq[i]) { //template does not match the wt sequence, throw error
                printf("get_vbases: error - mismatch between variant template sequence (%c) and wt sequence (%c) at a non-randomized, non-gap position (position %d). aborting...\n", bmap->nts[i], bmap->wt->sq[i], bmap->wt->pos[i]);
                abort();
            }
            
        } else if (bmap->wt->sq[i] == '.' && isIUPACbase(bmap->rS[i])) { //nt is part of insertion or leader sequence

            if (!set_pos1) {                        //if reading leader sequence (pos 1 not yet set)
                rdng_ldr = 1;                       //set reading leader flag to true
                last_ref_pos = bmap->wt->pos[i];    //set the last_ref_pos to the current wt seq position
                ins_contig = 0;                     //make sure insertion contig length is 0
            } else {                                //if reading insertion sequence
                ins_contig++;                       //increment insertion contig length
            }
            
            bmap->pos[i] = last_ref_pos;            //set basemap reference position to the last ref position
            
            if (bmap->nts[i] == '*') {
                crrnt_ins = &bmap->v_ins[i];        //set crrnt_ins pointer to variable insertion array index i
                bmap->vi_ix[bmap->vi_cnt++] = i;    //store index of variable insertion, increment vi_cnt
            } else if (isIUPACbase(bmap->nts[i])) {
                crrnt_ins = &bmap->c_ins[i];        //set crrnt_ins pointer to constant insertion array index i
                bmap->ci_ix[bmap->ci_cnt++] = i;    //store index of constant insertion, increment ci_cnt
            } else {
                printf("get_vbases: error - unexpected character in basemap nts array. aborting...\n");
                abort();
            }
            
            crrnt_ins->pos = last_ref_pos; //set basemap reference position to the last ref position
            crrnt_ins->ins_contig = (set_pos1) ? ins_contig : 0; //set ins_contig. leader nts have ins_contig=zero
            crrnt_ins->seq = bmap->rS[i];  //set sequence
        }
    }
    
    bmap->nts[i] = '\0'; //terminate variable nts string
    
    if (vb_cnt > SEQ2BIN_MAX_KEY) {
        printf("get_vbases: variant template %s contains %d variable bases. the maximum number of variable bases is %d\n", bmap->nm, vb_cnt, SEQ2BIN_MAX_KEY);
        abort();
    }
    
    return;
}

/* set_basemap_seq: set basemap structure sequence values */
void set_basemap_seq(basemap * bmap, char * p_nm, char * p_vb, wt_source * wt)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    if (wt->sq == NULL) { //check that wt sequence pointer was set, which indicates wt source struct was set
        printf("set_basemap_seq: error - wt source sequence was not set before setting basemap. aborting...\n");
        abort();
    }
    
    if (strlen(p_nm) > MAXLEN) { //check that name does not exceed maximum length
        printf("set_basemap_seq: error - name length (%lu) exceeds maximum length (%d). aborting...\n", strlen(p_nm), MAXLEN);
        abort();
    }
    
    if (strlen(p_vb) > MAXLEN) { //check that sequence does not exceed maximum length
        printf("set_basemap_seq: error - variable base sequence length (%lu) exceeds maximum length (%d). aborting...\n", strlen(p_vb), MAXLEN);
        abort();
    }
    
    if (strlen(p_vb) != strlen(wt->sq)) { //check that wt and variable bases sequences are the same length
        printf("set_basemap_seq: error - variable base sequence length (%lu) is not equal to wt sequence length (%lu). aborting...\n", strlen(p_vb), strlen(wt->sq));
        abort();
    }
    
    init_vtmp_mem(bmap, strlen(p_nm), strlen(p_vb)); //initialize basemap memory
    bmap->wt = wt;                                   //point basemap wt member to wt source sequence
    strcpy(bmap->nm, p_nm);                          //set basemap name
    
    for (i = 0; p_vb[i] && i < MAXLEN; i++) {  //for every character in the variable base sequence string
        if (isIUPACbase(p_vb[i]) ||            //check that the character is an IUPAC base
            p_vb[i] == '.'       ||            //or a spacer character
            p_vb[i] == '-') {                  //or a deletion character
            bmap->rS[i] = toupper(p_vb[i]);    //copy uppercase variable base sequence to basemap reference sequence
        } else {
            printf("set_basemap_seq: error - variable base sequence contains invalid character %c. aborting...\n", p_vb[i]);
            abort();
        }
    }
    bmap->rS[i] = '\0'; //terminate reference sequence string

    get_vbases(bmap, NO_MASK); //parse variable bases in variable base reference sequence for storage as basemap
    set_lookup(&bmap->lkp0, &bmap->lkp_start, bmap->lkp, bmap->wt->pos, bmap->rS); //construct basemap lookup array
    
    return;
}

/* add_pairs_to_basemap: add pairing constraint to basemap struct */
void add_pairs_to_basemap(basemap * bmap, char * p_pr, char * p_id)
{
    int i = 0; //general purpose index
    
    if (bmap->rP_cnt == MAX_SEC_STRUCT) { //check that max number of pairing constraints will not be exceeded
        printf("add_pairs_to_basemap: error - number of pair constraints for variant template %s exceeds the maximum (%d). aborting...\n", bmap->nm, MAX_SEC_STRUCT);
        abort();
    }
    
    //check that the secondary structure length matches source sequence length
    if (strlen(p_pr) != strlen(bmap->wt->sq)) {
        printf("add_pairs_to_basemap: error - variant template %s pair constraint %s length (%lu) does not match wt sequence length (%lu). aborting...\n", bmap->nm, p_id, strlen(p_pr), strlen(bmap->wt->sq));
        abort();
    }
    
    //check that pair constraint contains only valid characters
    for (i = 0; p_pr[i] && i < MAXLEN; i++) {
        
        //check that pair constraint string char is valid and consistent with reference sequence
        if (p_pr[i] == '.') { //character is a spacer
            
            //do nothing if character is a spacer because spacers can correspond
            //to any character in the variable base reference sequence
            
        } else if (p_pr[i] == '(' || p_pr[i] == ')') { //character specifies a base pair
            
            if (!isIUPACbase(bmap->rS[i])) { //check that corresponding variable base seq char is a base
                printf("add_pairs_to_bmap: error - '%c' character in variant template %s pair constraint %s corresponds to non-base character '%c' in variable base reference sequence. aborting...\n", p_pr[i], bmap->nm, p_id, bmap->rS[i]);
                abort();
            }
            
        } else { //unrecognized character
            printf("add_pairs_to_bmap: error - variant template pair constraint contains unrecognized character '%c'. pair constraints must contain only '.', '(', and ')' characters. aborting...\n", p_pr[i]);
            abort();
        }
    }

    //copy pair constraint to basemap structure, increment rP_cnt
    strcpy(bmap->rP[bmap->rP_cnt++], p_pr);
    
    return;
}

/* set_basemap_pairs: set pair parameters in the basemap struct. */
void set_basemap_pairs(basemap * bmap, int mode)
{
    /* store simple pair parameters in the basemap struct. each time the first member of a simple
     pair is observed (by the presence of a '(' char in the rP seq), store its index in the pair_bfr
     array. when the second member of a pair is observed, (by the presence of a ')' char in the rP
     seq), establish it as a pair with the most recent addition to the pair_bfr array by performing
     the following operations:
     
     1. set the pairs array entry for the first mate to the index of the second mate
     2. set the pairs array entry for the second mate to the index of the first mate
     3. if the sequence is a barcode template, set the variable nucleotides entry for
        the second mate to 'P' to indicate it must pair with the first mate
     */
 
    if (bmap->rS == NULL) { //check that set_basemap_seq was run
        printf("set_basemap_pairs: error - set_basemap_pairs cannot be run unless set_basemap_seq was run first\n. aborting...");
        abort();
    }
    
    if (bmap->rP_cnt < 1) {
        printf("set_basemap_pairs: error - basemap must contain at least 1 reference pair constraint. aborting...\n");
        abort();
    }
    
    if (mode == NO_MASK) { //running NO_MASK mode
        //no checks needed for NO_MASK mode
        
    } else if (mode == MASK_PAIRS) { //running MASK_PAIRS mode
        if (bmap->rP_cnt != 1) {     //requires 1 pair constraint
            printf("set_basemap_pairs: error - basemap pairs cannot be set in MASK_PAIRS mode if the basemap does not contain exactly one pair constraint. aborting...\n");
            abort();
        }
        
    } else { //unrecognized mode
        printf("set_basemap_pairs: unrecognized mode. aborting...\n");
        abort();
    }
    
    int i = 0; //general purpose index
    int p = 0; //index for pair constraints entries
    
    //variables for tracking simple pairs
    int * pair_bfr = NULL;    //buffer for storing location of the first mate in simple pairs
    int bfr_indx = 0;         //index for pair buffer
    int tot_smpl_1st = 0;     //total number of simple pair first mates found
    int tot_smpl_2nd = 0;     //total number of simple pair second mates found
    
    if ((pair_bfr = calloc((strlen(bmap->rS)+1), sizeof(*pair_bfr))) == NULL) { //allocate pair buffer memory
        printf("get_vbases: error - pair buffer memory allocation failed\n");
        abort();
    }
    
    for (p = 0; bmap->rP[p][0] && p < MAX_SEC_STRUCT; p++) { //for every secondary structure constraint string
        
        for (i = 0; bmap->rP[p][i]; i++) {
            if (bmap->rP[p][i] == '(') {  //found first mate of a pair
                tot_smpl_1st++;           //increment count of 1st mates
                pair_bfr[bfr_indx++] = i; //store index of first mate in pair_bfr array
    
            } else if (bmap->rP[p][i] == ')') { //nucleotide is the second mate in a simple pair
                tot_smpl_2nd++;                 //increment count of second mates
                
                if (mode == MASK_PAIRS) {       //in MASK_PAIRS mode, mask the second mate of a pair as 'P' to signal
                    bmap->nts[i] = 'P';         //that the second mate should be set as a WC pair with the first mate
                }
                
                //decrement the pair buffer index to that of the last observed first mate. this
                //entry in the pair buffer will be over written the next time a '(' is observed.
                if (--bfr_indx < 0) { //check that not decrementing beyond array bounds
                    printf("get_vbases: error - found second mate of pair before first mate (pair buffer was empty). aborting...\n");
                    abort();
                }
                
                //check that mates can form a pair and that at least 1 mate is a variable base, then establish
                //the pair. this test implicitly confirms that both mates are IUPAC DNA base characters
                if (test_possible_pairs(bmap->rS[pair_bfr[bfr_indx]], bmap->rS[i], TEST_BASEPAIR, NO_TYPE_TEST) && //mates can pair
                    (bmap->nts[pair_bfr[bfr_indx]] == '*' || bmap->nts[i] == '*')) {                   //>=1 mate is variable
                    
                    bmap->prs[p][pair_bfr[bfr_indx]] = i; //set pr_indx entry for the 1st mate to the index of the 2nd mate
                    bmap->prs[p][i] = pair_bfr[bfr_indx]; //set pr_indx entry for the 2nd mate to the index of the 1st mate

                } else { //invalid pair
                    printf("get_vbases: error - pair between bases %d%c and %d%c is invalid. confirm that these nucleotides are IUPAC DNA characters that can form a pair and that at least 1 is a variable base. aborting...\n", pair_bfr[bfr_indx]+1, bmap->rS[pair_bfr[bfr_indx]], i+1, bmap->rS[i]);
                    abort();
                }
            }
        }
    }
    
    //check that the number of observed first pair mates
    //is equal to the number of observed second pair mates
    if (tot_smpl_1st != tot_smpl_2nd) {
        printf("get_vbases: error - found %d simple pair 1st mate(s) and %d simple pair 2nd mate(s). these values should be equal. aborting...\n", tot_smpl_1st, tot_smpl_2nd);
        abort();
    }
    
    free(pair_bfr); //free pair buffer memory
    
    return;
}

/* mk_vbase_nm: construct variable base name */
void mk_vbase_nm(basemap * bmap, int i, char *vb_nm, int array_size, char vbase)
{
    char tmp[MAXLEN] = {0}; //temp array for storing variable base name
    
    /* leader nucleotides are not considered insertions
     for the purpose of variant names and are
     filtered by their negative position number */
    
    if (vbase && bmap->typ[i] != VARIABLE) {
        printf("mk_vbase_nm: error - variable base was supplied for non-variable position %d. aborting...", i);
        abort();
    }
    
    if (bmap->typ[i] == VARIABLE) {   //nucleotide is a variable base
        
        if (bmap->v_ins[i].pos > 0) { //variable nucleotide is an insertion
            
            //use the format:
            //<position of last non-insertion base>i<contiguous insertion number><base identity>
            
            sprintf(tmp, "%di%d%c", bmap->v_ins[i].pos, bmap->v_ins[i].ins_contig, (vbase) ? vbase : toupper(bmap->rS[i]));
            
        } else { //variable nucleotide is not an insertion
            
            //use the format:
            //<position><base identity>
            
            sprintf(tmp, "%d%c", bmap->pos[i], (vbase) ? vbase : toupper(bmap->rS[i]));
        }
        
    } else if (bmap->typ[i] == CONSTANT) { //nucleotide is a constant base
    
        if (bmap->c_ins[i].pos > 0) {    //constant nucleotide is an insertion
            
            //use the format:
            //c<position of last non-insertion base>i<contiguous insertion number><base identity>
            
            sprintf(tmp, "c%di%d%c", bmap->c_ins[i].pos, bmap->c_ins[i].ins_contig, toupper(bmap->rS[i]));
            
        } else { //constant nucleotide is not an insertion
            
            sprintf(tmp, "c%d%c", bmap->pos[i], toupper(bmap->rS[i]));
        }
        
    } else if (bmap->typ[i] == DELETION) { //nucleotide is a deletion
        
        //use the format:
        //d<position of deletion>i<base identity>
        
        sprintf(tmp, "d%d%c", bmap->dels[i].pos, bmap->dels[i].seq);
        
    } else {
        tmp[0] = '\0';
    }

    if (strlen(tmp) <= MAX_VBASE_NAME_FIELD && array_size >= MAX_VBASE_NAME_FIELD+1) {
        strcpy(vb_nm, tmp);
    }
    
    return;
}

/* print_basemap: print basemap sequence/structure information */
void print_basemap(basemap * bmap) {
    int i = 0; //general purpose index
    int j = 0; //general purpose index

    char tmp[MAX_VBASE_NAME_FIELD+1] = {0}; //temporary array for storing indel position/insertion string
    int field_len = 6;  //size of indel position/insertion field
    int spaces = 0;     //number of spaces to print to fill indes position/insertion field
        
    for (i = 0; bmap->rS[i]; i++) {
        
        //print wt pos/seq, vbase reference pos/seq, and nts code
        printf("%3d%c %3d%c %c ", bmap->wt->pos[i], bmap->wt->sq[i], bmap->pos[i], bmap->rS[i], bmap->nts[i]);
       
        tmp[0] = '\0';
        
        if (bmap->typ == VARIABLE || bmap->v_ins[i].pos || bmap->c_ins[i].pos || bmap->dels[i].pos) {
            mk_vbase_nm(bmap, i, tmp, MAX_VBASE_NAME_FIELD+1, '\0'); //assemble vbase name
            printf("%s", tmp);                                //print vbase name
        }
        
        //fill in indel position/insertion field width with spaces
        for (j = 0, spaces = field_len - strlen(tmp); j < spaces; j++) {
            printf(" ");
        }
        
        //print pair constraint indices
        if (bmap->rP_cnt > 0) {
            for (j = 0; j < bmap->rP_cnt; j++) {
                printf(" |  %c %3d %3d", bmap->rP[j][i], i, bmap->prs[j][i]);
            }
        }
        
        //print type
        switch (bmap->typ[i]) {
            case VARIABLE: printf("\tvariable"); break;
            case CONSTANT: printf("\tconstant"); break;
            case DELETION: printf("\tdeletion"); break;
            case SPACER: printf("\tspacer"); break;
            default:
                printf("\t");
                break;
        }
        
        printf("\n");
    }
    printf("\n");
    
    return;
}

/* print_vbases: print table of variable base information */
void print_vbases(basemap * bmap)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    printf("\nnum     pos     ins     seq  lookup\n");
    for (i = 0, j = 0; bmap->nts[i]; i++) {
        if (bmap->nts[i] == '*') {
            printf("%3d  |  %3d  |  %3d  |   %c  |   %c  |  VBASE:", i, bmap->pos[i], bmap->v_ins[i].ins_contig, bmap->rS[i], bmap->lkp0[bmap->pos[i]][bmap->v_ins[i].ins_contig]);

            if (bmap->wt->pos[i]) {
                printf("substitution\n");
            } else if (bmap->v_ins[i].pos) {
                printf("insertion\n");
            }
            j++;
        }
    }
    printf("\n"); return;
}

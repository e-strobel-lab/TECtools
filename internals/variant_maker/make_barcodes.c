//
//  barcode_variants.c
//  
//
//  Created by Eric Strobel on 3/2/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../utils/io_management.h"
#include "../utils/gen_utils.h"
#include "../seq_utils/ispair.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/basemap.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"
#include "./count_variants.h"
#include "./expand_variant_template.h"

#include "make_barcodes.h"

/* NOTE: the barcode linker sequence is appended in get_rndm_brcds.
 any change to the barcode linker sequence must be made there. */

//barcode storage
barcode_bank bc_root = {0};    //root barcode bank
barcode_bank * bc_crnt = NULL; //pointer to current barcode bank
uint64_t barcode_count = 0;    //number of barcodes that passed filter

/* mk_brcds: coordinates barcode generation */
int mk_brcds(int brcds2mk) {
    
    printf("running in MAKE_BARCODES mode\n");
    
    extern fasta *vrnts;    //pointer to fasta structures used when generating variants
    extern uint64_t v_indx; //index of current variant
    
    extern FILE * prcs_ofp;            //output file pointer for processing messages
    extern char prcs_out_nm[MAX_LINE]; //name of processing message output file
    extern char out_msg[MAX_LINE];     //output message
    
    int i = 0;                  //general purpose index
    uint64_t passed_filter = 0; //number of barcodes that passed the parameter filter
    target * pf_brcds = NULL;   //pointer to target structures, used when picking random barcodes
    
    wt_source brcd_src[BARCODE_TEMPLATES] = {0}; //wt_source array for storing barcode template sequences
    basemap brcd_bmap[BARCODE_TEMPLATES] = {0};  //basemap array for storing barcode template sequences
    
    char out_dir[512] = {0};          //array to store output directory name
    sprintf(out_dir, "barcodes_out"); //construct output directory name
    mk_out_dir(out_dir);              //make output directory
    
    //generate input details file
    sprintf(prcs_out_nm, "./%s/barcodes_processing.txt", out_dir);
    if ((prcs_ofp = fopen(prcs_out_nm, "w")) == NULL) {
        printf("mk_brcds: ERROR - could not open processing messages file. Aborting program...\n");
        abort();
    }
    
    init_brcd_tmplts(&brcd_src[0], &brcd_bmap[0]);   //initialize variant templates for barcode generation
    passed_filter = xpnd_brcd_tmplts(&brcd_bmap[0]); //expand barcode variant templates
    
    printf("passed filter: %llu\n", (long long unsigned int)passed_filter);
    merge_barcode_bank(&bc_root, passed_filter);
    
    //allocate memory to store passed filter barcodes in targets structures
    if ((pf_brcds = calloc(passed_filter, sizeof(*pf_brcds))) == NULL) {
        printf("mk_brcds: error - passed filter barcodes memory allocation failed. aborting...\n");
        abort();
    }
    store_pf_brcds(pf_brcds, passed_filter);  //store passed filter barcodes in targets structures
    
    
    //allocate target structure pointers for output barcodes
    target ** brcd_out = NULL; //pointers to output barcode structures
    if ((brcd_out = calloc(brcds2mk, sizeof(*brcd_out))) == NULL) {
        printf("mk_brcds: error - passed filter barcodes memory allocation failed. aborting...\n");
        abort();
    }
    
    //randomly select barcodes to output
    get_rndm_brcds(brcd_out, pf_brcds, passed_filter, brcds2mk);
    
    //close output file
    if (fclose(prcs_ofp) == EOF) {
        printf("main: error - error occurred when closing processing messages output file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

/* init_brcd_tmplts: initialize variant templates for barcode sequences */
int init_brcd_tmplts(wt_source * brcd_src, basemap * brcd_bmap)
{
    printf("initializing barcode templates\n");
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int vars = 0; //total number of variants encoded by barcode templates
    
    const char brcd_nm[17] = {"barcode_template"}; //name placeholder for barcode variant templates
    const char brcd_sq[BARCODE_TEMPLATES][45] ={   //barcode variant template sequences
                             //NNNNNNNNSloopSNNNNNNNNtKKtKYYgRYgRtcNNloopNN
                              "NNNNNNNNSTKCGSNNNNNNNNTKKTKYYGRYGRTCNNTKCGNN",
                              "NNNNNNNNSTKCGSNNNNNNNNTKKTKYYGRYGRTCNNGKGANN",
                              "NNNNNNNNSGKGASNNNNNNNNTKKTKYYGRYGRTCNNTKCGNN",
                              "NNNNNNNNSGKGASNNNNNNNNTKKTKYYGRYGRTCNNGKGANN"};
    const char brcd_pr[45] = {"(((((((((....)))))))))..............((....))"}; //barcode variant template structure
    
    //initialize a variant template for each barcode sequence
    
    for (i = 0; i < BARCODE_TEMPLATES; i++) {
        
        //allocate memory for sequence name in wt_source struct
        if ((brcd_src[i].nm = malloc((strlen(brcd_nm)+1) * sizeof(*(brcd_src[i].nm)))) == NULL) {
            printf("init_brcd_tmplts: error - memory allocation for sequence name failed. aborting...\n");
            abort();
        }
        
        //allocate memory for source sequence in wt_source struct
        if ((brcd_src[i].sq = malloc((strlen(brcd_sq[i])+1) * sizeof(*(brcd_src[i].sq)))) == NULL) {
            printf("init_brcd_tmplts: error - memory allocation for source sequence failed. aborting...\n");
            abort();
        }
        
        //allocate memory for source sequence positions in wt_source struct
        if ((brcd_src[i].pos = malloc((strlen(brcd_sq[i])+1) * sizeof(*(brcd_src[i].pos)))) == NULL) {
            printf("init_brcd_tmplts: error - memory allocation for source sequence failed. aborting...\n");
            abort();
        }
        
        strcpy(brcd_src[i].nm, brcd_nm);          //set barcode name in wt_source
        strcpy(brcd_src[i].sq, brcd_sq[i]);       //set barcode sequence in wt_source
        for (j = 0; brcd_src[i].sq[j]; j++) {
            if (isIUPACbase(brcd_src[i].sq[j])) { //sanity check
                brcd_src[i].pos[j] = j+1;         //set barcode seq positions numbers
            } else {                              //sanity check failed?
                printf("why did you put a non-IUPAC base character (%c) in the barcode?\n", brcd_src[i].sq[j]);
                abort();
            }
        }
        brcd_bmap[i].wt = &brcd_src[i]; //set basemap wt sequence to point to barcode sequence wt_source
        
        //initialize memory for barcode variant template
        init_vtmp_mem(&brcd_bmap[i], strlen(brcd_nm), strlen(brcd_sq[i]));
        strcpy(brcd_bmap[i].nm, brcd_nm);    //set barcode template name
        strcpy(brcd_bmap[i].rS, brcd_sq[i]); //set barcode template sequence
        strcpy(brcd_bmap[i].rP[0], brcd_pr); //set barcode template pair constraints
        brcd_bmap[i].rP_cnt = 1;             //set pair constraint count to 1
        
        printf("%s\n", brcd_bmap[i].rS);
    }
    
    return 1;
}

/* xpnd brcd_tmplts: expand barcode templates to generate all possible barcode variants */
uint64_t xpnd_brcd_tmplts(basemap * brcd_bmap)
{
    printf("\nexpanding barcode templates\n");
    
    extern barcode_bank bc_root;   //root barcode bank
    extern barcode_bank * bc_crnt; //pointer to current barcode bank
    extern uint64_t barcode_count; //total number of passed filter barcodes
    
    int i = 0;    //general purpose index
    int vars = 0; //total number of variants encoded by barcode templates
    
    char var_outpt[MAXLEN] = {0};       //array to store sequences during barcode expansion
    struct trgt_vb var_lcl_bases = {0}; //array of variable bases that were set for a given barcode
    
    init_barcode_bank(&bc_root, NULL);  //initialize root barcode bank
    bc_crnt = &bc_root;                 //set current barcode bank pointer to root barcode bank
    
    //expand barcode templates
    for (i = 0; i < BARCODE_TEMPLATES; i++) {
        printf("expanding ref %s\n", brcd_bmap[i].rS);
        get_vbases(&brcd_bmap[i], MASK_PAIRS);          //convert variant template to a basemap
        set_basemap_pairs(&brcd_bmap[i], MASK_PAIRS);   //set pairing information in MASK_PAIRS mode
        expand_variant_template(&brcd_bmap[i], 0, var_outpt, var_lcl_bases, MAKE_BARCODES); //expand variant template
    }
    
    printf("\nexpanded variant templates yielding %llu pf barcodes\n", (long long unsigned int)barcode_count);
    
    return barcode_count; //return number of barcodes generated
}

/* filter_barcode: filter barcodes that do not meet the parameters specified in make_barcodes.h */
int filter_barcode(char * sq)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int at_cnt = 0;        //total number of AT nucleotides in barcode sequence
    int gc_cnt = 0;        //total number of GC nucleotides in barcode sequence
    
    int max_homopol = 0;   //maximum observed homopolymer length
    int crnt_homopol = 0;  //variable to track length of a homopolymer
    char prev_base = '\0'; //identity of previous base, used for homopolymer tracking
    
    int s1bot_gc = 0;      //number of GC nucleotides in positions 1 to 3
    int s1top_gc = 0;      //number of GC nucleotides in positions 4 to 6
    int s1all_gc = 0;      //number of GC nucleotides in positions 1 to 6
    int s2bot_gc = 0;      //number of GC nucleotides in positions 27 and 28
    int s2top_gc = 0;      //number of GC nucleotides in positions 36 and 37
        
    for (j = 0; sq[j]; j++) {
        
        //track base identity
        switch (sq[j]) {
            case 'A':
                at_cnt++;
                break;
                
            case 'T':
                at_cnt++;
                break;
                
            case 'G':
                gc_cnt++;
                if (j >= STEM1BOT_N_STRT && j <= STEM1BOT_N_END) {s1bot_gc++;} //track bottom stem1 gc content
                if (j >= STEM1TOP_N_STRT && j <= STEM1TOP_N_END) {s1top_gc++;} //track top stem1 gc content
                if (j >= STEM1BOT_N_STRT && j <= STEM1TOP_N_END) {s1all_gc++;} //track entire stem1 gc content
                if (j >= STEM2BOT_STRT && j <= STEM2BOT_END) {s2bot_gc++;}  //track s2bot gc content
                
                if (j >= STEM2TOP_STRT && j <= STEM2TOP_END) {s2top_gc++;}  //track s2top gc content
                break;
            case 'C':
                gc_cnt++;
                if (j >= STEM1BOT_N_STRT && j <= STEM1BOT_N_END) {s1bot_gc++;} //track bottom stem1 gc content
                if (j >= STEM1TOP_N_STRT && j <= STEM1TOP_N_END) {s1top_gc++;} //track top stem 1 gc content
                if (j >= STEM1BOT_N_STRT && j <= STEM1TOP_N_END) {s1all_gc++;} //track entire stem1 gc content
                if (j >= STEM2BOT_STRT && j <= STEM2BOT_END) {s2bot_gc++;}  //track s2bot gc content
                
                if (j >= STEM2TOP_STRT && j <= STEM2TOP_END) {s2top_gc++;}  //track s2top gc content
                
                break;
                
            default:
                printf("filter_barcodes: error - unrecognized base. aborting...\n");
                abort();
                break;
        }
        
        //track max homopolymer length
        if (sq[j] == prev_base) { //found homopolymer
            crnt_homopol++;       //increment crnt_homopolymer length
            
            //check if current homopolymer exceeds max observed homopolymer
            if (crnt_homopol > max_homopol) {
                max_homopol = crnt_homopol;
            }
            
        } else { //not in homopolymer
            crnt_homopol = 1; //set crnt_homopol to 1
        }
        
        prev_base = sq[j]; //set previous base to current base
    }
    
    //check if barcode sequence meets all specifications
    if (max_homopol > HOMOPOL_LIMIT ||
        gc_cnt < BRCD_GC_MIN        ||
        gc_cnt > BRCD_GC_MAX        ||
        s1bot_gc < S1BOT_GC_MIN     ||
        s1bot_gc > S1BOT_GC_MAX     ||
        s1top_gc < S1TOP_GC_MIN     ||
        s1top_gc > S1TOP_GC_MAX     ||
        s1all_gc < S1ALL_GC_MIN     ||
        s1all_gc > S1ALL_GC_MAX     ||
        s2top_gc > S2TOP_GC_LIMIT /*||
        pos25to26_gc > POS25TO26_GC_LIMIT*/) {
        return 1; //return fail
    } else {
        return 0; //return pass
    }
}

/* store_pf_brcds: store barcodes that passed filter in a targets struct*/
int store_pf_brcds(target * pf_brcds, uint64_t passed_filter)
{
    extern fasta *vrnts; //pointer to fasta structures used when generating variant
    extern uint64_t v_indx;   //index of current variant
    
    int i = 0;           //general purpose index
    int b_indx = 0;      //barcode index
    
    for (i = 0; i < v_indx; i++) { //for every sequence in the vrnts array
        
        //allocate memory for storing barcode sequence in target structure
        if ((pf_brcds[b_indx].sq = malloc((strlen(vrnts[i].sq)+1) * sizeof(*(pf_brcds[b_indx].sq)))) == NULL) {
            printf("store_pf_brcds: error - failed to allocate memory for barcode sequence\n");
            abort();
        }
        
        strcpy(pf_brcds[b_indx].sq, vrnts[i].sq); //copy barcode sequence to target structure
        free(vrnts[i].sq); //free memory for current barcode in vrnts array to keep memory use static
        b_indx++;
    }
     
    //check that number of sequences stored in pf_brcds matches the number of sequences that passed filter
    if (b_indx != passed_filter) {
        printf("store_pf_barcodes: error - variants stored in pf_brcds array (%llu) does not equal the number of passed filter variants(%d). aborting...\n", (long long unsigned int)passed_filter, b_indx);
        abort();
    }
    
    return 1;
    
}

/* get_rndom_brcds: select random passed filter barcodes to output*/
int get_rndm_brcds(target ** brcd_out, target * pf_brcds, uint64_t passed_filter, int brcds2mk)
{
    srand(time(NULL)); //seed random number generation
    
    int i = 0;  //general purpose index
    int j = 0;  //general purpose index
    int r = 0;  //randomized index for selecting barcodes
    
    const char brcd_lnkr[10] = "AAACCAACT"; //barcode linker sequence //TODO: confirm finalized
    
                        //NNNNNNNNSloopSNNNNNNNNtKKtKYYgRYgRtcNNloopNN
    char rPos[45] =     {"*************..........**.***.**.*..******.."}; //randomized positions
    char brcd_sec[45] = {"(((((((((....))))))))).....(((((((((((....))"}; //secondary structure
    //NOTE: the K in the 'KYY' motife at the beginning of stem 2 will not always be paired
    
    int dstnc = 0;         //variable to track variable position differences during sequence comparision
    int dstnc_thrshld = 7; //minimum number of differences between barcodes
    int too_close = 0;     //flag that current barcode sequence is too closely related to a previous sequence
    int brcd_cnt = 0;      //number of barcodes selected for output
        
    int tot_seen = 0;      //number of barcodes tested
    int cnsctv_seen = 0;   //number of barcodes seen without finding an untested barcode
    
    FILE * out_fp = NULL;  //output file pointer
    
    //open output file
    if ((out_fp = fopen("barcodes.txt", "w")) == NULL) {
        printf("get_rndm_brcds: ERROR - could not generate barcodes output file. Aborting program...\n");
        abort();
    }
    
    fprintf(out_fp, "%d barcodes\n", brcds2mk); //first line of file is number of barcodes
    fprintf(out_fp, "linker=%s\n", brcd_lnkr);  //second line of file is the linker sequence that will precede barcode
    fprintf(out_fp, "bStruct=%s\n", brcd_sec);  //third line of file is the barcode secondary structure
    
    //select first barcode and print to file
    brcd_out[brcd_cnt] = &pf_brcds[rand() % passed_filter];         //set output pointer for 1st barcode
    printf("%d\t%s\n", brcd_cnt+1, brcd_out[brcd_cnt]->sq);         //print barcode to screen
    fprintf(out_fp, "%d\t%s\n",brcd_cnt+1, brcd_out[brcd_cnt]->sq); //print barcode to file
    brcd_cnt++;                                                     //increment barcode count
    
    //select barcodes until brcds2mk barcodes have been output
    //or until 100000 barcodes have been tested without outputting a barcode
    //TODO: add linear search that is performed after cnsctv_seen limit is reached?
    while ((brcd_cnt < brcds2mk) && cnsctv_seen < 100000) {
        
        too_close = 0;              //set too_close flag to zero
        r = rand() % passed_filter; //select barcode index
        
        if (!pf_brcds[r].mul) {     //if barcode was not tested yet
            pf_brcds[r].mul = 1;    //set flag that barcode was tested
            tot_seen++;             //increment number of seen barcodes
            cnsctv_seen = 0;        //set cnsctv_seen to 0
        } else {
            too_close = 1;          //barcode was already tested, set too_close to true to skip barcode
            cnsctv_seen++;          //increment cnsctv_seen
        }
        
        //check current barcode is sufficiently different than all output barcodes
        for (i = 0; brcd_out[i] != NULL && !too_close; i++) {
            
            //compare asterisked variable positions for current and output barcodes
            for (j = 0, dstnc = 0; pf_brcds[r].sq[j] && brcd_out[i]->sq[j] && !too_close; j++) {
                if (rPos[j] == '*' && (pf_brcds[r].sq[j] != brcd_out[i]->sq[j])) {
                    dstnc++;
                }
            }
            
            if (dstnc < dstnc_thrshld) { //too few differences betwen current and output barcode
                too_close = 1;           //set too_close to true
            }
        }
        
        if (!too_close) { //barcode is sufficiently different than all previously output barcodes
            brcd_out[brcd_cnt] = &pf_brcds[r];                            //set output pointer to current barcode
            printf("%d\t%s\n",brcd_cnt+1, brcd_out[brcd_cnt]->sq);        //print barcode to screen
            fprintf(out_fp, "%d\t%s\n",brcd_cnt, brcd_out[brcd_cnt]->sq); //print barcode to file
            brcd_cnt++;                                                   //increment barcode count
        }
    }
    
    printf("%d barcodes seen during random selection of %d barcodes\n", tot_seen, brcd_cnt);
    
    if (fclose(out_fp) == EOF) {
        printf("get_rndm_brcds: error - error occurred when closing barcode output file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

/* init_barcode_bank: initialize new barcode bank */
void init_barcode_bank (barcode_bank * bc_bnk, barcode_bank * prev)
{
    bc_bnk->fa = calloc(BC_BLOCK_SIZE, sizeof(*(bc_bnk->fa))); //allocate memory for barcode sequences
    bc_bnk->cnt = 0;     //initialize bank barcode count to zero
    bc_bnk->nxt = NULL;  //initialize pointer to next bank to NULL
    bc_bnk->prev = prev; //initialize pointer to previous bank to the 'prev' argument
}

/* add_barcode_bank: allocate a new barcode bank within the barcode bank linked list */
void add_barcode_bank (barcode_bank * bc_bnk)
{
    bc_bnk->nxt = calloc(1, sizeof(*(bc_bnk->nxt))); //allocate memory for the new barcode bank
    init_barcode_bank(bc_bnk->nxt, bc_bnk);          //initialize the new barcode bank
}

/* merge_barcode_bank: merge the barcode bank linked list into a fasta array */
void merge_barcode_bank(barcode_bank * bc_bnk, uint64_t cnt)
{
    extern fasta * vrnts;          //pointer to fasta structures used when generating variants
    extern uint64_t v_indx;        //index of current variant
    extern barcode_bank * bc_crnt; //pointer to current barcode bank
    
    barcode_bank * lcl_crnt_bnk = bc_bnk; //set local current barcode bank to bc_bnk argument
    
    int i = 0; //general purpose index
    
    //allocate memory to store barcode sequences
    if ((vrnts = calloc(cnt+1, sizeof(*vrnts))) == NULL) {
        printf("xpnd_brcd_tmplts: error - variant memory allocation failed. aborting...\n");
        abort();
    }
        
    int found_last_bank = 0; //flag that last barcode bank was found
    
    while (!found_last_bank) {
        for (i = 0; i < lcl_crnt_bnk->cnt; i++) { //for every variant in the current barcode bank
            
            //allocate barcode sequence storage in the vrnts array
            if ((vrnts[v_indx].sq = malloc((strlen(lcl_crnt_bnk->fa[i].sq)+1) * sizeof(*(vrnts[v_indx].sq)))) == NULL) {
                printf("merge_barcode_bank: variant memory allocation failed. aborting...\n");
                abort();
            }
            
            strcpy(vrnts[v_indx].sq, lcl_crnt_bnk->fa[i].sq); //store variant sequence in vrnts array
            v_indx++;                                         //increment variant index
            free(lcl_crnt_bnk->fa[i].sq);                     //free bc bank variant sequence memory
        }
        
        free(lcl_crnt_bnk->fa); //free bc_bank fasta memory
        
        if ((lcl_crnt_bnk = lcl_crnt_bnk->nxt) == NULL) { //if there is not another bank to merge
            found_last_bank = 1;                          //set flag that last bank was found
        }
    }
    
    lcl_crnt_bnk = bc_crnt; //set local current barcode bank to last barcode bank
    
    while (lcl_crnt_bnk->prev != NULL) {   //until the root barcode bank is reached
        lcl_crnt_bnk = lcl_crnt_bnk->prev; //set the local current bc bank to the previous bc bank
        free(lcl_crnt_bnk->nxt);           //free next barcode bank memory
    }
        
    return;
}

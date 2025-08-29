//
//  barcode_variants.c
//  
//
//  Created by Eric Strobel on 3/2/23.
//

//CRITICAL NOTE: check all indices to make sure uint64_t is used where needed

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
    
    int i = 0;                         //general purpose index
    uint64_t passed_filter = 0;        //number of barcodes that passed the parameter filter
    barcode_target * pf_brcds = NULL;  //pointer to target structures, used when picking random barcodes
    
    wt_source brcd_src[BARCODE_TEMPLATES] = {0}; //wt_source array for storing barcode template sequences
    basemap brcd_bmap[BARCODE_TEMPLATES] = {0};  //basemap array for storing barcode template sequences
    
    char out_dir[512] = {0};           //array to store output directory name
    sprintf(out_dir, "barcodes_out");  //construct output directory name
    mk_out_dir(out_dir);               //make output directory
    
    //NOTE: If multiple custom barcode template usage is enabled in the future, will need to
    //include check that all barcode length/randomization schemes match, or edit code to accommodate
    //usage of distinct barcode strategies simultaneously.
    barcode_template brcd_tmplt[BARCODE_TEMPLATES] = {"barcode_template",
                                                      "NNNNNNNNNNNNNNNN",
                                                      "................",
                                                      "****************"};
    
    if (strlen(brcd_tmplt->sq) > MAX_BARCODE_LEN) {
        printf("make_barcodes: error - barcode template exceeds maximum length (%d nucleotides). aborting...\n", MAX_BARCODE_LEN);
        abort();
    }
    
    //generate input details file
    sprintf(prcs_out_nm, "./%s/barcodes_processing.txt", out_dir);
    if ((prcs_ofp = fopen(prcs_out_nm, "w")) == NULL) {
        printf("mk_brcds: ERROR - could not open processing messages file. Aborting program...\n");
        abort();
    }
    
    init_brcd_tmplts(&brcd_src[0], &brcd_bmap[0], &brcd_tmplt[0]); //initialize variant templates for barcode generation
    printf("generating and filtering possible barcodes\n");
    passed_filter = xpnd_brcd_tmplts(&brcd_bmap[0]); //expand barcode variant templates
    printf("passed filter: %llu\n", (long long unsigned int)passed_filter);
    
    //allocate memory to store passed filter barcodes in targets structures
    if ((pf_brcds = calloc(passed_filter, sizeof(*pf_brcds))) == NULL) {
        printf("mk_brcds: error - passed filter barcodes memory allocation failed. aborting...\n");
        abort();
    }
    
    //merge barcode bank by storing barcodes in barcode_targets structure array
    printf("storing passed filter barcodes in targets structures\n");
    merge_barcode_bank(pf_brcds, &bc_root, passed_filter);
    
    //allocate target structure pointers for output barcodes
    barcode_target * brcd_out = NULL; //pointers to output barcode structures
    if ((brcd_out = calloc(brcds2mk, sizeof(*brcd_out))) == NULL) {
        printf("mk_brcds: error - selected barcode memory allocation failed. aborting...\n");
        abort();
    }
    
    //randomly select barcodes to output
    uint8_t tbl[256] = {0};     //array for storing hamming distance lookup table
    mk_8bit_dstnc_tbl(&tbl[0]); //generate hamming distance lookup table
    get_rndm_brcds(brcd_out, pf_brcds, passed_filter, brcds2mk, &brcd_tmplt[0], &tbl[0]); //select barcodes to output
    
    //close output file
    if (fclose(prcs_ofp) == EOF) {
        printf("main: error - error occurred when closing processing messages output file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

/* init_brcd_tmplts: initialize variant templates for barcode sequences */
int init_brcd_tmplts(wt_source * brcd_src, basemap * brcd_bmap, barcode_template brcd_tmplt[BARCODE_TEMPLATES])
{
    printf("initializing barcode templates\n");
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int vars = 0; //total number of variants encoded by barcode templates
    
    //initialize a variant template for each barcode sequence
    for (i = 0; i < BARCODE_TEMPLATES; i++) {
        
        //allocate memory for sequence name in wt_source struct
        if ((brcd_src[i].nm = malloc((strlen(brcd_tmplt[i].nm)+1) * sizeof(*(brcd_src[i].nm)))) == NULL) {
            printf("init_brcd_tmplts: error - memory allocation for sequence name failed. aborting...\n");
            abort();
        }
        
        //allocate memory for source sequence in wt_source struct
        if ((brcd_src[i].sq = malloc((strlen(brcd_tmplt[i].sq)+1) * sizeof(*(brcd_src[i].sq)))) == NULL) {
            printf("init_brcd_tmplts: error - memory allocation for source sequence failed. aborting...\n");
            abort();
        }
        
        //allocate memory for source sequence positions in wt_source struct
        if ((brcd_src[i].pos = malloc((strlen(brcd_tmplt[i].sq)+1) * sizeof(*(brcd_src[i].pos)))) == NULL) {
            printf("init_brcd_tmplts: error - memory allocation for source sequence failed. aborting...\n");
            abort();
        }
        
        strcpy(brcd_src[i].nm, brcd_tmplt[i].nm); //set barcode name in wt_source
        strcpy(brcd_src[i].sq, brcd_tmplt[i].sq); //set barcode sequence in wt_source
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
        init_vtmp_mem(&brcd_bmap[i], strlen(brcd_tmplt[i].nm), strlen(brcd_tmplt[i].sq));
        strcpy(brcd_bmap[i].nm, brcd_tmplt[i].nm);    //set barcode template name
        strcpy(brcd_bmap[i].rS, brcd_tmplt[i].sq);    //set barcode template sequence
        strcpy(brcd_bmap[i].rP[0], brcd_tmplt[i].pr); //set barcode template pair constraints
        brcd_bmap[i].rP_cnt = 1;                      //set pair constraint count to 1
        
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
    
    char var_outpt[MAXLEN+1] = {0};       //array to store sequences during barcode expansion
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
    int max_strong_heteropol = 0;
    int crnt_strong_heteropol = 0;
    char prev_base = '\0'; //identity of previous base, used for homopolymer tracking
        
    for (j = 0; sq[j]; j++) {
        
        //track base identity
        switch (sq[j]) {
            case 'A':
                at_cnt++;
                if (crnt_strong_heteropol > max_strong_heteropol) {
                    max_strong_heteropol = crnt_strong_heteropol;
                }
                crnt_strong_heteropol = 0;
                break;
                
            case 'T':
                at_cnt++;
                if (crnt_strong_heteropol > max_strong_heteropol) {
                    max_strong_heteropol = crnt_strong_heteropol;
                }
                crnt_strong_heteropol = 0;
                break;
                
            case 'G':
                gc_cnt++;
                crnt_strong_heteropol++;
                break;
                
            case 'C':
                gc_cnt++;
                crnt_strong_heteropol++;
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
        max_strong_heteropol > STRONG_HETEROPOL_LIMIT) {
        return 1; //return fail
    } else {
        return 0; //return pass
    }
}

/* get_rndom_brcds: select random passed filter barcodes to output*/
int get_rndm_brcds(barcode_target * brcd_out, barcode_target * pf_brcds, uint64_t passed_filter, int brcds2mk, barcode_template brcd_tmplt[BARCODE_TEMPLATES], uint8_t * tbl)
{
    extern int debug;  //flag for running debug mode
    
    srand(time(NULL)); //seed random number generation
    
    int i = 0;       //general purpose index
    int j = 0;       //general purpose index
    uint64_t r = 0;  //randomized index for selecting barcodes
        
    int brcd_cnt = 0;       //number of barcodes selected for output
    
    uint64_t diff_bits = 0; //difference (XOR) between two binary sequences
    int b_dstnc = 0;        //binary-encoded sequence hamming distance
    int min_dstnc = 6;      //minimum number of differences between barcodes
    int too_close = 0;      //flag that current barcode sequence is too closely related to a previous sequence
    barcode_target * fRef = NULL; //pointer to barcode target for performing failed match trimming
    
    uint64_t crrnt_cnt = passed_filter;  //number of variants in current pool
    uint64_t crrnt_tested = 0;           //number of variants in current pool that have been tested
    uint64_t tot_seen = 0;               //number of barcodes tested
    uint64_t prox_brcds_flagged = 0;     //number of barcodes flagged by a proximal search
    uint64_t sMatch_trim = 0;            //total number of proximal barcodes trimmed after a successful search
    uint64_t fMatch_trim = 0;            //total number of proximal barcodes trimmed after a failed search
    
    int sntnl_cnt = 0;              //number of sentinels
    uint32_t min_sntnl_dstnc = 0;   //minimum observed distance from sentinels
    uint64_t xcldd_by_sntnls = 0;   //number of variants that were excluded by sentinels
    barcode_target * sntnls = NULL; //storage for sentinel sequences
    if ((sntnls = calloc(brcds2mk, sizeof(*sntnls))) == NULL) {
        printf("get_rndm_brcds: error - sentinel target pointer memory allocation failed. aborting...\n");
        abort();
    }

    //open output file
    FILE * out_fp = NULL;  //output file pointer
    if ((out_fp = fopen("barcodes.txt", "w")) == NULL) {
        printf("get_rndm_brcds: ERROR - could not generate barcodes output file. Aborting program...\n");
        abort();
    }
    
    //select first barcode
    r = ((uint64_t)rand()) % passed_filter;  //select random barcode index
    pf_brcds[r].mul = 1;                     //set flag that barcode was "tested"
    store_brcd(&brcd_out[0], &brcd_cnt, &pf_brcds[r], sntnls, &sntnl_cnt, NULL, min_dstnc+1, min_dstnc, brcds2mk);
    
    tot_seen++;     //increment total number of barcodes seen
    crrnt_tested++; //increment barcode tested counter
    
    //flag proximal, too closely related barcodes for trimming
    prox_brcds_flagged = flag_proximal_vrnts(&pf_brcds[r], pf_brcds, crrnt_cnt, r, min_dstnc, &brcd_tmplt[0], tbl);
    tot_seen += prox_brcds_flagged;     //increment total barcodes seen by the number of flagged barcodes
    crrnt_tested += prox_brcds_flagged; //increment current tested barcodes by number of flagged barcodes
    sMatch_trim += prox_brcds_flagged;  //track number of proximal barcodes flagged after a successful match
    
    //select barcodes until brcds2mk barcodes have been output
    //or until 100000000 barcodes have been tested without outputting a barcode
    while (crrnt_cnt && (brcd_cnt < brcds2mk) && crrnt_tested < 100000000) {
        
        too_close = 0;                       //set too_close flag to zero
        r = ((uint64_t)rand()) % crrnt_cnt;  //select barcode index
                
        if (!pf_brcds[r].mul) {  //if barcode was not tested yet
            
            pf_brcds[r].mul = 1; //set flag that barcode was tested
            tot_seen++;          //increment number of seen barcodes
            crrnt_tested++;      //increment count of barcodes tested
            fRef = NULL;         //set failed match reference to NULL
            
            //check selected barcode against sentinels first. if too close, point failed match reference to too close sentinel
            //if not too close to any sentinel, check against all output barcodes
            if ((fRef = test_sentinels(&pf_brcds[r].bsq, sntnls, sntnl_cnt, min_dstnc, tbl, &min_sntnl_dstnc)) == NULL) {
                
                //check current barcode is sufficiently different than all output barcodes
                for (i = 0; brcd_out[i].bsq.sq != NULL && !too_close; i++) {
                    
                    b_dstnc = calc_binseq_dstnc(&pf_brcds[r].bsq, &brcd_out[i].bsq, tbl); //calculate distance
                    if (b_dstnc < min_dstnc) { //too few differences betwen current and output barcode
                        too_close = 1;         //set too_close to true
                        fRef = &brcd_out[i];   //point failed match reference to too close barcode
                    }
                    
                    //if running debug mode, also assess distance using character string and compare to b_dstnc
                    if (debug) {test_c_dstnc(pf_brcds[r].sq, brcd_out[i].sq, &brcd_tmplt[0], b_dstnc, "Mtest");}
                }
            } else {
                too_close = 1;     //set too_close to true
                xcldd_by_sntnls++; //track number of barcodes excluded by sentinels
            }
        } else {
            too_close = -1; //barcode was already tested, set too_close to -1 to skip barcode storage/proximal barcode flagging
        }
        
        //if the current barcode was not previously seen, perform storage and/or proximal barcode trimming
        //depending on whether the barcode was sufficiently different from previously selected barcodes
        
        prox_brcds_flagged = 0; //set proximal barcodes flagged to zero
        
        if (!too_close) {       //barcode is sufficiently different than all previously output barcodes
                                //store barcode and flag barcodes that are proximal to the newly selected barcode
            
            store_brcd(&brcd_out[0], &brcd_cnt, &pf_brcds[r], sntnls, &sntnl_cnt, fRef, min_sntnl_dstnc, min_dstnc, brcds2mk);
            prox_brcds_flagged = flag_proximal_vrnts(&pf_brcds[r], pf_brcds, crrnt_cnt, r, min_dstnc, &brcd_tmplt[0], tbl);
            sMatch_trim += prox_brcds_flagged;  //track number of proximal barcodes flagged after a successful match
            
        } else if (too_close == 1) {  //barcode is too close to previously selected barcode
                                      //flag barcodes that are proximal to the too close previously selected barcode
            
            prox_brcds_flagged = flag_proximal_vrnts(fRef, pf_brcds, crrnt_cnt, r, min_dstnc, &brcd_tmplt[0], tbl);
            fMatch_trim += prox_brcds_flagged;  //track number of proximal barcodes flagged after a failed match
        }
        
        tot_seen += prox_brcds_flagged;     //increment total barcodes seen by the number of flagged barcodes
        crrnt_tested += prox_brcds_flagged; //increment current tested barcodes by number of flagged barcodes
        
        
        //when 5% of the current number of candidate barcodes is tested, trim tested barcodes
        if ((((double)crrnt_tested) / ((double)crrnt_cnt)) > 0.05) {
            crrnt_cnt = trim_barcodes(pf_brcds, crrnt_cnt, crrnt_cnt - crrnt_tested);
            crrnt_tested = 0;
        }
    }
    
    if (debug) { //in debug mode, print barcode exclusion metrics
        printf("%llu excluded by successful match proximal trim\n", (long long unsigned int)sMatch_trim);
        printf("%llu excluded by failed match proximal trim\n", (long long unsigned int)fMatch_trim);
        printf("%d sentinels\n%llu excluded by sentinels\n", sntnl_cnt, (long long unsigned int)xcldd_by_sntnls);
    }
    
    //print output
    fprintf(out_fp, "%d barcodes\n", brcd_cnt);        //first line of file is number of barcodes
    fprintf(out_fp, "bStruct=%s\n", brcd_tmplt[0].pr); //second line of file is the barcode secondary structure
    //TODO: if barcode design flexibility is added, need to require barcode structure to be idential for all BC templates
    
    for (i = 0; i < brcd_cnt; i++) { //print barcodes to file
        fprintf(out_fp, "%d\t%s\n", i+1, brcd_out[i].sq);
    }
    
    if (fclose(out_fp) == EOF) { //close file
        printf("get_rndm_brcds: error - error occurred when closing barcode output file. Aborting program...\n");
        abort();
    }
    
    printf("%llu barcodes seen during random selection of %d barcodes\n", (long long unsigned int)tot_seen, brcd_cnt);
    
    return 1;
}

/* store_brcd: store selected barcode in output barcodes array */
void store_brcd(barcode_target * brcd_out, int * brcd_cnt, barcode_target * bc2set, barcode_target * sntnls, int * sntnl_cnt, barcode_target * fRef, uint32_t min_sntnl_dstnc, int min_dstnc, int brcds2mk)
{
    copy_binary_seq(&brcd_out[*brcd_cnt].bsq, &bc2set->bsq); //copy binary seq struct to output barcodes array
    
    //store character-encoded barcode sequence in output barcodes array
    if ((brcd_out[*brcd_cnt].sq = malloc((int)(brcd_out[*brcd_cnt].bsq.ln+1) * sizeof(*(brcd_out[*brcd_cnt].sq)))) == NULL) {
        printf("get_rndm_brcds: error - failed to allocate memory for barcode sequence\n");
        abort();
    }
    bin2seq(brcd_out[*brcd_cnt].sq, &brcd_out[*brcd_cnt].bsq, brcd_out[*brcd_cnt].bsq.ln+1);
       
    //if new barcode is sufficiently different from other sentinels, add new sentinel
    if (fRef == NULL && min_sntnl_dstnc > min_dstnc && *sntnl_cnt < brcds2mk) {
        
        copy_binary_seq(&sntnls[*sntnl_cnt].bsq, &brcd_out[*brcd_cnt].bsq); //copy binary seq struct
        
        //store character-encoded barcode sequence in sentinels array
        if ((sntnls[*sntnl_cnt].sq = malloc((int)(strlen(brcd_out[*brcd_cnt].sq)+1) * sizeof(*(sntnls[*sntnl_cnt].sq)))) == NULL) {
            printf("get_rndm_brcds: error - failed to allocate memory for barcode sequence\n");
            abort();
        }
        strcpy(sntnls[(*sntnl_cnt)++].sq, brcd_out[*brcd_cnt].sq); //copy character-encoded sequence and increment sntnl_cnt
    }
    
    printf("%d\t%s\n", *brcd_cnt+1, brcd_out[*brcd_cnt].sq); //print barcode to screen
    (*brcd_cnt)++;                                           //increment barcode count
    
    return;
}

/* calc_binseq_dstnc: calculate hamming distance between two binary-encoded sequences using a lookup table */
int calc_binseq_dstnc(binary_seq * bsq1, binary_seq * bsq2, uint8_t * tbl)
{
    uint64_t dstnc = 0;                         //hamming distance
    uint64_t diff_bits = *bsq1->sq ^ *bsq2->sq; //XOR to identifiy different bits
    
    //sum distance of each byte to calculate total distance
    //Note: this calculation is set up for 16 nt long barcodes only
    dstnc = tbl[(uint8_t)((diff_bits & 0xFF00000000000000) >> 56)] +
            tbl[(uint8_t)((diff_bits & 0x00FF000000000000) >> 48)] +
            tbl[(uint8_t)((diff_bits & 0x0000FF0000000000) >> 40)] +
            tbl[(uint8_t)((diff_bits & 0x000000FF00000000) >> 32)];
        
    return dstnc; //return hamming distance
}

/* test_sentinels: test barcode against sentinel barcodes */
barcode_target * test_sentinels(binary_seq * bsq, barcode_target * sntnls, int sntnl_cnt, int min_dstnc, uint8_t * tbl, uint32_t * min_sntnl_dstnc)
{
    int i = 0;                     //general purpose index
    int crrnt_dstnc = 0;           //distance between candidate barcode and current sentinel
    int too_close = 0;             //flag that the barcode is too close to sentinel
    barcode_target * fRef = NULL;  //pointer to sentinel to which candidate barcode was too close
    
    *min_sntnl_dstnc = 0xFFFFFFFF; //set minimum sentinel distance to max unsigned 32-bit integer value
    
    //test whether candidate barcode is too close to sentinel sentinel barcodes
    for (i = 0, too_close = 0; i < sntnl_cnt && !too_close; i++) {
        
        
        if ((crrnt_dstnc = calc_binseq_dstnc(bsq, &sntnls[i].bsq, tbl)) < min_dstnc) { //if distance is less than min distance
            too_close = 1;     //turn too_close on
            fRef = &sntnls[i]; //set failed match reference pointer to current sentinel
        }
        
        if (crrnt_dstnc < *min_sntnl_dstnc) { //if current distance is less than min sentinel distance
            *min_sntnl_dstnc = crrnt_dstnc;   //set min sentinel distance to current distance
        }
    }
    
    return fRef; //return pointer to too close sentinel barcode (NULL if candidate was not too close to any sentinel)
}

/* flag_proximal_vrnts: perform upstream and downstream searches for barcodes
   that are too close to a reference barcode and flag hits for trimming */
uint64_t flag_proximal_vrnts(barcode_target * ref, barcode_target * pf_brcds, int crrnt_cnt, int r, int min_dstnc, barcode_template brcd_tmplt[BARCODE_TEMPLATES], uint8_t * tbl)
{
    extern int debug; //flag for running debug mode
    
    uint64_t i = 0;   //general purpose index
    int j = 0;        //general purpose index
    
    int b_dstnc = 0;  //hamming distance between binary-encoded sequences
      
    uint64_t diff_bits = 0;      //bitwise difference (XOR) between current barcode and reference barcode
    uint64_t dnstrm_flagged = 0; //number of downstream barcodes trimmed
    uint64_t upstrm_flagged = 0; //number of upstream barcodes trimmed
    
    //TODO: add check that barcodes are same length here and elsewhere? for potential future uses? or add something to disallow barcode templates with different lengths
    
    //flag downstream variants
    for (i = 1, b_dstnc = 0; (r+i) < crrnt_cnt && b_dstnc < min_dstnc; i++) {
        
        b_dstnc = calc_binseq_dstnc(&ref->bsq, &pf_brcds[r+i].bsq, tbl); //calculate hamming distance
        if (b_dstnc < min_dstnc && !pf_brcds[r+i].mul) { //if distance is too low and barcode was not yet flagged
            pf_brcds[r+i].mul = 1;                       //flag barcode for trimming
            dnstrm_flagged++;                            //increment downstream flagged counter
        }
        
        //if running debug mode, also assess distance using character string and compare to b_dstnc
        if (debug) {test_c_dstnc(ref->sq, pf_brcds[r+i].sq, &brcd_tmplt[0], b_dstnc, "Dtrim");}
    }
    
    //flag upstream variants
    for (i = 1, b_dstnc = 0; i <= r && b_dstnc < min_dstnc; i++) {
        
        b_dstnc = calc_binseq_dstnc(&ref->bsq, &pf_brcds[r-i].bsq, tbl); //calculate hamming distance
        if (b_dstnc < min_dstnc && !pf_brcds[r-i].mul) { //if distance is too low and barcode was not yet flagged
            pf_brcds[r-i].mul = 1;                       //flag barcode for trimming
            upstrm_flagged++;                            //increment upstream flagged counter
        }
        
        //if running debug mode, also assess distance using character string and compare to b_dstnc
        if (debug) {test_c_dstnc(ref->sq, pf_brcds[r-i].sq, &brcd_tmplt[0], b_dstnc, "Utrim");}
    }
    
    return  dnstrm_flagged + upstrm_flagged; //return total barcodes flagged
}

/* trim_barcodes: trim tested barcodes from barcodes array */
uint64_t trim_barcodes(barcode_target * pf_brcds, uint64_t crrnt_brcd_cnt, uint64_t brcds2keep)
{
    printf("\ntrimming %llu of %llu barcodes...\n", (long long unsigned int)(crrnt_brcd_cnt - brcds2keep), (long long unsigned int)crrnt_brcd_cnt);
 
    extern int debug; //flag for running debug mode
    
    uint64_t i = 0;   //general purpose index
    uint64_t j = 0;   //general purpose index
    
    uint64_t trimmed_count = 0; //number of trimmed barcodes
    
    uint64_t * freed_indices = NULL;              //array of indices that have been freed
    uint64_t nxt_free_lkup = 0xFFFFFFFFFFFFFFFF;  //freed indices array location of next free index
    
    uint64_t kept_bcs = 0; //number of kept barcodes
    
    //allocate memory for freed indices
    if ((freed_indices = calloc(crrnt_brcd_cnt, sizeof(*freed_indices))) == NULL) {
        printf("trim_barcodes: error - failed to allocate memory tracking freed indices\n");
        abort();
    }
    
    //initialize freed indices to max value of an unsigned 64-bit int
    for (i = 0; i < crrnt_brcd_cnt; i++) {
        freed_indices[i] = 0xFFFFFFFFFFFFFFFF;
    }
    
    for (i = 0, j = 0; i < crrnt_brcd_cnt; i++) { //until all current barcodes have been assessed
        
        if (pf_brcds[i].mul) {         //if current barcode was already tested
            trimmed_count++;           //increment number of trimmed barcodes
            free(pf_brcds[i].bsq.sq);  //free trimmed barcode binary-encoded sequence
            free(pf_brcds[i].sq);      //free trimmed barcode character-encoded sequence
            pf_brcds[i].bsq.sq = NULL; //set binary-encoded sequence pointer to NULL
            pf_brcds[i].sq = NULL;     //set character-encoded sequence pointer to NULL
            pf_brcds[i].mul = 0;       //set already-tested flag to false
            freed_indices[j] = i;      //add freed index to freed indices array
            if (nxt_free_lkup == 0xFFFFFFFFFFFFFFFF) { //if nxt_free_lkup has not yet been set,
                nxt_free_lkup = j;                     //set to j
            }
            j++; //increment j
            
        } else if (nxt_free_lkup != 0xFFFFFFFFFFFFFFFF) { //if current bc has not yet been tested and there is a free indx
            
            //shift barcode
            if (pf_brcds[freed_indices[nxt_free_lkup]].mul == 0 &&    //check that next free index is free
                pf_brcds[freed_indices[nxt_free_lkup]].sq  == NULL) {
                
                kept_bcs++; //increment kept barcodes counter,
                copy_binary_seq(&pf_brcds[freed_indices[nxt_free_lkup]].bsq, &pf_brcds[i].bsq); //copy binary seq struct
                free(pf_brcds[i].bsq.sq);   //free copied binary_seq sequence memory
                pf_brcds[i].bsq.sq = NULL;  //set copied binary_seq pointer to null
                
                if (debug) { //in debug mode, also copy character-encoded sequence
                    //allocate memory for character-encoded sequence
                    if ((pf_brcds[freed_indices[nxt_free_lkup]].sq = malloc((strlen(pf_brcds[i].sq)+1) * sizeof(*(pf_brcds[freed_indices[nxt_free_lkup]].sq)))) == NULL) {
                        printf("trim_barcodes: error - failed to allocate memory for barcode sequence\n");
                        abort();
                    }
                    strcpy(pf_brcds[freed_indices[nxt_free_lkup]].sq, pf_brcds[i].sq); //copy character-encoded sequence
                    free(pf_brcds[i].sq);                                              //free previous barcode target memory
                    pf_brcds[i].sq = NULL;                                             //and set pointer to NULL
                }
                
                nxt_free_lkup++;        //increment next free lookup index
                freed_indices[j++] = i; //add index i to freed indices and increment freed_indices counter
                
            } else { //next free index is not free, throw error and abort
                printf("trim_barcodes: error - target marked as freed is not free. aborting...\n");
                abort();
            }
        }  else { //current bc not tested and no free indx
            //printf("barcode was not tested and there is no free index, proceeding to next barcode\n");
            kept_bcs++;
        }
    }
      
    free(freed_indices); //free freed_indices array
    
    printf("trimmed  %llu barcodes, %llu barcodes remain\n\n", (long long unsigned int)trimmed_count, (long long unsigned int)kept_bcs);
    
    return kept_bcs;
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
void merge_barcode_bank(barcode_target * pf_brcds, barcode_bank * bc_bnk, uint64_t cnt)
{
    extern int debug;              //flag for running debug mode
    extern barcode_bank * bc_crnt; //pointer to current barcode bank
    
    barcode_bank * lcl_crnt_bnk = bc_bnk; //set local current barcode bank to bc_bnk argument
    
    int i = 0;      //general purpose index
    int b_indx = 0; //barcode index
    
    int found_last_bank = 0;   //flag that last barcode bank was found
    
    char test_array[33] = {0};
    
    while (!found_last_bank) {                    //until the last bank is found
        for (i = 0; i < lcl_crnt_bnk->cnt; i++) { //for every variant in the current barcode bank
            
            //allocate memory for 2-bit sequence in target structure
            if ((pf_brcds[b_indx].bsq.sq = calloc(BSQ_ARRAY_MAX, sizeof(*(pf_brcds[b_indx].bsq.sq)))) == NULL) {
                printf("merge_barcode_bank: error - failed to allocate memory for binary-encoded barcode sequence\n");
                abort();
            }
            seq2bin_long(lcl_crnt_bnk->fa[i].sq, &pf_brcds[b_indx].bsq, BSQ_ARRAY_MAX); //store 2-bit sequence
            
            //bin2seq conversion test
            if (debug) {
                bin2seq(test_array, &pf_brcds[b_indx].bsq, 17);
                if (strcmp(lcl_crnt_bnk->fa[i].sq, test_array)) {
                    printf("conversion error\n");
                }
            }
            
            //if running debug mode, allocate memory for storing barcode sequence in target structure
            if (debug) {
                if ((pf_brcds[b_indx].sq = malloc((strlen(lcl_crnt_bnk->fa[i].sq)+1) * sizeof(*(pf_brcds[b_indx].sq)))) == NULL) {
                    printf("merge_barcode_bank: error - failed to allocate memory for barcode sequence\n");
                    abort();
                }
                strcpy(pf_brcds[b_indx].sq, lcl_crnt_bnk->fa[i].sq); //copy barcode sequence to target structure
            }
            
            b_indx++;                     //increment barcode index
            free(lcl_crnt_bnk->fa[i].sq); //free memory for current barcode in vrnts array to keep memory use staticS
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
        
    //check that number of sequences stored in pf_brcds matches the number of sequences that passed filter
    if (b_indx != cnt) {
        printf("store_pf_barcodes: error - variants stored in pf_brcds array (%llu) does not equal the number of passed filter variants(%d). aborting...\n", (long long unsigned int)cnt, b_indx);
        abort();
    }
    
    return;
}

/* mk_8bit_dstnc_tbl: make table for looking up hamming distance between two 8-bit binary encoded sequences */
void mk_8bit_dstnc_tbl(uint8_t * tbl)
{
    int i = 0;     //general purpose index, necesary because b overflows after last loop iteration
    uint8_t b = 0; //8-bit int for calculating hamming distance
    
    uint8_t bits_7_6 = 0xC0; //1100 0000 mask for bits 7 and 6
    uint8_t bits_5_4 = 0x30; //0011 0000 mask for bits 5 and 4
    uint8_t bits_3_2 = 0x0C; //0000 1100 mask for bits 3 and 2
    uint8_t bits_1_0 = 0x03; //0000 0011 mask for bits 1 and 0
    
    //for each integer of an 8-bit number, sum the number
    //of bit pairs in which one or more bits are on
    for (i = 0, b = 0; i <= 255; i++, b++) {
        tbl[i] += (b & bits_7_6) ? 1 : 0;
        tbl[i] += (b & bits_5_4) ? 1 : 0;
        tbl[i] += (b & bits_3_2) ? 1 : 0;
        tbl[i] += (b & bits_1_0) ? 1 : 0;
    }
    
    return;
}

/* calc_string_dstnc: calculate hamming distance between two character-encoded sequences */
int calc_string_dstnc(char * str1, char * str2, barcode_template brcd_tmplt[BARCODE_TEMPLATES])
{
    int i = 0;     //general purpose index
    int dstnc = 0; //hamming distance
    
    if (strlen(str1) != strlen(str2)) { //check that string lengths match
        printf("calc_string_dstnc: input strings are not the same length. aborting...\n");
        abort();
    }
    
    for (i = 0, dstnc = 0; str1[i] && str2[i]; i++) {              //for each character of strings 1 and 2
        //NOTE: line below only accommodates one barcode template
        if (brcd_tmplt[0].rnd[i] == '*' && (str1[i] != str2[i])) { //if i is a randomized position and the chars don't match
            dstnc++;                                               //increment hamming distance
        }
    }
    
    return dstnc; //return hamming distance
}

/* test_c_dstnc: check that character-encoded seq difference matches binary-encoded seq difference */
void test_c_dstnc(char * sq1, char * sq2, barcode_template brcd_tmplt[BARCODE_TEMPLATES], int b_dstnc, char * test_loc)
{
    int c_dstnc = 0; //character-encoded sequence hamming distance
    
    c_dstnc = calc_string_dstnc(sq1, sq2, &brcd_tmplt[0]); //calculate character-encoded sequence hamming distance
    
    if (b_dstnc != c_dstnc) { //if c_dstnc does not match b_dstnc, throw error
        printf("\n%s:  b_dstnc does not match c_dstnc: b=%d c=%d\n", test_loc, b_dstnc, c_dstnc);
        printf("%s\n%s\n", sq1, sq2);
    }
    
    return;
}

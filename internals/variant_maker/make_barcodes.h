//
//  barcode_variants.h
//  
//
//  Created by Eric Strobel on 3/2/23.
//

#ifndef make_barcodes_h
#define make_barcodes_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../seq_utils/ispair.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/seq2bin_long.h"
#include "../seq_utils/basemap.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"
#include "./expand_variant_template.h"

#define BC_BLOCK_SIZE 65536      //barcode memory allocation block size

#define BARCODE_TEMPLATES 1      //number of variant templates used for barcode generation
#define BSQ_ARRAY_MAX 1          //bacode sequence array maximum - curently use one barcode template,
                                 //but code is mostly set up to accommodate more

#define MAX_BARCODE_LEN 16

//barcode parameter definitions
#define HOMOPOL_LIMIT 2          //maximum length of a homopolymeric stretch
#define STRONG_HETEROPOL_LIMIT 3 //maximum length of strong base pair stretch
#define BRCD_GC_MIN 6            //minimum GC content ( 6 for 16N)
#define BRCD_GC_MAX 10           //maximum GC content (10 for 16N)

/* barcode_bank: storage for passed filter barcodes during barcode generation. bank is a doubly linked list */
typedef struct barcode_bank {
    fasta *fa;                 //fasta struct pointer for barcode sequence storage allocation
    struct barcode_bank *nxt;  //pointer to next barcode bank
    struct barcode_bank *prev; //pointer to previous barcode bank
    int cnt;                   //number of barcodes in bank
} barcode_bank;

/* barcode_template: storage for barcode generation parameters */
typedef struct barcode_template {
    char nm[128];
    char sq[128];
    char pr[128];
    char rnd[128];
} barcode_template;

/* barcode_target: storage for compressed barcode and, in debug mode, uncompressed barcode */
typedef struct barcode_target {
    char * sq;       //charcter-encoded target sequence
    binary_seq bsq;  //2-bit-encoded target sequence
    uint8_t mul;     //flag that sequence is identical to a prior target
} barcode_target;

/* mk_brcds: coordinates barcode generation */
int mk_brcds(int brcds2mk);

/* init_brcd_tmplts: initialize variant templates for barcode sequences */
int init_brcd_tmplts(wt_source * brcd_src, basemap * brcd_bmap, barcode_template brcd_tmplt[BARCODE_TEMPLATES]);

/* xpnd brcd_tmplts: expand barcode templates to generate all possible barcode variants */
uint64_t xpnd_brcd_tmplts(basemap * brcd_bmap);

/* filter_barcode: filter barcodes that do not meet the parameters specified in make_barcodes.h */
int filter_barcode(char * sq);

/* get_rndom_brcds: select random passed filter barcodes to output*/
int get_rndm_brcds(barcode_target * brcd_out, barcode_target * pf_brcds, uint64_t passed_filter, int brcds2mk, barcode_template brcd_tmplt[BARCODE_TEMPLATES], uint8_t * tbl);

/* store_brcd: store selected barcode in output barcodes array */
void store_brcd(barcode_target * brcd_out, int * brcd_cnt, barcode_target * bc2set, barcode_target * sntnls, int * sntnl_cnt, barcode_target * fRef, uint32_t min_sntnl_dstnc, int min_dstnc, int brcds2mk);

/* calc_binseq_dstnc: calculate hamming distance between two binary-encoded sequences using a lookup table */
int calc_binseq_dstnc(binary_seq * bsq1, binary_seq * bsq2, uint8_t * tbl);

/* test_sentinels: test barcode against sentinel barcodes */
barcode_target * test_sentinels(binary_seq * bsq, barcode_target * sntnls, int sntnl_cnt, int min_dstnc, uint8_t * tbl, uint32_t * min_sntnl_dstnc);

/* flag_proximal_vrnts: perform upstream and downstream searches for barcodes
   that are too close to a reference barcode and flag hits for trimming */
uint64_t flag_proximal_vrnts(barcode_target * ref, barcode_target * pf_brcds, int crrnt_cnt, int r, int min_dstnc, barcode_template brcd_tmplt[BARCODE_TEMPLATES], uint8_t * tbl);

/* trim_barcodes: trim tested barcodes from barcodes array */
uint64_t trim_barcodes(barcode_target * pf_brcds, uint64_t crrnt_brcd_cnt, uint64_t brcds2keep);

/* init_barcode_bank: initialize new barcode bank */
void init_barcode_bank(barcode_bank * bc_bnk, barcode_bank * prev);

/* add_barcode_bank: allocate a new barcode bank within the barcode bank linked list */
void add_barcode_bank(barcode_bank * bc_bnk);

/* merge_barcode_bank: merge the barcode bank linked list into a fasta array */
void merge_barcode_bank(barcode_target * pf_brcds, barcode_bank * bc_bnk, uint64_t cnt);

/* mk_8bit_dstnc_tbl: make table for looking up hamming distance between two 8-bit binary encoded sequences */
void mk_8bit_dstnc_tbl(uint8_t * tbl);

/* calc_string_dstnc: calculate hamming distance between two character-encoded sequences */
int calc_string_dstnc(char * str1, char * str2, barcode_template brcd_tmplt[BARCODE_TEMPLATES]);

/* test_c_dstnc: check that character-encoded seq difference matches binary-encoded seq difference */
void test_c_dstnc(char * sq1, char * sq2, barcode_template brcd_tmplt[BARCODE_TEMPLATES], int b_dstnc, char * test_loc);

#endif /* barcode_variants_h */

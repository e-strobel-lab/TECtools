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

#define BC_BLOCK_SIZE 65536

#define BARCODE_TEMPLATES 4   //number of variant templates used for barcode generation
#define STEM1BOT_N_STRT 0     //start of stem 1 bottom N stretch
#define STEM1BOT_N_END  3     //  end of stem 1 bottom N stretch
#define STEM1TOP_N_STRT 4     //start of stem 1 top N stretch
#define STEM1TOP_N_END  7     //  end of stem 1 top N stretch
#define STEM2BOT_STRT   27
#define STEM2BOT_END    28
#define STEM2TOP_STRT   36
#define STEM2TOP_END    37

//barcode parameter definitions
#define HOMOPOL_LIMIT 2       //maximum length of a homopolymeric stretch
#define S1BOT_GC_MIN 1        //maximum number of AT pairs allowed at the base of stem 1
#define S1BOT_GC_MAX 3        //maximum number of GC pairs allowed at the base of stem 1
#define S1TOP_GC_MIN 1        //maximum number of AT pairs allowed at the top of stem 1
#define S1TOP_GC_MAX 3        //maximum number of GC pairs allowed at the top of stem 1
#define S1ALL_GC_MIN 3        //minimum number of GC pairs allowed in stem 1
#define S1ALL_GC_MAX 4        //maximum number of GC pairs allowed in stem 1
#define S2BOT_GC_LIMIT 1      //maximum number of GC pairs allowed at the base of stem 2
#define S2TOP_GC_LIMIT 1      //maximum number of GC pairs allowed at the top of stem 2
#define BRCD_GC_MIN 20        //minimum GC content
#define BRCD_GC_MAX 24        //maximum GC content

/* barcode bank: storage for passed filter barcodes during barcode generation. bank is a doubly linked list */
typedef struct barcode_bank {
    fasta *fa;                 //fasta struct pointer for barcode sequence storage allocation
    struct barcode_bank *nxt;  //pointer to next barcode bank
    struct barcode_bank *prev; //pointer to previous barcode bank
    int cnt;                   //number of barcodes in bank
} barcode_bank;

/* mk_brcds: coordinates barcode generation */
int mk_brcds(int brcds2mk);

/* init_brcd_tmplts: initialize variant templates for barcode sequences */
int init_brcd_tmplts(wt_source * brcd_src, basemap * brcd_bmap);

/* xpnd brcd_tmplts: expand barcode templates to generate all possible barcode variants */
uint64_t xpnd_brcd_tmplts(basemap * brcd_bmap);

/* filter_barcode: filter barcodes that do not meet the parameters specified in make_barcodes.h */
int filter_barcode(char * sq);

/* store_pf_brcds: store barcodes that passed filter in a targets struct*/
int store_pf_brcds(target * pf_brcds, uint64_t passed_filter);

/* get_rndom_brcds: select random passed filter barcodes to output*/
int get_rndm_brcds(target ** brcd_out, target * pf_brcds, uint64_t passed_filter, int brcds2mk);

/* init_barcode_bank: initialize new barcode bank */
void init_barcode_bank(barcode_bank * bc_bnk, barcode_bank * prev);

/* add_barcode_bank: allocate a new barcode bank within the barcode bank linked list */
void add_barcode_bank(barcode_bank * bc_bnk);

/* merge_barcode_bank: merge the barcode bank linked list into a fasta array */
void merge_barcode_bank(barcode_bank * bc_bnk, uint64_t cnt);

#endif /* barcode_variants_h */

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

#define BARCODE_TEMPLATES 1      //number of variant templates used for barcode generation

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

/* mk_brcds: coordinates barcode generation */
int mk_brcds(int brcds2mk);

/* init_brcd_tmplts: initialize variant templates for barcode sequences */
int init_brcd_tmplts(wt_source * brcd_src, basemap * brcd_bmap, barcode_template brcd_tmplt[BARCODE_TEMPLATES]);

/* xpnd brcd_tmplts: expand barcode templates to generate all possible barcode variants */
uint64_t xpnd_brcd_tmplts(basemap * brcd_bmap);

/* filter_barcode: filter barcodes that do not meet the parameters specified in make_barcodes.h */
int filter_barcode(char * sq);

/* store_pf_brcds: store barcodes that passed filter in a targets struct*/
int store_pf_brcds(target * pf_brcds, uint64_t passed_filter);

/* get_rndom_brcds: select random passed filter barcodes to output*/
int get_rndm_brcds(target ** brcd_out, target * pf_brcds, uint64_t passed_filter, int brcds2mk, barcode_template brcd_tmplt[BARCODE_TEMPLATES]);

/* init_barcode_bank: initialize new barcode bank */
void init_barcode_bank(barcode_bank * bc_bnk, barcode_bank * prev);

/* add_barcode_bank: allocate a new barcode bank within the barcode bank linked list */
void add_barcode_bank(barcode_bank * bc_bnk);

/* merge_barcode_bank: merge the barcode bank linked list into a fasta array */
void merge_barcode_bank(barcode_bank * bc_bnk, uint64_t cnt);

#endif /* barcode_variants_h */

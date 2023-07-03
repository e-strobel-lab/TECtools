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

#define BARCODE_TEMPLATES 4   //number of variant templates used for barcode generation
#define STEM1BOT_N_STRT 1     //start of stem 1 bottom N stretch
#define STEM1BOT_N_END  3     //  end of stem 1 bottom N stretch
#define STEM1TOP_N_STRT 4     //start of stem 1 top N stretch
#define STEM1TOP_N_END  6     //  end of stem 1 top N stretch

//barcode parameter definitions
#define HOMOPOL_LIMIT 2       //maximum length of a homopolymeric stretch
#define POS1TO3_GC_MIN 1      //maximum number of AT pairs allowed at the base of stem 1
#define POS1TO3_GC_MAX 2      //maximum number of GC pairs allowed at the base of stem 1
#define POS4TO6_GC_MIN 1      //maximum number of AT pairs allowed at the top of stem 1
#define POS4TO6_GC_MAX 2      //maximum number of GC pairs allowed at the top of stem 1
#define POS1TO6_GC_MIN 3      //minimum number of GC pairs allowed in stem 1
#define POS1TO6_GC_MAX 4      //maximum number of GC pairs allowed in stem 1
#define POS25TO26_GC_LIMIT 1  //maximum number of GC pairs allowed at the base of stem 2
#define BRCD_GC_MIN 19        //minimum GC content
#define BRCD_GC_MAX 23        //maximum GC content

/* mk_brcds: coordinates barcode generation */
int mk_brcds(int brcds2mk);

/* init_brcd_tmplts: initialize variant templates for barcode sequences */
int init_brcd_tmplts(wt_source * brcd_src, basemap * brcd_bmap);

/* xpnd brcd_tmplts: expand barcode templates to generate all possible barcode variants */
int xpnd_brcd_tmplts(basemap * brcd_bmap);

/* filter_barcode: filter barcodes that do not meet the parameters specified in make_barcodes.h */
int filter_barcode(char * sq);

/* store_pf_brcds: store barcodes that passed filter in a targets struct*/
int store_pf_brcds(target * pf_brcds, int passed_filter);

/* get_rndom_brcds: select random passed filter barcodes to output*/
int get_rndm_brcds(target ** brcd_out, target * pf_brcds, int passed_filter, int brcds2mk);

#endif /* barcode_variants_h */

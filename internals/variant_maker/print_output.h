//
//  print_output.h
//  
//
//  Created by Eric Strobel on 8/29/24.
//

#ifndef print_output_h
#define print_output_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"

#include "read_bcFile.h"

void print_output(names * nm, basemap * bmap, int vTmpCnt, int varCnt, char * out_dir, int append_priming, int append_barcode, FILE * fp_brcd, int first_bc_2_use, char * lnkr, int make_fasta);

/* print_standard_variant: print variant sequence that does not contain a barcode */
void print_standard_variant(FILE * out_fp, FILE * fasta_fp, int append_priming, int crrnt_var, int make_fasta);

void print_barcoded_variant(FILE * fp_brcd, FILE * out_fp, FILE * fasta_fp, int append_priming, char * cstm_lnkr, int vTmpCnt, int crrnt_var, int bcs_per_var, int first_bc_2_use, int make_fasta);

#endif /* print_output_h */

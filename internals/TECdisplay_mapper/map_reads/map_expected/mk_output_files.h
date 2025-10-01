//
//  mk_output_files.h
//  
//
//  Created by Eric Strobel on 5/3/23.
//

#ifndef mk_output_files_h
#define mk_output_files_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../TECdisplay_mapper_structs.h"
#include "../../TECdisplay_mapper_defs.h"

#include "../../TECdisplay_output_column_headers.h"

#include "../../../utils/io_management.h"
#include "../../../utils/gen_utils.h"

/* print_output_header: print data output file header line */
void print_output_header(FILE * out_fp, char * out_nm);

/* print_output: print data output file */
void print_output(target * trgts, target_params * trg_prms, TDSPLY_names * nm);

/* print_metrics: print read mapping metrics */
void print_metrics(target * trgts, target_params * trg_prms, TDSPLY_metrics * met, TDSPLY_names * nm);

/* print_navigator_template: make template files for TECdisplay_navigator */
void print_navigator_template(target *refs, TDSPLY_fasta * wt, target_params * trg_prms);

#endif /* mk_output_files_h */

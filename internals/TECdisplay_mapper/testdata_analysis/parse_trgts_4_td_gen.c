//
//  parse_trgts_4_td_gen.c
//  
//
//  Created by Eric Strobel on 7/9/26.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"

#include "../TECdisplay_mapper_defs.h"
#include "../TECdisplay_mapper_structs.h"

#include "../map_reads/map_expected/parse_vmt_trgts.h"

#include "parse_trgts_4_td_gen.h"

/* parse_trgts_4_td_gen: parse targets input file and store data in targets structures to be used for testdata generation. this function is used to circumvent the need to write a separate testdata generation function that is compatible with target sequences being stored in compact_target structures when using compact_target structures for read mapping */
void parse_trgts_4_td_gen(TDSPLY_names * nm, int trgt_ftype, target ** td_refs, opt_ref ** td_ref_val, target ** td_trgts, opt_mx_trg ** td_trg_val, target_params * td_trg_prms, TDSPLY_fasta * td_wt)
{
    FILE * fp_trgs = NULL; //pointer to targets file
    
    get_file(&(fp_trgs), nm->trgs);                  //open targets file
    parse_header_lines(fp_trgs, td_trg_prms, td_wt); //parse targets file header lines
    
    //allocate memory
    if (((*td_refs) = calloc(MAXREF, sizeof(**td_refs))) == NULL) {
        printf("parse_trgts_4_td_gen: error - reference target memory allocation failed\n");
        abort();
    }
    
    if (((*td_ref_val) = calloc(MAXREF, sizeof(**td_ref_val))) == NULL) {
        printf("parse_trgts_4_td_gen: error - reference target value memory allocation failed\n");
        abort();
    }
    
    if (((*td_trgts) = calloc(td_trg_prms->xpctd, sizeof(**td_trgts))) == NULL) {
        printf("parse_trgts_4_td_gen: error - target memory allocation failed\n");
        abort();
    }
    
    if (((*td_trg_val) = calloc(td_trg_prms->xpctd, sizeof(**td_trg_val))) == NULL) {
        printf("parse_trgts_4_td_gen: error - target value memory allocation failed\n");
        abort();
    }
    
    //store targets as standard targets and close targets file
    parse_vmt_trgts(fp_trgs, trgt_ftype, *td_refs, *td_ref_val, *td_trgts, *td_trg_val, td_trg_prms, td_wt, TDSPLY, STD_TRGT_STRUCT);
    fclose(fp_trgs); //close targets file
    
    return;
}

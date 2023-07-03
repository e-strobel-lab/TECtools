//
//  parse_constraints.h
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#ifndef parse_constraints_h
#define parse_constraints_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"
#include "../seq_utils/test_possible_pairs.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/basemap.h"
#include "../seq_utils/is_dgnrt_mtch.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"
#include "./read_vbase.h"
#include "./search_4_vbase_match.h"

/* parse_constraints: read constraints file and store the settings for each constraint in a constraints struct */
int parse_constraints(FILE * ipt, constraints * cons, basemap * bmap, char * cnstnt_indels);

/* set_cnstnt_indel_constraints: add constant indel constraints to the constraints structure base constraints */
void set_cnstnt_indel_constraints(constraints * cons, char * cnstnt_indels, int * bi, basemap * bmap);

/* chk_vbase_complete: check that there is a base constraint for every reference vbase */
void chk_vbase_complete(constraints * cons, basemap * bmap);

/* chk_vbase_duplicates: check that constraints do not contain duplicate vbase entries */
void chk_vbase_duplicates(constraints * cons);

/* chk_pair_vbases: check that bases in pair constraints match vbases
 in base constraints or constant bases in the reference sequence */
void chk_pair_vbases(constraints * cons, basemap * bmap);

/* chk_pair_duplicates: check that constraints do not contain duplicate pair entries */
void chk_pair_duplicates(constraints * cons);

/* validate_pairs: validate base pair constraints */
void validate_pairs(constraints * cons, basemap * bmap);

/* check for conflicting pair attributes in overlapping pairs */
//TODO: add pair overlap test later
//void chk_pair_overlap(constraints * cons);

/* print_vbase: print formatted vbase */
void print_vbase(base_params * vbase);

/* print_base_constraints: print list of all base constraints */
void print_base_constraints(constraints * cons);

#endif /* parse_constraints_h */

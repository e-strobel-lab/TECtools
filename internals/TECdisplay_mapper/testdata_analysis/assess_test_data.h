//
//  assess_test_data.h
//  
//
//  Created by Eric Strobel on 5/3/23.
//

#ifndef assess_test_data_h
#define assess_test_data_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../TECdisplay_mapper_structs.h"
#include "../TECdisplay_mapper_defs.h"

#include "../../utils/gen_utils.h"

/* parse_testdata_id: parse test data id line and store attributes */
void parse_testdata_id(testdata_vars * testdata, char **td_trg_id, int * crnt_mut_cd, char * id_line);

/* evaluate_testdata_mtch: check that testdata read mapped to expected target */
int eval_testdata_mtch(testdata_vars * testdata, char * td_trg_id, int crnt_mut_cd, char * end5p, h_node **p_rdnd);

/* print_testdata_analysis: assess testdata analysis
 outcomes and print results to a file and to the screen */
void print_testdata_analysis(metrics *met, testdata_vars * testdata, target_params * trg_prms, target * trgts);

#endif /* assess_test_data_h */

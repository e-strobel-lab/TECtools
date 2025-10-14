//
//  testdata3pEnd_analysis.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef testdata3pEnd_analysis_h
#define testdata3pEnd_analysis_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/gen_utils.h"
#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/mapping_metrics.h"

#include "testdata3pEnd_analysis.h"

#define TESTDATA_NAT 0	//flag that test data reads contain native 3' ends
#define TESTDATA_RND 1	//flag that test data reads contain 3' ends with a 1 nt substitution/indel

struct testdata_3pEnd_vars {	//structure containing variables for testdata metrics
    int run;					//flag indicating whether to perform testdata analysis
    int mode;					//flag indicating whether test data reads contain native or randomized ends
    int end3p;					//current 3' end, used during testdata analysis
    int perfect[END_MAX];		//tracks reads in which the mapped 3' end matches the expected value
    int distance[END_MAX];		//tracks distribution of mapped:expected 3' end distances
    int max_distance;			//maximum mapped:expected 3' end distance
};

/* check_testdata_ipt: verify testdata analysis input files */
int check_testdata_ipt(TPROBE_names * nm);

/* get_testdata_3pEND: get 3' end info from testdata read id */
void get_testdata_3pEND(char * id);

void compare_testdata_3pEnd(int end, char id[MAX_LINE], char sq[MAX_LINE]);

/* print_3pEnd_testdata_analysis: assess testdata analysis
 outcomes and print results to a file and to the screen */
void print_3pEnd_testdata_analysis(mapping_metrics * met, target3p_params trg_prms, target * trgts); //print report of test data analysis

#endif /* testdata3pEnd_analysis_h */

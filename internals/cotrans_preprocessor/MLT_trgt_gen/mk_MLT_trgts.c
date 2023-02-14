//
//  mk_MLT_trgts.c
//  
//
//  Created by Eric Strobel on 3/16/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/io_management.h"
#include "../../seq_utils/parse_fasta.h"
#include "mk_intermed_trgts.h"
#include "mk_3pEnd_trgts.h"
#include "mk_3pEnd_testdata.h"

#include "mk_MLT_trgts.h"

/* mk_MLT_trgts: manage intermediate and 3' end target generation */
int mk_MLT_trgts(FILE * fp_fasta, int usr_min_len, int usr_end_len, int make_test_data, int mltplir, int randomize_end){
    
    char fa_name[MAX_LINE] = {0};	//name from fasta file
    char seq[MAX_LINE] = {0};		//input sequence
    
    parse_fasta(fp_fasta, &fa_name[0], &seq[0]); //parse fasta file
    
    //test that target sequence length is permissible
    if (strlen(seq) > END_MAX) {
        printf("mk_MLT_trgts: error - target sequence length (%lu) exceeds allowed maximum (%d). aborting...\n", strlen(seq), END_MAX);
    }
    
    target3p_genVals ends = {{{0}}};	//structure containing end target variables
    int min_len = 0;					//minimum transcript intermediate-1
    
    //print target name and sequence to screen
    printf("\nname: %s\n", fa_name);
    printf(" seq: %s\n", seq);
    
    set_min_len(usr_min_len, &min_len);	//set minimum target length
    set_end_len(usr_end_len, &ends);	//set minimum end length
    
    
    /******** check that minimum length end length values are permissible *******/
    if (ends.len > min_len) {
        printf("mk_ends_file: error - end target length (%d) cannot exceed minimum transcript intermediate target length (%d). aborting...\n", ends.len, min_len);
        abort();
    }
    
    if (min_len > strlen(seq)) {
        printf("mk_ends_file: error - minimum target length (%d) cannot exceed input sequence length (%lu). aborting...\n", min_len, strlen(seq));
        abort();
    }
    /****** end check that minimum length end length values are permissible ******/
    
    
    /******** make output directories *******/
    //+35 accounts for extra chars printed, +1 accounts for terminating null
    char out_dir[MAX_LINE+35+1] = {0};
    
    sprintf(out_dir, "./%s_targets", fa_name);
    mk_out_dir(out_dir);
    
    sprintf(out_dir, "./%s_targets/intermediate_transcripts", fa_name);
    mk_out_dir(out_dir);
    
    sprintf(out_dir, "./%s_targets/metrics", fa_name);
    mk_out_dir(out_dir);
    /****** end make output directories ******/
    
    mk_intermed_trgts(fa_name, seq, min_len);		//generate fasta files for intermediate transcripts
    mk_3pEnd_trgts(&ends, fa_name, seq, min_len);	//generate 3' ends target file
    
    //generate test data
    if (make_test_data) {
        sprintf(out_dir, "./%s_targets/test_data", fa_name);
        mk_out_dir(out_dir);
        
        srand(time(NULL));
        mk_3pEnd_testdata(fa_name, seq, &ends, randomize_end, mltplir, min_len);
    }
    
    return 1;
}



/* set_min_len: the length of the shortest target*/
int set_min_len(int usr_min_len, int * min_len)
{
    if (usr_min_len) {
        (*min_len) = usr_min_len-1; //min_len is 1 less than the shortest target that will be generated
        printf("using user-specified minimum transcript length: %d nts\n\n", *min_len+1);
    } else {
        (*min_len) = DFLT_MIN;
        printf("using default minimum transcript length: %d nts\n\n", *min_len+1);
    }
    return 1;
}



/* set_end_len: set the length of 3' end target sequences */
int set_end_len(int usr_end_len, target3p_genVals * ends)
{
    if (usr_end_len) {				//check whether end length was specified
        ends->len = usr_end_len;	//set to user specified end length
        printf("\nusing user-specified end target length: %d nts\n", ends->len);
    } else {
        ends->len = END_LEN_INIT;	//set to default end length (14)
        printf("\nusing default end target length: %d nts\n", ends->len);
    }
    return 1;
}





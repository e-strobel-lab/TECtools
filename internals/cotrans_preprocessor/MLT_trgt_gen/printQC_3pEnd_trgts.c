//
//  printQC_3pEnd_trgts.c
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/gen_utils.h"

#include "printQC_3pEnd_trgts.h"

/* print_native_end_QC: prints table that contains QC information for each native end target */
int print_native_end_QC(FILE * metrics_fp, target3p_genVals * ends,
                        int len, int fndMtch, int mtchPos, int mode)
{
    int i = 0;
    
    char out_str[8192] = {0}; //string for storing output line
    
    //print native end QC table title  header line
    if (mode == HEADER) {
        
        //print title and start of table header line
        sprintf(out_str, "*** Native End QC Table ***\n3' end\t| sequence");
        printf2_scrn_n_fl(metrics_fp, out_str);
        
        //print spaces to account for variation in sequence length
        for (i = 0; i < (ends->len-8); i++) {
            printf(" ");
            fprintf(metrics_fp, " ");
        }
        
        //print the rest of the table header line
        sprintf(out_str, "\t| matches | exact | pass?\n");
        printf2_scrn_n_fl(metrics_fp, out_str);
        
    } else if (mode == LINE) { //print native end QC table data line
        
        //print 3' end value, target sequence, match count, and match position information
        sprintf(out_str, "%6d\t| %s\t| %7d | %5d | ", len, ends->sq[ends->nat], fndMtch, mtchPos);
        printf2_scrn_n_fl(metrics_fp, out_str);
        
        //print pass/fail outcome
        if (mtchPos == len) { //test that match was at the expected length
            if (fndMtch == 1) {				  //only 1 perfect match, pass
                sprintf(out_str, "    Y\n");
                
            } else if (fndMtch > 1) {		  //more than 1 perfect match, fail
                sprintf(out_str, " FAIL, try increasing end target length to >14\n");
                
            } else if (fndMtch == 0) {		  //no match to target in input sequence
                printf("\nprint_native_end_QC: error - cannot find target in input sequence, a serious error has occurred.\naborting...\n");
                abort();
                
            } else {
                printf("\nprint_native_end_QC: error - this code should be unreachable, a serious error has occurred.\naborting...\n");
                abort();
            }
        } else {
            sprintf(out_str, " FAIL, perfect match is not at the expected length\n");
        }
        printf2_scrn_n_fl(metrics_fp, out_str);
        
    } else {
        printf("print_native_end_QC: error - invalid mode. aborting...\n");
        abort();
    }
    
    return 1;
}

/* print_match_distance_table: prints table that contains summary data
 describing distance of each end target from all other end targets*/
int print_match_distance_table(FILE * metrics_fp, target3p_genVals * ends,  int total_tested)
{
    int i = 0;
    
    char out_str[8192] = {0}; //string for storing output line
    
    //print table title and header line
    sprintf(out_str, "\n*** Match Distance Table ***\ndistance| counts\n");
	printf2_scrn_n_fl(metrics_fp, out_str);
    
    for (i = 0; i <= ends->len; i++) {
        
        //print distance and counts columns (all data lines)
        sprintf(out_str, "%8d: %6d", i, ends->dstnc[i]);
        printf2_scrn_n_fl(metrics_fp, out_str);
        
        //print expected number of distance zero end targets
        //in distance zero table line
        if (i == 0) {
            sprintf(out_str, " (expected %d)\n",total_tested);	//print expected 0 distance
        } else {
            sprintf(out_str, "\n");	//all other data lines just get a newline
        }
        printf2_scrn_n_fl(metrics_fp, out_str);
    }
    
    return 1;
}

/* print_end_target_table: prints table that contains the number of end targets generated */
int print_end_target_table(FILE * metrics_fp, target3p_genVals * ends)
{
    char out_str[8192] = {0}; //string for storing output line
    
    sprintf(out_str, "\n*** End Target Metrics Table ***\n    end type | generated\n"); //title and header
    printf2_scrn_n_fl(metrics_fp, out_str);
    
    sprintf(out_str, "      native | %9d\n", ends->act[NAT]); //number of native seqs
    printf2_scrn_n_fl(metrics_fp, out_str);
    
    sprintf(out_str, "substitution | %9d\n", ends->act[SUB]); //number of subsitution seqs
    printf2_scrn_n_fl(metrics_fp, out_str);
    
    sprintf(out_str, "   insertion | %9d\n", ends->act[INS]); //number of insertion seqs
    printf2_scrn_n_fl(metrics_fp, out_str);
    
    sprintf(out_str, "    deletion | %9d\n", ends->act[DEL]); //number of deletion seqs
    printf2_scrn_n_fl(metrics_fp, out_str);
    
    return 1;
}

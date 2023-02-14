//
//  mk_3pEnd_trgts.c
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "printSubs.h"
#include "printInsrts.h"
#include "printDels.h"
#include "printQC_3pEnd_trgts.h"

#include "mk_3pEnd_trgts.h"

/* mk_3pEnd_trgts: generate a file containing every possible 3' end
 of length len, including native sequences, 1 nt substitutions,
 1 nt insertions, and 1 nt deletions
 */
int mk_3pEnd_trgts(target3p_genVals * ends, char name[MAX_LINE], char seq[MAX_LINE], int min_len)
{
    //open output file
    FILE * out_fp = {NULL};			//output file pointer
    char out_name[MAX_LINE] = {0};	//output file name
    
    sprintf(out_name, "./%s_targets/%s_ends_%02dntLong.txt", name, name, ends->len);
    if ((out_fp = fopen(out_name,"w")) == NULL) {
        printf("mk_3pEnd_trgts: error - could not open ends output file. Aborting program...\n");
        abort();
    }
    
    //open metrics file
    FILE * metrics_fp = {NULL};
    char metrics_name[MAX_LINE] = {0};
    
    sprintf(metrics_name, "./%s_targets/metrics/%s_metrics.txt", name, name);
    if ((metrics_fp = fopen(metrics_name,"w")) == NULL) {
        printf("mk_3pEnd_trgts: error - could not open metrics file. Aborting program...\n");
        abort();
    }
    
    //print input name and sequence
    printf("\n\n>inputs\nname: %s\nseq: %s\n", name, seq);
    fprintf(metrics_fp, ">inputs\nname: %s\nseq: %s\n", name, seq);
    
    //general variables
    int i = 0;
    int j = 0;
    int k = 0;
    
    int ipt_len = strlen(seq);	//length of the input sequence
    int MMCnt = 0;				//used to assess native 3' end target ambiguity
    int fndMtch = 0;			//flag to indicate end target match when searching seq string
    int mtchPos = 0;			//position of end target match in seq string
    
    //print input sequence length
    printf("ipt len = %d\n\n", ipt_len);
    fprintf(metrics_fp, "ipt len = %d\n\n", ipt_len);
    
    //calculate expected end counts
    //see end target generation functions for an explanation
    //of the expected end count calculations
    ends->exp[NAT] = ipt_len - min_len;
    ends->exp[SUB] = ends->exp[NAT] * ends->len * 3;
    ends->exp[INS] = (ends->exp[NAT] * 6) + (ends->exp[NAT] * (ends->len-2) * 4) + 1;
    ends->exp[DEL] = ends->exp[NAT] * (ends->len - 1);
    
    
    /* generate native end targets and assess distance
     from all other possible native end targets
     
     native matches are simply the last 14 nt of
     each transcript length starting at length 30
     */
    print_native_end_QC(metrics_fp, ends, i, fndMtch, mtchPos, HEADER);
    for (i = min_len, ends->nat = 0; i <= ipt_len; i++, ends->nat++) {
        
        /* copy all posible native ends to ends array */
        for (j = (i-ends->len), k = 0; j < i && k < MAX_END_LEN && seq[j]; j++, k++) {
            ends->sq[ends->nat][k] = toupper(seq[j]);    //ensure all sequences are uppercase
        }
        ends->sq[ends->nat][k] = '\0';
        
        /* scan end sequence accross entire target sequence to assess
         whether it is sufficiently different from all other possible ends */
        //Note: the SC1 sequenced is not included in this check because the maximum
        //allowed target length (20) is the same as the minimum transcript length (20).
        //Only transcripts shorter than 20 nt require a 3' end target that overlaps SC1
        //to be mappable.
        for (j = 0, fndMtch = 0, mtchPos = 0; j <= ipt_len-ends->len; j++) {
            for (k = 0, MMCnt = 0; k < ends->len; k++) {
                if (ends->sq[ends->nat][k] != toupper(seq[j+k])) {
                    MMCnt++;
                }
            }
            ends->dstnc[MMCnt]++; //increment distance table
            
            //check for perfect match
            if (MMCnt == 0) {
                if (mtchPos != 0) { //more than one perfect match
                    printf("error - more than one perfect match for %d nt end target\n", i);
                    fprintf(metrics_fp, "error - more than one perfect match for %d nt end target\n", i);
                    mtchPos = 99999; //error code for more than one perfect match
                } else { //first perfect match
                    mtchPos = (j+ends->len);
                }
            }
            
            //check if distance from perfect match is <MATCHDIF threshold
            if (MMCnt < MATCHDIF) {	//if MMCnt is less than MATCHDIF
                fndMtch++;			//increment fndMtch
            }
        }
        
        //print native end QC table line
        print_native_end_QC(metrics_fp, ends, i, fndMtch, mtchPos, LINE);
    }
    
    //print match distance table. i-min_len is total tested variants
    print_match_distance_table(metrics_fp, ends, i-min_len);
    
    //print native end targets to output file
    //skip min_len end, which is only used when making deletion variants
    for (i = 1; i < ends->nat; i++) {
        fprintf(out_fp, "%s_%03d_NAT\t%s\n", name, min_len+i, ends->sq[i]);
    }
    ends->act[NAT] = i-1; //subtract 1 since i started at 1 because min_len is omitted
    
    //print sub, ins, and del end targets to output file
    ends->act[SUB] = printSubs(ends, name, out_fp, min_len);
    ends->act[INS] = printInsrts(ends, name, out_fp, min_len);
    ends->act[DEL] = printDels(ends, name, out_fp, min_len);
    
    //check that actual and expected end target counts match
    if (ends->act[NAT] != ends->exp[NAT] ||
        ends->act[SUB] != ends->exp[SUB] ||
        ends->act[INS] != ends->exp[INS] ||
        ends->act[DEL] != ends->exp[DEL] ) {
        printf("mk_3pEnd_trgts: error - the actual number of 3' ends generated does not match the expected value. aborting...");
        abort();
    }
    
    //print end target table
    print_end_target_table(metrics_fp, ends);
    
    //close output files
    if ((fclose(out_fp)) == EOF) {
        printf("mk_3pEnd_trgts: error - error occurred when closing ends file. Aborting program...\n");
        abort();
    }
    
    if ((fclose(metrics_fp)) == EOF) {
        printf("mk_3pEnd_trgts: error - error occurred when closing metrics file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

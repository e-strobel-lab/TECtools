//
//  printDels.c
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

#include "printDels.h"

/* printDels: generate all possible 1 nt deletions for each end target
 
 starting at length 30, delete each nucleotide and
 append the first character of the preceding target
 to the start of the current target
 
 it is not possible to delete the last nucleotide in
 a target because this would appear identical to the
 preceding target
 
 calculated value: (ipt_len - min_len) * (len - 1)
 */
int printDels(target3p_genVals * ends, char name[MAX_LINE], FILE * out_fp, int min_len)
{
    int i = 0;
    int j = 0;
    int k = 0;
    
    int pos2Var = 0;					//current position in target to be deleted
    int del_count = 0;					//total deletion targets generated
    char var[MAX_END_LEN+2] = {0};		//array for constructing deletion targets
    
    
    for (i = 1; i < ends->nat; i++) {	//i starts at 1 to skip first length
        
        //each iteration of the pos2var loop appends the 1st nt of the
        //preceding target to the current target, then copies the current
        //target sequence but skips the nt with index pos2var
        
        //the upper bound for the pos2var loop is len-1 because deleting
        //the last nucleotide will generate a target with an identical
        //sequence to that of the preceding target
        
        for (pos2Var = 0; pos2Var < ends->len-1 && pos2Var < MAX_END_LEN; pos2Var++) {
            var[0] = ends->sq[i-1][0]; //set var[0] as the 1st nt of the preceding target
            
            //copy end target, but skip pos2var to make deletion
            for (j = 0, k= 1; j < ends->len && k < ends->len && j < MAX_END_LEN; j++) {
                if (j != pos2Var) { //copy ends->sq[i][j] to var[k] except when j == pos2var
                    var[k++] = ends->sq[i][j];
                }
            }
            var[k] = '\0';
            
            //check that j and k are equal to ends->len, which indicates that pos2var was skipped
            if (j != ends->len || k != ends->len) {
                printf("printDels: error - deletion variant was not generated correctly. aborting...\n");
                abort();
            }
            
            //print deletion end target to output file
            fprintf(out_fp, "%s_%03d_DEL%02d_%d\t%s\n", name, min_len+i, pos2Var, 0, var);
            del_count++;
        }
    }
    
    return del_count;
}



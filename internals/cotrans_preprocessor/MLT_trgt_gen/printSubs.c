//
//  printSubs.c
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

#include "printSubs.h"

/* printSubs: generate all possible 1 nt substitutions for each end target
 
 starting at length 30, substitute each nucleotide
 with the other three possible nucleotides
 
 calculated value: (ipt_len - min_len) * len * 3
 */
int printSubs(target3p_genVals * ends, char name[MAX_LINE], FILE * out_fp, int min_len)
{
    int i = 0;
    int j = 0;
    
    int pos2Var = 0;					//current position in target to be varied
    int sub_count = 0;					//total substitution targets generated
    char var[3][MAX_END_LEN+2] = {{0}};	//array for constructing substitution targets
    
    for (i = 1; i < ends->nat; i++) {	//i starts at 1 to skip first length
        
        //each iteration of the pos2var loop generates three substitution
        //variants of the current target at position pos2var
        
        for (pos2Var = 0; pos2Var < ends->len && pos2Var < MAX_END_LEN; pos2Var++) {
            for (j = 0; j < ends->len && j < MAX_END_LEN; j++) {	//for each nt in the end target sequence
                if (j == pos2Var) {				//if index j is being varied in this interation
                    switch (ends->sq[i][j]) {	//substitute sq[i][j] with the other 3 bases in var array
                        case 'A':
                            var[0][j] = 'T';
                            var[1][j] = 'G';
                            var[2][j] = 'C';
                            break;
                        case 'T':
                            var[0][j] = 'A';
                            var[1][j] = 'G';
                            var[2][j] = 'C';
                            break;
                        case 'G':
                            var[0][j] = 'A';
                            var[1][j] = 'T';
                            var[2][j] = 'C';
                            break;
                        case 'C':
                            var[0][j] = 'A';
                            var[1][j] = 'T';
                            var[2][j] = 'G';
                            break;
                        default:
                            printf("printSubs: error - unexpected character %c in native 3' end sequence. aborting...\n", ends->sq[i][j]);
                            abort();
                            break;
                    };
                } else { //not varying index j, copy native base to var array
                    var[0][j] = ends->sq[i][j];
                    var[1][j] = ends->sq[i][j];
                    var[2][j] = ends->sq[i][j];
                }
            }
            var[0][j] = '\0';
            var[1][j] = '\0';
            var[2][j] = '\0';
            
            //print substitution end targets to output file
            for (j = 0; j < 3; j++) {
                fprintf(out_fp, "%s_%03d_SUB%02d_%d\t%s\n", name, min_len+i, pos2Var, j, var[j]);
                sub_count++;
            }
        }
    }
    
    return sub_count;
}

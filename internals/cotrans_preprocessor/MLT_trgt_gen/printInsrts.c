//
//  printInsrts.c
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

/* printInsrts: generate all possible 1 nt insertions for each end target
 
 internal insertions can be any of the four DNA bases because there is not
 a native nucleotide at this position (term: (ipt_len - min_len) * 12 * 4).
 leading and trailing insertions can only be the 3 non-native DNA bases
 because inserting the native base at the leading position will generate a
 sequence that appears the same as the native target and inserting the native
 base at the trailing position will appear the same as the next native target
 (term: (ipt_len - min_len) * 6). the exception to this rule is the trailing
 insertion of the last 3' end target (term: + 1).
 
 calculated value: ((ipt_len - min_len) * 6) + ((ipt_len - min_len) * 12 * 4) + 1
 */
int printInsrts(target3p_genVals * ends, char name[MAX_LINE], FILE * out_fp, int min_len)
{
    int i = 0;
    int j = 0;
    int k = 0;
    
    int pos2Var = 0;					//current position in target to add insertion
    int ins_count = 0;					//total insertion targets generated
    char var[4][MAX_END_LEN+2] = {{0}};	//array for constructing insertion targets
    
    int varied_first_pos = 0;
    int varied_nonFirst_pos = 0;
    
    for (i = 1; i < ends->nat; i++) {	//i starts at 1 to skip first length
        for (pos2Var = 0; pos2Var < ends->len && pos2Var < MAX_END_LEN; pos2Var++) {
            
            //j is initialized to 1 because the first nucleotide of the
            //end target will always be skipped due to the insertion
            
            //except for the first position of an end target, the native end
            //target nucleotide is copied to the variant arrays and then an
            //insertion is made if the j index is equal to pos2var
            
            //when varying the first position of a target, the insertion is
            //handled as a substitution (see below) before copying the native
            //end target nucleotide to the variant arrays
            
            for (j = 1, k= 0; j < ends->len && k < ends->len && j < MAX_END_LEN; j++, k++) {
                varied_first_pos = varied_nonFirst_pos = 0;
                
                //if varying the first position of the target:
                //the insertion is effectively a substitution
                //and can be treated as such
                
                //logic for generating first nucleotide insertion variants
                //this code must precede copying a native nt to the var arrary;
                //all other insertions occur after the native nt is copied
                if (pos2Var == 0 && j == 1) {	//only use when insertion is at first position of the target
                    var[3][0] = '\0';			//only three variants will result, set fourth to zero
                    switch (ends->sq[i][0]) {
                        case 'A':
                            var[0][k] = 'T';
                            var[1][k] = 'G';
                            var[2][k] = 'C';
                            break;
                        case 'T':
                            var[0][k] = 'A';
                            var[1][k] = 'G';
                            var[2][k] = 'C';
                            break;
                        case 'G':
                            var[0][k] = 'A';
                            var[1][k] = 'T';
                            var[2][k] = 'C';
                            break;
                        case 'C':
                            var[0][k] = 'A';
                            var[1][k] = 'T';
                            var[2][k] = 'G';
                            break;
                        default:
                            printf("printInsrts: error - unexpected character %c in native 3' end sequence. aborting...\n", ends->sq[i][0]);
                            abort();
                            break;
                    };
                    k++; //increment to next position in var array
                    varied_first_pos = 1; //set flag that insertion was made at first position
                }
                
                //copy native nucleotide to variant arrays
                //when varying the first position, this code executes after the insertion was made
                //when varying all other positions, this code executes before the insertion is made
                //
                //sequence is copied to var[3] even if it will be discarded later. when var[3] will
                //be discarded, var[3][0]='\0'
                var[0][k] = ends->sq[i][j];
                var[1][k] = ends->sq[i][j];
                var[2][k] = ends->sq[i][j];
                var[3][k] = ends->sq[i][j];
                
                //logic for generating non-first insertion variants. this
                //should never be executed in the same loop interation as
                //the logic for first nucleotide insertion variants. an
                //error is thrown if this occurs.
                
                if (j == pos2Var) {	//at the position where the insertion will be made
                    k++;			//increment the variant arrary
                    
                    //if the insertion occurs at the last target nucleotide
                    //and the target is not the last possible variant, only
                    //three insertion variants will be generated because
                    //generating the fourth is guaranteed to yield an identical
                    //sequence to the next native target
                    
                    if (((j+1) == ends->len) && ((i+1) != ends->nat) ) {
                        var[3][0] = '\0';
                        
                        //make insertion that doesn't match the
                        //last position of the next native target
                        switch (ends->sq[i+1][j]) {
                            case 'A':
                                var[0][k] = 'T';
                                var[1][k] = 'G';
                                var[2][k] = 'C';
                                break;
                            case 'T':
                                var[0][k] = 'A';
                                var[1][k] = 'G';
                                var[2][k] = 'C';
                                break;
                            case 'G':
                                var[0][k] = 'A';
                                var[1][k] = 'T';
                                var[2][k] = 'C';
                                break;
                            case 'C':
                                var[0][k] = 'A';
                                var[1][k] = 'T';
                                var[2][k] = 'G';
                                break;
                            default:
                                printf("printInsrts: error - unexpected character %c in native 3' end sequence. aborting...\n", ends->sq[i+1][j]);
                                abort();
                                break;
                        };
                    } else {
                        //the insertion is internal to the end target,
                        //so four insertion variants are generated
                        
                        var[0][k] = 'A';
                        var[1][k] = 'T';
                        var[2][k] = 'G';
                        var[3][k] = 'C';
                    }
                    varied_nonFirst_pos = 1;
                }
                
                //this error should never occur
                if (varied_first_pos && varied_nonFirst_pos) {
                    printf("printInserts: error - varied first and non-first position within the same loop iteration. aborting...");
                    abort();
                }
            }
            var[0][j] = '\0';
            var[1][j] = '\0';
            var[2][j] = '\0';
            var[3][j] = '\0';
            
            //check that j and k are equal to ends->len, which indicates that an insertion was made
            if (j != ends->len || k != ends->len) {
                printf("printInsrts: error - insertion variant was not generated correctly. aborting...\n");
                abort();
            }
            
            //print insertion variants to file
            for (j = 0; j < 3; j++) {
                fprintf(out_fp, "%s_%03d_INS%02d_%d\t%s\n", name, min_len+i, pos2Var, j, var[j]);
                ins_count++;
            }
            if (var[3][0]) {
                fprintf(out_fp, "%s_%03d_INS%02d_%d\t%s\n", name, min_len+i, pos2Var, j, var[3]);
                ins_count++;
            }
        }
    }
    
    return ins_count;
}

#include "printInsrts.h"

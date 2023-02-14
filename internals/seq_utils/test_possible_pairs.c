//
//  test_possible_pairs.c
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#include <stdio.h>
#include <stdlib.h>

#include "test_possible_pairs.h"

/* test_possible_pairs: tests possible interactions between bs1 and bs2 using a lookup table.
 can be run in BASEPAIR mode to test for possible pairing interactions, or in MISMATCH mode
 to test whether mismatches are possible.
 */
int test_possible_pairs(char bs1, char bs2, int mode)
{
    char *p_ipt[2] = {NULL}; //pointer array for testing bs1 and bs2 identity in a loop
    
    p_ipt[0] = &bs1; //set index 0 pointer to point to bs1
    p_ipt[1] = &bs2; //set index 1 pointer to point to bs2
    
    int index[2] = {-1, -1}; //array for setting indices based on bs1 and bs2 identity
    
    int i = 0;
    
    /*** set indices based on the identity of bs1 and bs2 ***/
    //indices are used in the lookup tables below
    for (i = 0; i < 2; i++) {
        switch (*(p_ipt[i])) {
            case 'A': index[i] = 0; break;
            case 'T': index[i] = 1; break;
            case 'G': index[i] = 2; break;
            case 'C': index[i] = 3; break;
            case 'R': index[i] = 4; break;
            case 'Y': index[i] = 5; break;
            case 'M': index[i] = 6; break;
            case 'K': index[i] = 7; break;
            case 'S': index[i] = 8; break;
            case 'W': index[i] = 9; break;
            case 'B': index[i] = 10; break;
            case 'D': index[i] = 11; break;
            case 'H': index[i] = 12; break;
            case 'V': index[i] = 13; break;
            case 'N': index[i] = 14; break;
                
            default:
                printf("test_for_pair: error - unexpected input %c. aborting...\n", *(p_ipt[i]));
                abort();
                break;
        }
    }
    
    //TODO: complete check of table generation scripts
    int pair_lkup[15][15] = {
        //lookup table for base pairing potential. value indicates
        //which (if any) pairing interactions are possible
        //0 = NO PAIR
        //1 = WEAK PAIR
        //2 = STRONG PAIR
        //3 = WEAK AND STRONG PAIR
        
        //      A  T  G  C  R  Y  M  K  S  W  B  D  H  V  N
        /*A*/ { 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1 },
        /*T*/ { 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*G*/ { 0, 1, 0, 2, 0, 3, 2, 1, 2, 1, 3, 1, 3, 2, 3 },
        /*C*/ { 0, 0, 2, 0, 2, 0, 0, 2, 2, 0, 2, 2, 0, 2, 2 },
        /*R*/ { 0, 1, 0, 2, 0, 3, 2, 1, 2, 1, 3, 1, 3, 2, 3 },
        /*Y*/ { 1, 0, 3, 0, 3, 0, 1, 3, 3, 1, 3, 3, 1, 3, 3 },
        /*M*/ { 0, 1, 2, 0, 2, 1, 0, 3, 2, 1, 3, 3, 1, 2, 3 },
        /*K*/ { 1, 1, 1, 2, 1, 3, 3, 1, 3, 1, 3, 1, 3, 3, 3 },
        /*S*/ { 0, 1, 2, 2, 2, 3, 2, 3, 2, 1, 3, 3, 3, 2, 3 },
        /*W*/ { 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*B*/ { 1, 1, 3, 2, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3, 3 },
        /*D*/ { 1, 1, 1, 2, 1, 3, 3, 1, 3, 1, 3, 1, 3, 3, 3 },
        /*H*/ { 1, 1, 3, 0, 3, 1, 1, 3, 3, 1, 3, 3, 1, 3, 3 },
        /*V*/ { 0, 1, 2, 2, 2, 3, 2, 3, 2, 1, 3, 3, 3, 2, 3 },
        /*N*/ { 1, 1, 3, 2, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3, 3 }
    };
    
    int MM_lkup[15][15] = {
        //lookup table for mismatch potential.
        //0 = CANNOT MISMATCH
        //1 = CAN MISMATCH

        //      A  T  G  C  R  Y  M  K  S  W  B  D  H  V  N
        /*A*/ { 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*T*/ { 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*G*/ { 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*C*/ { 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*R*/ { 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*Y*/ { 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*M*/ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*K*/ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*S*/ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*W*/ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*B*/ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*D*/ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*H*/ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*V*/ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        /*N*/ { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
    };
    
    if (mode == TEST_MISMATCH) {
        return MM_lkup[index[0]][index[1]];
    } else if (mode == TEST_BASEPAIR) {
        return pair_lkup[index[0]][index[1]];
    } else {
        printf("test_possible_pairs: error - unrecognized mode\n. aborting...");
        abort();
    }
}

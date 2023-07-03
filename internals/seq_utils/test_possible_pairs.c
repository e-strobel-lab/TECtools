//
//  test_possible_pairs.c
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "test_possible_pairs.h"

/* test_possible_pairs: tests possible interactions between bs1 and bs2 using a lookup table.
 can be run in BASEPAIR mode to test for possible pairing interactions, or in MISMATCH mode
 to test whether mismatches are possible.
 */
int test_possible_pairs(char bs1, char bs2, int mode, int type_test)
{
    if ((mode == TEST_BASEPAIR || mode == TEST_MISMATCH) && type_test != NO_TYPE_TEST) {
        printf("test_possible_pairs: error - type_test variable was not set to NO_TYPE_TEST in TEST_BASEPAIR or TEST_MISMATCH mode. aborting...\n");
        abort();
    }
    
    if (mode == TEST_PAIRTYPE && (type_test < 0 || type_test > 7)) {
        printf("test_possible_pairs: error - test value (%d) is out of bounds. test value must be >= %d or <= %d. aborting...\n", type_test, TEST_NO_PAIR, TEST_GT_AT_GC);
        abort();
    }
    
    char *p_ipt[2] = {NULL}; //pointer array for testing bs1 and bs2 identity in a loop
    
    p_ipt[0] = &bs1; //set index 0 pointer to point to bs1
    p_ipt[1] = &bs2; //set index 1 pointer to point to bs2
    
    int index[2] = {-1, -1}; //array for setting indices based on bs1 and bs2 identity
    
    int i = 0;
    
    /*** set indices based on the identity of bs1 and bs2 ***/
    //indices are used in the lookup tables below
    for (i = 0; i < 2; i++) {
        switch (toupper(*(p_ipt[i]))) {
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
                printf("test_possible_pairs: error - unexpected input %c. aborting...\n", *(p_ipt[i]));
                abort();
                break;
        }
    }
    
    int pair_type_lkup[15][15] = {
        //lookup table for base pairing potential. value indicates
        //which (if any) pairing interactions are possible
        //0 = NO PAIR
        //1 = GT PAIR
        //2 = AT_PAIR
        //3 = GT and AT PAIRS
        //4 = GC PAIR
        //5 = GT and GC PAIRS
        //6 = AT and GC PAIRS
        //7 = GT, AT, and GC PAIRS
        
        //      A  T  G  C  R  Y  M  K  S  W  B  D  H  V  N
        /*A*/ { 0, 2, 0, 0, 0, 2, 0, 2, 0, 2, 2, 2, 2, 0, 2 },
        /*T*/ { 2, 0, 1, 0, 3, 0, 2, 1, 1, 2, 1, 3, 2, 3, 3 },
        /*G*/ { 0, 1, 0, 4, 0, 5, 4, 1, 4, 1, 5, 1, 5, 4, 5 },
        /*C*/ { 0, 0, 4, 0, 4, 0, 0, 4, 4, 0, 4, 4, 0, 4, 4 },
        /*R*/ { 0, 3, 0, 4, 0, 7, 4, 3, 4, 3, 7, 3, 7, 4, 7 },
        /*Y*/ { 2, 0, 5, 0, 7, 0, 2, 5, 5, 2, 5, 7, 2, 7, 7 },
        /*M*/ { 0, 2, 4, 0, 4, 2, 0, 6, 4, 2, 6, 6, 2, 4, 6 },
        /*K*/ { 2, 1, 1, 4, 3, 5, 6, 1, 5, 3, 5, 3, 7, 7, 7 },
        /*S*/ { 0, 1, 4, 4, 4, 5, 4, 5, 4, 1, 5, 5, 5, 4, 5 },
        /*W*/ { 2, 2, 1, 0, 3, 2, 2, 3, 1, 2, 3, 3, 2, 3, 3 },
        /*B*/ { 2, 1, 5, 4, 7, 5, 6, 5, 5, 3, 5, 7, 7, 7, 7 },
        /*D*/ { 2, 3, 1, 4, 3, 7, 6, 3, 5, 3, 7, 3, 7, 7, 7 },
        /*H*/ { 2, 2, 5, 0, 7, 2, 2, 7, 5, 2, 7, 7, 2, 7, 7 },
        /*V*/ { 0, 3, 4, 4, 4, 7, 4, 7, 4, 3, 7, 7, 7, 4, 7 },
        /*N*/ { 2, 3, 5, 4, 7, 7, 6, 7, 5, 3, 7, 7, 7, 7, 7 }
    };
    

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
    
    
    if (mode == TEST_MISMATCH) {              //in TEST_MISMATCH mode
        return MM_lkup[index[0]][index[1]];   //return MM_lkp value
        
    } else if (mode == TEST_BASEPAIR) {       //in TEST_BASEPAIR mode
        return pair_lkup[index[0]][index[1]]; //return pair_lkup value
        
    } else if (mode == TEST_PAIRTYPE) {       //in TEST_PAIRYPE mode
        
        if (pair_type_lkup[index[0]][index[1]] == type_test || //if lookup val equals type test (required for NO_PAIR)
            pair_type_lkup[index[0]][index[1]] & type_test) {  //or type_test masked lookup value is TRUE
            return 1;                                          //return 1 (pair type match)
        } else {                                               //otherwise
            return 0;                                          //return 0 (no pair type match)
        }
        
    } else { //unrecognized mode
        printf("test_possible_pairs: error - unrecognized mode\n. aborting...");
        abort();
    }
    
    return -1; //this is unreachable
}

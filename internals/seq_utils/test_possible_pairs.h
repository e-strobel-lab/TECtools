//
//  test_possible_pairs.h
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#ifndef test_possible_pairs_h
#define test_possible_pairs_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define TEST_MISMATCH 0
#define TEST_BASEPAIR 1
#define TEST_PAIRTYPE 2

#define NO_TYPE_TEST -1
#define TEST_NO_PAIR  0
#define TEST_GT       1
#define TEST_AT       2
#define TEST_GT_AT    3
#define TEST_GC       4
#define TEST_GT_GC    5
#define TEST_AT_GC    6
#define TEST_GT_AT_GC 7

int test_possible_pairs(char bs1, char bs2, int mode, int type_test);

#endif /* test_possible_pairs_h */

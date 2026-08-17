//
//  mutual_inf.h
//  
//
//  Created by Eric Strobel on 7/31/26.
//

#ifndef mutual_inf_h
#define mutual_inf_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include "data_dist.h"


#define B_A 0
#define B_T 1
#define B_G 2
#define B_C 3
#define B_POSS 4

#define P_AA 0
#define P_AT 1
#define P_AG 2
#define P_AC 3
#define P_TA 4
#define P_TT 5
#define P_TG 6
#define P_TC 7
#define P_GA 8
#define P_GT 9
#define P_GG 10
#define P_GC 11
#define P_CA 12
#define P_CT 13
#define P_CG 14
#define P_CC 15
#define P_POSS 16

double mutual_inf(char * str1, char * str2);
void get_marginal_probs(double * P, char * STR, int len);
void get_joint_probs(double * J, char * STR1, char * STR2, int len);
double calculate_MI(double * J, double * P1, double * P2);

#endif /* mutual_inf_h */

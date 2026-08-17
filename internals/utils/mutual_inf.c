//
//  mutual_inf.c
//  
//
//  Created by Eric Strobel on 7/31/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include "data_dist.h"

#include "mutual_inf.h"

double mutual_inf(char * str1, char * str2)
{
    int i = 0;
    int len = 0;
    
    len = strlen(str1);
    if (strlen(str2) != len) {
        printf("mutual_inf: error - input string lengths are not equal. aborting...\n");
    }
    
    char * STR1 = NULL;
    char * STR2 = NULL;
    
    if ((STR1 = malloc(len * sizeof(*str1))) == NULL) {
        printf("mutual_inf: error - memory allocation failed. aborting...\n");
        abort();
    }
    
    if ((STR2 = malloc(len * sizeof(*str2))) == NULL) {
        printf("mutual_inf: error - memory allocation failed. aborting...\n");
        abort();
    }
    
    for (i = 0; i < len; i++) {
        
        if (str1[i] == 'U' || str1[i] == 'u') {
            STR1[i] = 'T';
        } else {
            STR1[i] = toupper(str1[i]);
        }
        
        if (str2[i] == 'U' || str2[i] == 'u') {
            STR2[i] = 'T';
        } else {
            STR2[i] = toupper(str2[i]);
        }
        
    }
    STR1[i] = '\0';
    STR2[i] = '\0';
    
    
    
    double P1[B_POSS] = {0};
    double P2[B_POSS] = {0};
    double J[P_POSS]  = {0};
    
    double MI = 0;
    
    get_marginal_probs(&P1[0], STR1, len);
    get_marginal_probs(&P2[0], STR2, len);
    
    get_joint_probs(&J[0], STR1, STR2, len);
    
    return calculate_MI(&J[0], &P1[0], &P2[0]);
}

void get_marginal_probs(double * P, char * STR, int len)
{
    int i = 0;
    
    double c[B_POSS] = {0};
    
    for (i = 0; i < len; i++) {
        switch (STR[i]) {
            case 'A': c[B_A]++; break;
            case 'a': c[B_A]++; break;
            case 'T': c[B_T]++; break;
            case 't': c[B_T]++; break;
            case 'G': c[B_G]++; break;
            case 'g': c[B_G]++; break;
            case 'C': c[B_C]++; break;
            case 'c': c[B_C]++; break;
                
            default:
                printf("mutual_inf: error - unrecognized base %c. aborting...\n", STR[i]);
                break;
        }
    }
    
    for (i = 0; i < B_POSS; i++) {
        P[i] = c[i]/((double)len);
        printf("%f\n", P[i]);
    }
    printf("\n");
    
    return;
}

void get_joint_probs(double * J, char * STR1, char * STR2, int len)
{
    int i = 0;
    
    double c[P_POSS] = {0};
    
    for (i = 0; i < len; i++) {
        if      (STR1[i] == 'A' && STR2[i] == 'A') {c[P_AA]++;}
        else if (STR1[i] == 'A' && STR2[i] == 'T') {c[P_AT]++;}
        else if (STR1[i] == 'A' && STR2[i] == 'G') {c[P_AG]++;}
        else if (STR1[i] == 'A' && STR2[i] == 'C') {c[P_AC]++;}
        
        else if (STR1[i] == 'T' && STR2[i] == 'A') {c[P_TA]++;}
        else if (STR1[i] == 'T' && STR2[i] == 'T') {c[P_TT]++;}
        else if (STR1[i] == 'T' && STR2[i] == 'G') {c[P_TG]++;}
        else if (STR1[i] == 'T' && STR2[i] == 'C') {c[P_TC]++;}
        
        else if (STR1[i] == 'G' && STR2[i] == 'A') {c[P_GA]++;}
        else if (STR1[i] == 'G' && STR2[i] == 'T') {c[P_GT]++;}
        else if (STR1[i] == 'G' && STR2[i] == 'G') {c[P_GG]++;}
        else if (STR1[i] == 'G' && STR2[i] == 'C') {c[P_GC]++;}
        
        else if (STR1[i] == 'C' && STR2[i] == 'A') {c[P_CA]++;}
        else if (STR1[i] == 'C' && STR2[i] == 'T') {c[P_CT]++;}
        else if (STR1[i] == 'C' && STR2[i] == 'G') {c[P_CG]++;}
        else if (STR1[i] == 'C' && STR2[i] == 'C') {c[P_CC]++;}
        
        else {
            printf("get_joint_probs: error - unrecognized condition. aborting...\n");
            abort();
        }
    }
    
    for (i = 0; i < P_POSS; i++) {
        J[i] = c[i]/len;
        printf("%f\n", J[i]);
    }
    printf("\n");
}

double calculate_MI(double * J, double * P1, double * P2)
{
    int i = 0;
    
    double MI = 0;
    
    for (i = 0; i < P_POSS; i++) {
        
        if (J[i] != 0.0) {
            switch (i) {
                case P_AA: MI += J[P_AA] * log2(J[P_AA]/(P1[B_A] * P2[B_A])); break;
                case P_AT: MI += J[P_AT] * log2(J[P_AT]/(P1[B_A] * P2[B_T])); break;
                case P_AG: MI += J[P_AG] * log2(J[P_AG]/(P1[B_A] * P2[B_G])); break;
                case P_AC: MI += J[P_AC] * log2(J[P_AC]/(P1[B_A] * P2[B_C])); break;
                
                case P_TA: MI += J[P_TA] * log2(J[P_TA]/(P1[B_T] * P2[B_A])); break;
                case P_TT: MI += J[P_TT] * log2(J[P_TT]/(P1[B_T] * P2[B_T])); break;
                case P_TG: MI += J[P_TG] * log2(J[P_TG]/(P1[B_T] * P2[B_G])); break;
                case P_TC: MI += J[P_TC] * log2(J[P_TC]/(P1[B_T] * P2[B_C])); break;
                    
                case P_GA: MI += J[P_GA] * log2(J[P_GA]/(P1[B_G] * P2[B_A])); break;
                case P_GT: MI += J[P_GT] * log2(J[P_GT]/(P1[B_G] * P2[B_T])); break;
                case P_GG: MI += J[P_GG] * log2(J[P_GG]/(P1[B_G] * P2[B_G])); break;
                case P_GC: MI += J[P_GC] * log2(J[P_GC]/(P1[B_G] * P2[B_C])); break;
                    
                case P_CA: MI += J[P_CA] * log2(J[P_CA]/(P1[B_C] * P2[B_A])); break;
                case P_CT: MI += J[P_CT] * log2(J[P_CT]/(P1[B_C] * P2[B_T])); break;
                case P_CG: MI += J[P_CG] * log2(J[P_CG]/(P1[B_C] * P2[B_G])); break;
                case P_CC: MI += J[P_CC] * log2(J[P_CC]/(P1[B_C] * P2[B_C])); break;
                    
                default:
                    printf("calculate_MI: error - unrecognized condition. aborting...\n");
                    break;
            }
            
            if (J[i] != 0.0) {
                switch (i) {
                    case P_AA:
                        printf("AA\n");
                        printf("num = %f\n", J[P_AA]);
                        printf("1ba = %f\n", P1[B_A]);
                        printf("2ba = %f\n", P2[B_A]);
                        printf("den = %f\n", P1[B_A] * P2[B_A]);
                        printf("quo = %f\n", (J[P_AA]/(P1[B_A] * P2[B_A])));
                        printf("log = %f\n", log2(J[P_AA]/(P1[B_A] * P2[B_A])));
                        printf("PMI = %f\n", J[P_AA] * log2(J[P_AA]/(P1[B_A] * P2[B_A])));
                        break;
                    case P_AT:
                        printf("AT\n");
                        printf("num = %f\n", J[P_AT]);
                        printf("1ba = %f\n", P1[B_A]);
                        printf("2bt = %f\n", P2[B_T]);
                        printf("den = %f\n", P1[B_A] * P2[B_T]);
                        printf("quo = %f\n", (J[P_AT]/(P1[B_A] * P2[B_T])));
                        printf("log = %f\n", log2(J[P_AT]/(P1[B_A] * P2[B_T])));
                        printf("PMI = %f\n", J[P_AT] * log2(J[P_AT]/(P1[B_A] * P2[B_T])));
                        break;
                    case P_AG: printf("%f\n", log2(J[P_AG]/(P1[B_A] * P2[B_G]))); break;
                    case P_AC: printf("%f\n", log2(J[P_AC]/(P1[B_A] * P2[B_C]))); break;
                    
                    case P_TA:
                        printf("TA\n");
                        printf("num = %f\n", J[P_TA]);
                        printf("1bt = %f\n", P1[B_T]);
                        printf("2ba = %f\n", P2[B_A]);
                        printf("den = %f\n", P1[B_T] * P2[B_A]);
                        printf("quo = %f\n", (J[P_TA]/(P1[B_T] * P2[B_A])));
                        printf("log = %f\n", log2(J[P_TA]/(P1[B_T] * P2[B_A])));
                        printf("PMI = %f\n", J[P_TA] * log2(J[P_TA]/(P1[B_T] * P2[B_A])));
                        break;
                    case P_TT: printf("%f\n", log2(J[P_TT]/(P1[B_T] * P2[B_T]))); break;
                    case P_TG: printf("%f\n", log2(J[P_TG]/(P1[B_T] * P2[B_G]))); break;
                    case P_TC: printf("%f\n", log2(J[P_TC]/(P1[B_T] * P2[B_C]))); break;
                        
                    case P_GA: printf("%f\n", log2(J[P_GA]/(P1[B_G] * P2[B_A]))); break;
                    case P_GT: printf("%f\n", log2(J[P_GT]/(P1[B_G] * P2[B_T]))); break;
                    case P_GG: printf("%f\n", log2(J[P_GG]/(P1[B_G] * P2[B_G]))); break;
                    case P_GC:
                        printf("GC\n");
                        printf("num = %f\n", J[P_GC]);
                        printf("1bg = %f\n", P1[B_G]);
                        printf("2bc = %f\n", P2[B_C]);
                        printf("den = %f\n", P1[B_G] * P2[B_C]);
                        printf("quo = %f\n", (J[P_GC]/(P1[B_G] * P2[B_C])));
                        printf("log = %f\n", log2(J[P_GC]/(P1[B_G] * P2[B_C])));
                        printf("PMI = %f\n", J[P_GC] * log2(J[P_GC]/(P1[B_G] * P2[B_C])));
                        break;
                        
                    case P_CA: printf("%f\n", log2(J[P_CA]/(P1[B_C] * P2[B_A]))); break;
                    case P_CT: printf("%f\n", log2(J[P_CT]/(P1[B_C] * P2[B_T]))); break;
                    case P_CG:
                        printf("num = %f\n", J[P_CG]);
                        printf("1bc = %f\n", P1[B_C]);
                        printf("2bg = %f\n", P2[B_G]);
                        printf("den = %f\n", P1[B_C] * P2[B_G]);
                        printf("quo = %f\n", (J[P_CG]/(P1[B_C] * P2[B_G])));
                        printf("log = %f\n", log2(J[P_CG]/(P1[B_C] * P2[B_G])));
                        printf("PMI = %f\n", J[P_CG] * log2(J[P_CG]/(P1[B_C] * P2[B_G])));
                        break;
                    case P_CC: printf("%f\n", log2(J[P_CC]/(P1[B_C] * P2[B_C]))); break;
                        
                    default:
                        printf("calculate_MI: error - unrecognized condition. aborting...\n");
                        break;
                }
            printf("MI  = %f\n\n", MI);
            }
        }
    }
 
    printf("MI = %f\n", MI);
    return MI;
}

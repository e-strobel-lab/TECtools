//
//  make_pair_type_table.c
//  
//
//  Created by Eric Strobel on 6/8/23.
//

#include <stdio.h>
#include <stdlib.h>

#include "./ispair.h"

int main(void)
{
    char base_list[16] = "ATGCRYMKSWBDHVN";
    
    char base_table[16][5] = {
        {"A"},      //A
        {"T"},      //T
        {"G"},      //G
        {"C"},      //C
        {"AG"},     //R
        {"TC"},     //Y
        {"AC"},     //M
        {"GT"},     //K
        {"GC"},     //S
        {"AT"},     //W
        {"TGC"},    //B
        {"AGT"},    //D
        {"ATC"},    //H
        {"AGC"},    //V
        {"ATGC"},   //N
        {0}
    };
    
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    
    int pr_val = 0;
    
    uint8_t PR_bit = 0;
    uint8_t GT_bit = 1;
    uint8_t AT_bit = 2;
    uint8_t GC_bit = 4;
    
    printf("\n");
    
    printf("//      A  T  G  C  R  Y  M  K  S  W  B  D  H  V  N\n");
    
    for (i = 0; base_list[i]; i++) {
    
        printf("/*%c*/ { ", base_list[i]);
        
        for (j = 0; base_table[j][0]; j++) {
            
            pr_val = 0;
            PR_bit = 0;

            for (k = 0; base_table[i][k]; k++) {
                for (l = 0; base_table[j][l]; l++) {
                    pr_val = ispair(base_table[i][k], base_table[j][l]);
                    if (pr_val == GC_PAIR) {
                        PR_bit |= GC_bit;
                    } else if (pr_val == AT_PAIR) {
                        PR_bit |= AT_bit;
                    } else if (pr_val == GU_PAIR) {
                        PR_bit |= GT_bit;
                    }
                }
            }
                
            printf("%d", PR_bit);
            
            if (base_list[j+1]) {
                printf(", ");
            }

        }
        
        if (base_list[i+1]) {
            printf(" },\n");
        } else {
            printf(" }\n");
        }
    }
    
    return 1;
}

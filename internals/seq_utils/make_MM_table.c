//
//  make_pair_table.c
//  
//
//  Created by Eric Strobel on 7/7/22.
//

#include <stdio.h>
#include <stdlib.h>

#include "make_pair_table.h"
#include "./ispair.h"

#define WEAK 0
#define STRONG 1

int main(void)
{
    char base_list[16] = "ATGCRYMKSWBDHVN";
    
    char base_table[16][5] = {
        {"A"},		//A
        {"T"},		//T
        {"G"},		//G
        {"C"},		//C
        {"AG"},		//R
        {"TC"},		//Y
        {"AC"},		//M
        {"GT"},		//K
        {"GC"},		//S
        {"AT"},		//W
        {"TGC"},	//B
        {"AGT"},	//D
        {"ACT"},	//H
        {"ACG"},	//V
        {"ATGC"},	//N
        {0}
    };
    
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int pr_val = 0;
    int MM_found = 0;
    int pairs[2] = {0};
    
    printf("\n");
    
    printf("//      A  T  G  C  R  Y  M  K  S  W  B  D  H  V  N\n");
    
    for (i = 0; base_list[i]; i++) {
    
        printf("/*%c*/ { ", base_list[i]);
        
        for (j = 0; base_table[j][0]; j++) {
            
            MM_found = 0;

            for (k = 0; base_table[i][k]; k++) {
                
                for (l = 0; base_table[j][l]; l++) {
                    if (!ispair(base_table[i][k], base_table[j][l])) {
                        MM_found = 1;
                    }

                }
            }
            if (MM_found) {
                printf("1");
            } else {
                printf("0");
            }
            
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



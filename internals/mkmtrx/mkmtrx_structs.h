//
//  mkmtrx_structs.h
//  
//
//  Created by Eric Strobel on 10/12/22.
//

#ifndef mkmtrx_structs_h
#define mkmtrx_structs_h

#include <stdio.h>

#include "../global/global_defs.h"

typedef struct alignment_stats { //stores alignment stats from bowtie2 output
    int tot_C[2];   //calculated total reads (per read, not pairs)
    int algnd[2];   //calculated aligned reads (per read, not pairs)
    int tot_P[2];   //total reads analyzed
    int prd[2];       //paired reads
    int prd_c0[2];      //aligned concordantly 0 times
    int prd_c1[2];      //aligned concordantly 1 time
    int prd_cM[2];      //aligned concordantly more than 1 time
    int prd_d1[2];      //aligned discordantly 1 time
    int prd_0[2];       //aligned 0 times concorantly or discordantly
    int prd_0m[2];        //mates that make up these pairs
    int prd_0m0[2];         //aligned 0 times
    int prd_0m1[2];         //aligned 1 time
    int prd_0mM[2];         //aligned more than 1 time
    int unp[2];       //unpaired reads
    int unp_0[2];       //aligned 0 times
    int unp_1[2];       //aligned 1 time
    int unp_M[2];       //aligned more than 1 time
    double rate[2];    //overall alignment rate (parsed)
    double calc[2];    //overall alignment rate (calculated)
} alignment_stats;


#endif /* mkmtrx_structs_h */

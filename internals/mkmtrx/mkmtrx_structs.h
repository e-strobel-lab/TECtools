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

typedef struct rdat_metadata {   //stores metadata for rdat generation
    char fn[MAX_MD_STRING+1];    //output rdat file name
    char nm[MAX_MD_STRING+1];    //sample name
    char str[MAX_MD_STRING+1];   //structure
    int  offset;                 //sequence position offset
    int  smooth;                 //flag to indicate if neighboring transcript smoothing was applied
    char xtype[MAX_MD_STRING+1]; //experiment type
    char xmeth[MAX_MD_STRING+1]; //experiment method
    char probe[MAX_MD_STRING+1]; //chemical probe
    char tris[MAX_MD_STRING+1];  //tris concentration
    char KCl[MAX_MD_STRING+1];   //KCl concentration
    char EDTA[MAX_MD_STRING+1];  //EDTA concentration
    char DTT[MAX_MD_STRING+1];   //DTT concentration
    char MgCl2[MAX_MD_STRING+1]; //MgCl2 concentration
    char NTPs[MAX_MD_STRING+1];  //NTP concentration
    char BSA[MAX_MD_STRING+1];   //BSA concentration
    char other_chem[MAX_OTHER_CHEM][MAX_MD_STRING+1]; //array of other chemicals in reaction
    char temp[MAX_MD_STRING+1];  //temperature
    char comment[MAX_COMMENTS][MAX_LINE]; //array to store user-supplied comments
    int oc_cnt;                  //count of other chemical entries provided
    int cmnt_cnt;                //count of user-supplied comments provided
} rdat_metadata;


#endif /* mkmtrx_structs_h */

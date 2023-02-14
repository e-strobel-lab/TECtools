//
//  cotrans_mtrx.h
//  
//
//  Created by Eric Strobel on 10/11/22.
//

#ifndef cotrans_mtrx_h
#define cotrans_mtrx_h

#include <stdio.h>

#define MAX_NAME 128     //max array size for cotrans_matrix name
#define MAX_ROW 512      //max number of rows in cotrans_matrix
#define MAX_COL 512      //max number of columns in cotrans_matrix
#define MIN 0            //array index for transcript length window minimum
#define MAX 1            //array index for transcript length window maximum

#define UNT 0 //untreated read channel index
#define MOD 1 //modified read channel index
#define RCT 2 //reactivity index, used for val_cnt

/* structure declarations */
typedef struct cotrans_matrix {     //struct for storing cotranscriptional RNA structure probing matrices
    char nm[MAX_NAME];              //matrix name
    char fn[MAX_NAME];              //filename
    char sq[MAX_ROW];               //sequence
    char Nrchd[MAX_ROW];            //flag that transcript length is enriched
    char * vals[MAX_ROW][MAX_COL];  //multidimensional array of pointers for storing matrix value strings
    char * unt[MAX_ROW][MAX_COL];   //multidimensional array of pointers for storing untreated read count strings
    char * mod[MAX_ROW][MAX_COL];   //multidimensional array of pointers for stroing modified read count strings
    int row_cnt;                    //number of rows in matrix
    int col_cnt;                    //number of columns in matrix
    int tl[2];                      //transcript length min and max
    int nt[2];                      //nucleotide min and max
    int last_iStl;                  //position of the last internal stall site
    int frst_tStl;                  //position of the first terminal stall site
    int last_tStl;                  //position of the last terminal stall site
    int alias;                      //flag that alias was applied to matrix name
    struct cotrans_matrix * nxt;    //pointer to next matrix in linked list
} cotrans_matrix;

#endif /* cotrans_mtrx_h */

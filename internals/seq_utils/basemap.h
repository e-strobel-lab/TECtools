//
//  basemap.h
//  
//
//  Created by Eric Strobel on 5/30/23.
//

#ifndef basemap_h
#define basemap_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/test_possible_pairs.h"
#include "../seq_utils/seq2bin_hash.h"

#define MAXLEN 500             //maximum input sequence length
#define MAX_SEC_STRUCT 8       //maximum number of secondary structure lines
#define MAX_VBASE_NAME_FIELD 8 //maximum variable base name field length

#define NO_MASK 0    //indicates that the 2nd mate of a pair should NOT be masked in get_vbases
#define MASK_PAIRS 1 //indicates that the 2nd mate of pair should be masked in get_vbases

//variable base table indices
#define N 0
#define R 1
#define Y 2
#define M 3
#define K 4
#define S 5
#define W 6
#define B 7
#define D 8
#define H 9
#define V 10

//base type definitions
#define VARIABLE 0
#define DELETION 1
#define CONSTANT 2
#define SPACER 3

/* structure declarations */
typedef struct wt_source {   //structure to store fasta entries
    char * nm;           //sequence name
    char * sq;           //sequence with multiple sequence alignment characters preserved
    char ** lkp;         //pointers to wt sequence indexed by position number
    char ** lkp0;        //pos 0 in lkp array
    char ** lkp_start;   
    int * pos;           //position number of base
    int ex;              //flag to exclude variant
    int start_pos;       //position of first nucleotide in sequence
} wt_source;

//template expansion structures
static struct varbase {  //stores information for variable positions
    int crnt;            //current index
    char not;            //IUPAC notation
    char alt[5];         //alternate bases at variable position
} vb_tbl[11] = {
    0, 'N', "ATGC", //0
    0, 'R', "AG",   //1
    0, 'Y', "TC",   //2
    0, 'M', "AC",   //3
    0, 'K', "GT",   //4
    0, 'S', "GC",   //5
    0, 'W', "AT",   //6
    0, 'B', "TGC",  //7
    0, 'D', "ATG",  //8
    0, 'H', "ATC",  //9
    0, 'V', "ACG"   //10
    
};

typedef struct insertion {
    int pos;        //position of last non-insertion base that precedes the insertion
    int ins_contig; //position of the insertion in an insertion contig
    char seq;       //identity of the inserted base
} insertion;

typedef struct deletion {
    int pos;        //sequence position of the deletion (not index)
    char seq;       //identity of the deleted base
} deletion;

typedef struct basemap { //map of constant and variable positions in reference sequence
    char * nm;                         //variant template name
    char * rS;                         //variant template reference sequence
    char * nts;                        //sequence, * denotes variable base, P denotes WC pair
    char ** lkp;                       //pointers to reference sequence indexed by position number
    char ** lkp0;                      //pos 0 in lkp array
    char ** lkp_start;                 //pointer to start of reference sequence
    int  * pos;                        //position number of basemap nucleotides
    int  * typ;                        //nucleotide type
    char * rP[MAX_SEC_STRUCT];         //variant template reference pairs
    int  * prs[MAX_SEC_STRUCT];        //indices of pairing partners
    struct varbase ** p_vb;            //pointer to vb_tbl table entry
    struct insertion * v_ins;          //array of variable insertions
    struct insertion * c_ins;          //array of constant insertions
    struct deletion * dels;            //array of deletions
    int * vb_ix;                       //array of all variable base indices
    int * vi_ix;                       //array of variable insertion indices
    int * ci_ix;                       //array of constant insertion indices
    int * d_ix;                        //array of deletion indices
    int vb_cnt;                        //number of variable bases
    int vi_cnt;                        //number of variable insertions
    int ci_cnt;                        //number of constant insertions
    int d_cnt;                         //number of deletions
    int rP_cnt;                        //number of reference pairs
    wt_source * wt;                    //pointer to source sequence fasta
    int cnt[3];                        //number of variants encoded by the variant template
} basemap;

/* set_wt_seq: initialize wt sequence structure */
void set_wt_seq(wt_source * wt, char * nm, char * sq);

/* init_vtmp_mem: initialize memory in basemap struct for storing variant template */
int init_vtmp_mem(basemap * bmap, int nm_len, int sq_len);

/* set_lookup: set wt sequence positional lookup pointers */
void set_lookup(char *** lkp0, char ***lkp_start, char ** lkp, int * wt_pos, char * sq);

/* get_vbases: parses input variant template string and points variable bases to vbase table entry */
void get_vbases(struct basemap *bmap, int mode);

/* set_wt_seq: set basemap structure sequence values */
void set_basemap_seq(basemap * bmap, char * p_nm, char * p_vb, wt_source * wt);

/* add_pairs_to_basemap: add pairing constraint to basemap struct */
void add_pairs_to_basemap(basemap * bmap, char * p_pr, char * p_id);

/* set_basemap_pairs: set pair parameters in the basemap struct. */
void set_basemap_pairs(basemap * bmap, int mode);

/* mk_vbase_nm: construct variable base name */
void mk_vbase_nm(basemap * bmap, int i, char *vb_nm, int array_size, char vbase);

/* print_basemap: print basemap sequence/structure information */
void print_basemap(basemap * bmap);

/* print_vbases: print table of variable base information */
void print_vbases(basemap * bmap);

#endif /* basemap_h */

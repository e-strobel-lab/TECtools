//
//  filter_by_characteristics_structs.h
//  
//
//  Created by Eric Strobel on 11/20/25.
//

#ifndef filter_by_characteristics_structs_h
#define filter_by_characteristics_structs_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "./set_seq_attributes/parse_descriptor_file.h"
#include "./set_seq_attributes/parse_ct_file.h"

/* sequence_attributes: structure for storing input sequences and their respective attributes*/
typedef struct sequence_attributes {
    char * nm;
    char * sq0;
    char * sq1;
    descriptor * des;
    char * tecd;
    void ** att;
} sequence_attributes;

/* structure_characteristics */
typedef struct structProps {
    sequence_attributes * sq_att;
    void * att;
    int att_typ;
    char * db;
    char * db_an;
    double dG;
    struct structProps * nxt;
} structProps;

/* nucleotide_identity: nucleotide identity attribute structure */
typedef struct nucleotide_identity {
    descriptor * des;
    sequence_attributes * sq_att;
    char * sq;
} nucleotide_identity;

/* proximal_deltaG: proximal deltaG attribute structure */
typedef struct proximal_deltaG {
    descriptor * des;
    sequence_attributes * sq_att;
    char * sq;
    int sp_cnt;
    structProps sp;
} proximal_deltaG;

/* distal_deltaG: distal deltaG attribute structure */
typedef struct distal_deltaG {
    descriptor * des;
    sequence_attributes * sq_att;
    char * up;
    char * dn;
    int sp_cnt;
    structProps sp;
} distal_deltaG;

/* subsequence_length: subsequence length attribute structure */
typedef struct subsequence_length {
    descriptor * des;
    sequence_attributes * sq_att;
    char * sq;
    int len;
} subsequence_length;

/* descriptor_bank: attribute structure bank for all attribute types*/
typedef struct descriptor_bank {
    int cnt[DSCRPTR_TYPES];
    nucleotide_identity * nuc_id;
    proximal_deltaG * prx_dG;
    distal_deltaG * dst_dG;
    subsequence_length * ss_len;
} descriptor_bank;

#endif /* filter_by_characteristics_structs_h */

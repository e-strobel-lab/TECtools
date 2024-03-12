//
//  store_SM2_profile.h
//  
//
//  Created by Eric Strobel on 2/10/24.
//

#ifndef store_SM2_profile_h
#define store_SM2_profile_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../../global/global_defs.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "../../utils/gen_utils.h"
#include "../../seq_utils/isRNAbase.h"

#define PRFL_CLMNS 29   //number of columns in shapemapper 2 profile files

#define NUCLEOTIDE 0
#define SEQUENCE 1
#define MOD_MUTATIONS 2
#define MOD_READ_DEPTH 3
#define MOD_EFF_DEPTH 4
#define MOD_RATE 5
#define MOD_OFF_TARGET_DEPTH 6
#define MOD_LOW_MAPQ_DEPTH 7
#define MOD_MAPPED_DEPTH 8
#define UNT_MUTATIONS 9
#define UNT_READ_DEPTH 10
#define UNT_EFF_DEPTH 11
#define UNT_RATE 12
#define UNT_OFF_TARGET_DEPTH 13
#define UNT_LOW_MAPQ_DEPTH 14
#define UNT_MAPPED_DEPTH 15
#define DEN_MUTATIONS 16
#define DEN_READ_DEPTH 17
#define DEN_EFF_DEPTH 18
#define DEN_RATE 19
#define DEN_OFF_TARGET_DEPTH 20
#define DEN_LOW_MAPQ_DEPTH 21
#define DEN_MAPPED_DEPTH 22
#define REACTIVITY_PROFILE 23
#define STD_ERR 24
#define HQ_PROFILE 25
#define HQ_STDERR 26
#define NORM_PROFILE 27
#define NORM_STDERR 28

typedef struct channel_tracker {
    int mod;
    int unt;
    int den;
} channel_tracker;

typedef struct SM2_profile {
    int * nucleotide;
    char * sequence;
    int * mod_mutations;
    int * mod_read_depth;
    int * mod_eff_depth;
    double * mod_rate;
    int * mod_off_target_depth;
    int * mod_low_mapq_depth;
    int * mod_mapped_depth;
    int * unt_mutations;
    int * unt_read_depth;
    int * unt_eff_depth;
    double * unt_rate;
    int * unt_off_target_depth;
    int * unt_low_mapq_depth;
    int * unt_mapped_depth;
    int * den_mutations;
    int * den_read_depth;
    int * den_eff_depth;
    double * den_rate;
    int * den_off_target_depth;
    int * den_low_mapq_depth;
    int * den_mapped_depth;
    double * reactivity_profile;
    double * std_err;
    double * hq_profile;
    double * hq_stderr;
    double * norm_profile;
    double * norm_stderr;
    double * recalc_norm_profile;
    double * dataset_norm_profile;
    int trgt_start;
    int trg_nt_cnt;
    int tot_nt_cnt;
    channel_tracker chnls;
} SM2_profile;

/* store_SM2_profile: store shapemapper 2 profile in SM2_profile struct */
void store_SM2_profile(struct SM2_profile * prf, char * filepath);

/* get_line_local: get line from file, place into array, remove trailing newline, and return
 line length if successful. local version that allows files to end on non-newline characters */
int get_line_local(char *line, FILE *ifp);

/* validate_header: confirm that column headers match expectations */
void validate_header(char * hdr);

/* allocate_SM2_profile_memory: */
void allocate_SM2_profile_memory(struct SM2_profile * prf, int data_lines);

/* check_seq_str: check that sequence string is an RNA base and only a single character */
void check_seq_str(char * str);

#endif /* store_SM2_profile_h */

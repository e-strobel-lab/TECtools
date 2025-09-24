//
//  get_profiles.h
//  
//
//  Created by Eric Strobel on 10/11/22.
//

#ifndef get_profiles_h
#define get_profiles_h

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "./cotrans_mtrx.h"
#include "./mkmtrx_defs.h"
#include "./mkmtrx_structs.h"
#include "./find_nxt_dir_entry.h"
#include "./parse_log.h"

#include "../process_TECprobe_profiles/UNV/parse_sample_name.h"

/* get_profiles: read input directory to find shapemapper output files.
 then call parse_profile to store reactivity data in cotrans matrix struct*/
int get_profiles(char * prnt_dir, cotrans_matrix * mtrx, alignment_stats * algn, char * smpl_nm, int rct_typ, int preprocessed, int test_SM2_data);

/* parse_profile: shapemapper output file to validate file integrity
 and copy reactivity values into cotrans matrix struct */
int parse_profile(char * smo_loc, int len, cotrans_matrix * mtrx, int rct_typ);

/* vldt_SMO_hdr: validate shapemapper output file header by checking first column,
 columns that will be used in parse_profiles, and the number of columns */
int vldt_SMO_hdr(char * line);

/* check_contiguity: check that reactivity matrix is
 contiguous from minimum to maximum transcript lengths */
int check_contiguity(cotrans_matrix * mtrx);

int set_enriched_lengths(cotrans_matrix *mtrx, int incld_up2, int excld_trmnl);

#endif /* get_profiles_h */

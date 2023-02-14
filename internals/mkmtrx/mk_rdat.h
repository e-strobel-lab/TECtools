//
//  mk_rdat.h
//  
//
//  Created by Eric Strobel on 2/9/23.
//

#ifndef mk_rdat_h
#define mk_rdat_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "./cotrans_mtrx.h"
#include "./mkmtrx_defs.h"
#include "./mkmtrx_structs.h"

/* mk_rdat: manages rdata config parsing/checking and
 construction of output rdat file*/
int mk_rdat(FILE * fp_config, cotrans_matrix * mtrx, int mode);

/* parde_rdat_config: read rdat config file and store values in rdat_metadata struct */
int parse_rdat_config(FILE * fp_config, rdat_metadata * rdat_meta);

/* check_rdat_config: check that required rdat metadata has been supplied
 and is consistent with the shapemapper2 output */
int check_rdat_config(rdat_metadata * rdat_meta, cotrans_matrix * mtrx);

/* print_rdat: print rdat file */
int print_rdat(rdat_metadata * rdat_meta, cotrans_matrix * mtrx, int mode);

/* set_TF_value: set true or false value as 1 and 0, respectively */
void set_TF_value(char * tf, char * setting_name, int * config_val);

#endif /* mk_rdat_h */

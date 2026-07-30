//
//  run_RNAStructure.h
//  
//
//  Created by Eric Strobel on 2/27/26.
//

#ifndef run_RNAStructure_h
#define run_RNAStructure_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../seq_utils/mk_fasta.h"

#include "./parse_descriptor_file.h"

/* run_RNAStructure_Fold: predict RNA structures using the RNAStructure Fold algorithm */
int run_RNAStructure_Fold(char * path2RNAStructure, char * nm, char * in_dir, char * out_dir, char * fa_sffx, char * path2ct, int ct_path_maxlen, int ptyp);

#endif /* run_RNAStructure_h */

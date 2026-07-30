//
//  set_attributes.h
//  
//
//  Created by Eric Strobel on 12/1/25.
//

#ifndef set_attributes_h
#define set_attributes_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"
#include "../../seq_utils/mk_fasta.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"

#include "./parse_descriptor_file.h"
#include "./get_msa_subseq.h"
#include "./run_RNAStructure.h"
#include "./parse_ct_file.h"

/* set_attributes: manages sequence attribute setting */
void set_attributes(sequence_attributes * sq_att, descriptor * des, int seq_cnt, int des_cnt, char * path2RNAstructure);

/* set_nuc_id: set nucleotide_identity structure values */
void set_nuc_id(nucleotide_identity * nuc_id, descriptor * des, sequence_attributes * sq_att);

/* set_prx_dG: set proximal_deltaG structure values */
void set_prx_dG(proximal_deltaG * prx_dG, descriptor * des, sequence_attributes * sq_att, char * path2RNAstructure);

/* set_dst_dG: set distal_deltaG structure values */
void set_dst_dG(distal_deltaG * dst_dG, descriptor * des, sequence_attributes * sq_att, char * path2RNAstructure);

/* set_ss_len: set subsequence_length structure values */
void set_ss_len(subsequence_length * ss_len, descriptor * des, sequence_attributes * sq_att);

/* set_structProps: set structProps values */
void set_structProps(sequence_attributes * sq_att, void * att, int att_typ, structProps * sp, con_table * ct, int ct_cnt);

#endif /* set_attributes_h */

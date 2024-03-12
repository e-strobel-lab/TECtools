//
//  VL_input_validation.h
//  
//
//  Created by Eric Strobel on 2/22/24.
//

#ifndef VL_input_validation_h
#define VL_input_validation_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"

#include "../global/store_SM2_profile.h"

/* validate_VL_an_dir_contiguity: check that all transcript length directories
   contained within a TECprobe-VL analysis directory are contiguous and that
   no out-of-bounds transcript lengths are present. */
void validate_VL_an_dir_contiguity(SM2_analysis_directory * an_dir);

/* validate_VL_an_dir_compatibility: verify that TECprobe-VL analysis
   directories contain the same number and range of transcript lengths */
void validate_VL_an_dir_compatibility(SM2_analysis_directory * an_dir, int dir_count);

/* validate_channel_configuration: verify that profile contains a valid channel configuration */
void validate_channel_configuration(channel_tracker * chnls);

/* validate_channel_compatibility: confirm that all transcript lengths contain a modified
   channel, and that the presence of untreated and denatured channels is uniform */
void validate_channel_compatibility(channel_tracker * crrnt_chnls, channel_tracker * first_chnls);

/* validate_int_start_ix_compatibility: verify that all transcript lengths within an analysis directory have the same start index */
void validate_int_start_ix_compatibility(SM2_analysis_directory * an_dir);

/* validate_ext_start_ix_compatibility: verify that two analysis directories have the same start index */
void validate_ext_start_ix_compatibility(int ix1, int ix2);

/* validate_transcript substrings: verify that each transcript sequence of length n is
   a substring of the transcript with length n+1 */
void validate_transcript_substrings(SM2_analysis_directory * an_dir);

/* validate_trg_rct_cnt: verify that target reactivity count is the same */
void validate_trg_rct_cnt(int cnt1, int cnt2);

#endif /* VL_input_validation_h */

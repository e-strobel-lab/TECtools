//
//  VL_input_validation.c
//  
//
//  Created by Eric Strobel on 2/22/24.
//
#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"

#include "../global/store_SM2_profile.h"

#include "VL_input_validation.h"

/* validate_VL_an_dir_contiguity: check that all transcript length directories
   contained within a TECprobe-VL analysis directory are contiguous and that
   no out-of-bounds transcript lengths are present. */
void validate_VL_an_dir_contiguity(SM2_analysis_directory * an_dir)
{
    int i = 0; //general purpose index
    
    //check the transcript length contiguity of the opened profiles
    for (i = 0; i < MAX_ROW; i++) {
        
        //check that all transript lengths are within the min_tl/max_tl bounds
        if (an_dir->prf[i] != NULL && (i < an_dir->min_tl || i > an_dir->max_tl)) {
            printf("validate_VL_an_dir_contiguity: error - out of bounds transcript length in %s data set. aborting...\n", an_dir->prnt_dir_nm);
            abort();
        }
        
        //check that all transcript lengths within the min_tl/max_tl bounds have an associated profile
        if (an_dir->prf[i] == NULL && (i >= an_dir->min_tl && i <= an_dir->max_tl)) {
            printf("validate_VL_an_dir_contiguity: error - missing transcript length %d  in %s data set. aborting...\n", i, an_dir->prnt_dir_nm);
            abort();
        }
    }
    
    return;
}

/* validate_VL_an_dir_compatibility: verify that TECprobe-VL analysis
   directories contain the same number and range of transcript lengths */
void validate_VL_an_dir_compatibility(SM2_analysis_directory * an_dir, int dir_count)
{
    int i = 0; //general purpose index
    
    if (dir_count < 2) { //check that more than 1 analysis directory was supplied
        return;
        
    } else { //verify compatibility of all input directories
        for (i = 1; i < dir_count; i++) {
            
            //check that the number of profiles opened matches that of the first parent directory
            if (an_dir[i].prfs_opnd != an_dir[0].prfs_opnd) {
                printf("validate_VL_an_dir_compatibility: error - the number of reactivity profiles opened is not the same for all input directories (first = %d, current = %d). aborting...\n", an_dir[0].prfs_opnd, an_dir[i].prfs_opnd);
                abort();
            }
            
            //check that the minimum transcript length matches that of the first parent directory
            if (an_dir[i].min_tl != an_dir[0].min_tl) {
                printf("validate_VL_an_dir_compatibility: error - the minimum transcript length is not the same for all input directories. aborting...\n");
                abort();
            }
            
            //check that the maximum transcript length matches that of the first parent directory
            if (an_dir[i].max_tl != an_dir[0].max_tl) {
                printf("validate_VL_an_dir_compatibility: error - the maximum transcript length is not the same for all input directories. aborting...\n");
                abort();
            }
        }
    }
    
    return;
}

/* validate_channel_configuration: verify that profile contains a valid channel configuration */
void validate_channel_configuration(channel_tracker * chnls)
{
    //check that valid channel configuration was detected
    if (!(chnls->mod && !chnls->unt && !chnls->den) &&
        !(chnls->mod &&  chnls->unt && !chnls->den) &&
        !(chnls->mod &&  chnls->unt &&  chnls->den)) {
        printf("validate_channel_configuration: the following invalid channel configuration was detected:\n\n");
        printf("modified:  %s\n", (chnls->mod) ? "detcted" : "not detected");
        printf("untreated: %s\n", (chnls->unt) ? "detcted" : "not detected");
        printf("denatured: %s\n", (chnls->den) ? "detcted" : "not detected");
        printf("\naborting...\n");
        abort();
    }
}

/* validate_channel_compatibility: confirm that all transcript lengths contain a modified
   channel, and that the presence of untreated and denatured channels is uniform */
void validate_channel_compatibility(channel_tracker * crrnt_chnls, channel_tracker * first_chnls)
{
    if (!crrnt_chnls->mod) {
        printf("validate_channel_compatibility: reactivity profile contains no modified channel reads. aborting...\n");
        abort();
    }
    
    if (crrnt_chnls->unt != first_chnls->unt) {
        printf("validate_channel_compatibility: untreated channel reads were detected in some, but not all, reactivity profiles. aborting...\n");
        abort();
    }
    
    if (crrnt_chnls->den != first_chnls->den) {
        printf("validate_channel_compatibility: denatured channel reads were detected in some, but not all, reactivity profiles. aborting...\n");
        abort();
    }
    
    
}

/* validate_start_index_compatibility: verify that all transcript lengths have the same start index */
void validate_int_start_ix_compatibility(SM2_analysis_directory * an_dir)
{
    int i = 0;
    
    for (i = an_dir->min_tl; i < an_dir->max_tl; i++) {
        if (an_dir->data[i].trgt_start != an_dir->data[an_dir->min_tl].trgt_start) {
            printf("validate_start_index_compatibility: error - start index is not uniform across all transcript lengths. aborting\n");
            abort();
        }
    }
    
    return;
}

/* validate_ext_start_ix_compatibility: verify that two analysis directories have the same start index */
void validate_ext_start_ix_compatibility(int ix1, int ix2)
{
    if (ix1 != ix2) {
        printf("validate_ext_start_ix_compatibility: start index is not the same for all input directories (%d vs. %d). aborting...\n", ix1, ix2);
        abort();
    }
    
    return;
}

/* validate_transcript substring: verify that each transcript sequence of length n is
   a substring of the transcript with length n+1 */
void validate_transcript_substrings(SM2_analysis_directory * an_dir)
{
    int i = 0;
    
    char * p_prev = NULL; //pointer to previous transcript sequence substring

    for (i = an_dir->min_tl+1; i < an_dir->max_tl; i++) {
        
        //check that previous sequence is a substring of the current sequence
        if ((p_prev = strstr(an_dir->data[i].sequence, an_dir->data[i-1].sequence)) == NULL) {
            printf("validate_transcript_substrings: error - sequence of transcript %d is not a substring of transcript %d. aborting...\n", i-1, i);
            abort();
            
        }

        //check that the current and previous sequence string start at the same position
        if ((uint64_t)(&p_prev[0]) != (uint64_t)(&an_dir->data[i].sequence[0])) {
            printf("validate_transcript_substrings: error - sequence of transcripts %d and %d do not start at the same position. aborting...\n", i-1, i);
            abort();
            
        }
        
        //check that the current sequence string is one nt shorter than the previous sequence string
        if (strlen(an_dir->data[i-1].sequence) != (strlen(an_dir->data[i].sequence)-1)) {
            printf("validate_transcript_substrings: error - sequence of transcript %d is not 1 nt shorter than sequence of transcript %d. aborting...\n", i-1, i);
            abort();
        }
    }
    
    return;
}

/* validate_trg_rct_cnt: verify that target reactivity count is the same */
void validate_trg_rct_cnt(int cnt1, int cnt2)
{
    if (cnt1 != cnt2) {
        printf("validate_trg_rct_cnt: target RNA reactivity count is not the same for all input directories. aborting...\n");
        abort();
    }
    
    return;
}

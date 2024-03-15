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
        
        //check that all transcript lengths are within the min_tl/max_tl bounds
        if (an_dir->prf[i] != NULL && (i < an_dir->min_tl || i > an_dir->max_tl)) {
            printf("validate_VL_an_dir_contiguity: error - out of bounds transcript length (%d) in %s data set. aborting...\n", i, an_dir->prnt_dir_nm);
            abort();
        }
        
        //check that all transcript lengths within the min_tl/max_tl bounds have
        //an associated profile. throw warning rather than error, since it is
        //possible for a profile to be absent if no reads were assigned to that
        //transcrip length
        if (an_dir->prf[i] == NULL && (i >= an_dir->min_tl && i <= an_dir->max_tl)) {
            printf("validate_VL_an_dir_contiguity: warning - missing transcript length %d in %s data set.\n", i, an_dir->prnt_dir_nm);
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
            
            //check that the number of out directories opened matches that of the first parent directory
            if (an_dir[i].outs_cnt != an_dir[0].outs_cnt) {
                printf("validate_VL_an_dir_compatibility: warning - the number of shapemapper2 analysis directories opened is not the same for all input directories (first = %d, current = %d). aborting...\n", an_dir[0].outs_cnt, an_dir[i].outs_cnt);
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
        printf("validate_channel_compatibility: sample contains no modified channel reads. aborting...\n");
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

/* validate_trgt_start: verify that all transcript lengths have the same start index */
int validate_trgt_start(SM2_analysis_directory * an_dir)
{
    int i = 0; //general purpose index
    
    int first_trg_start = -1; //first target start index
    
    //check that target start index is the same for every profile
    for (i = an_dir->min_tl; i <= an_dir->max_tl; i++) { //for each transcript length
        
        if (an_dir->opnd[i]) {           //if the analysis directory was opened
            if (first_trg_start == -1) { //if the first target start index was not yet found
                first_trg_start = an_dir->data[i].trgt_start; //set the first target start index
                
            } else { //if the first target start index was already found
                
                //check that the current target start index
                //matches the first target start index
                if (an_dir->data[i].trgt_start != first_trg_start) {
                    printf("validate_trgt_start: error - start index is not uniform across all transcript lengths. aborting...\n");
                    abort();
                }
            }
        }
    }
    
    return first_trg_start;
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
    int i = 0;                 //general purpose index
    int skipped = 0;           //tracks how many transcript lengths were skipped
    
    char * p_last_seq = NULL;  //pointer to the last transcript sequence
    char * p_substring = NULL; //pointer to previous sequence substring in current sequence
    
    p_last_seq = &an_dir->data[an_dir->min_opnd].sequence[0]; //init p_last_seq to min opened sequence
    
    for (i = an_dir->min_opnd+1; i <= an_dir->max_opnd; i++) { //for every length after min opened
        
        if (an_dir->opnd[i]) { //if a profile was opened
            
            //check that previous sequence is a substring of the current sequence
            if ((p_substring = strstr(an_dir->data[i].sequence, p_last_seq)) == NULL) {
                printf("validate_transcript_substrings: error - sequence of transcript %d is not a substring of transcript %d. aborting...\n", i-1-skipped, i);
                abort();
                
            }
            
            //check that the current and previous sequence string start at the same position
            if ((uint64_t)(&p_substring[0]) != (uint64_t)(&an_dir->data[i].sequence[0])) {
                printf("validate_transcript_substrings: error - sequence of transcripts %d and %d do not start at the same position. aborting...\n", i-1-skipped, i);
                abort();
                
            }
                        
            //check that the current sequence string is the expected number
            //of nucleotides shorter than the previous sequence string
            if (strlen(p_last_seq) != (strlen(an_dir->data[i].sequence)-1-skipped)) {
                printf("validate_transcript_substrings: error - sequence of transcript %d is not %d nt shorter than sequence of transcript %d. aborting...\n", i-1-skipped, 1+skipped, i);
                abort();
            }
                        
            p_last_seq = &an_dir->data[i].sequence[0]; //set the last seq pointer to the current seq
            skipped = 0;                               //reset the skipped counter to zero
            
        } else {       //profile was not opened for current transcript length
            skipped++; //increment skipped counter
        }
    }
    
    return;
}

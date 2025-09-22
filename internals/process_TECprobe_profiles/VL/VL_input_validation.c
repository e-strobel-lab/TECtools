//
//  VL_input_validation.c
//  
//
//  Created by Eric Strobel on 2/22/24.
//
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

#include "../../global/global_defs.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"

#include "../global/store_SM2_profile.h"
#include "./read_VL_analysis_directories.h"

#include "VL_input_validation.h"

/* validate_VL_an_dir_contiguity: check that all transcript length directories
   contained within a TECprobe-VL analysis directory are contiguous and that
   no out-of-bounds transcript lengths are present. */
void validate_VL_an_dir_contiguity(SM2_analysis_directory * an_dir)
{
    extern const char empty_SM2out[6];
    
    int i = 0; //general purpose index
    
    int * ix = &an_dir->indx[0]; //set pointer to target indices
    
    //check the transcript length contiguity of the opened profiles
    for (i = 0; ix[i] < an_dir->max_id; i++) {
        
        //check that all transcript lengths are within the min_tl/max_tl bounds
        if (an_dir->loc[ix[i]] != NULL && (ix[i] < an_dir->min_id || ix[i] > an_dir->max_id)) {
            printf("validate_VL_an_dir_contiguity: error - out of bounds transcript length (%d) in %s data set. aborting...\n", ix[i], an_dir->prnt_dir_nm);
            abort();
        }
        
        //check for missing transcript lengths
        if (!strcmp(an_dir->loc[ix[i]], empty_SM2out)) {
            printf("validate_VL_an_dir_contiguity: warning - missing transcript length %d in %s data set.\n", ix[i], an_dir->prnt_dir_nm);
        }
    }
    
    return;
}

/* validate_VL_an_dir_compatibility: verify that TECprobe-VL analysis
   directories contain the same number and range of transcript lengths */
void validate_an_dir_compatibility(SM2_analysis_directory * an_dir, int dir_count)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int * ix1 = NULL; //pointer to target index for analysis directory zero
    int * ix2 = NULL; //pointer to target index for other analysis directories
    
    if (dir_count < 2) { //check that more than 1 analysis directory was supplied
        return;
        
    } else { //verify compatibility of all input directories
        
        ix1 = &an_dir[0].indx[0]; //point ix1 to analysis directory 0 index array
        
        for (i = 1; i < dir_count; i++) {
            
            //check that the number of out directories opened matches that of the first parent directory
            if (an_dir[i].outs_cnt != an_dir[0].outs_cnt) {
                printf("validate_an_dir_compatibility: warning - the number of shapemapper2 analysis directories opened is not the same for all input directories (first = %d, current = %d). aborting...\n", an_dir[0].outs_cnt, an_dir[i].outs_cnt);
                abort();
            }
            
            //check that the maximum target id matches that of the first parent directory
            if (an_dir[i].max_id != an_dir[0].max_id) {
                printf("validate_an_dir_compatibility: error - input analysis directories are incompatible. aborting...\n");
            }
            
            //check that file location/"empty" string are located at the same id indices
            for (j = 0; j <= an_dir[0].max_id; j++) {
                if ((an_dir[i].loc[j] == NULL && an_dir[0].loc[j] != NULL) ||
                    (an_dir[i].loc[j] != NULL && an_dir[0].loc[j] == NULL)) {
                    printf("validate_an_dir_compatibility: error - input analysis directories are incompatible. aborting...;");
                    abort();
                }
            }
            
            //check that target id indexing matches that of the first parent directory
            ix2 = &an_dir[i].indx[0]; //point ix2 to current analysis directory index array
            for (j = 0; (ix1[j] <= an_dir[0].max_id) && (ix2[j] <= an_dir[i].max_id); j++) {
                
                if (ix1[j] != ix2[j]) { //if target ids are not aligned, throw error and abort
                    printf("validate_an_dir_compatibility: error - target ids from input analysis directories 0 and %d are not aligned. aborting...\n", ix2[j]);
                    abort();
                }
            }
            
            if (ix1[j] != INT_MAX || ix2[j] != INT_MAX) { //if either loop did not terminate at INT_MAX sentinel
                printf("validate_an_dir_compatibility: error - target ids from input analysis directories 0 and %d are not aligned. aborting...\n", ix2[j]);
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
    
    int * ix = &an_dir->indx[0]; //set pointer to target indices
    
    //check that target start index is the same for every profile
    for (i = 0; ix[i] <= an_dir->max_id; i++) { //for each target
        
        if (an_dir->loc[ix[i]] != NULL) { //if the analysis directory was opened
            if (first_trg_start == -1) {  //if the first target start index was not yet found
                first_trg_start = an_dir->data[ix[i]].trgt_start; //set the first target start index
                
            } else { //if the first target start index was already found
                
                //check that the current target start index
                //matches the first target start index
                if (an_dir->data[ix[i]].trgt_start != first_trg_start) {
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

/* validate_transcript_substrings: verify that each transcript sequence of length n is
   a substring of the transcript with length n+1 */
void validate_transcript_substrings(SM2_analysis_directory * an_dir)
{
    extern const char empty_SM2out[6];
    
    int i = 0;                 //general purpose index
    int skipped = 0;           //tracks how many transcript lengths were skipped
    
    char * p_last_seq = NULL;  //pointer to the last transcript sequence
    char * p_substring = NULL; //pointer to previous sequence substring in current sequence
    
    int * ix = &an_dir->indx[0]; //set pointer to target indices
    
    //point p_last_seq to the sequence of the shortest target profile
    for (i = 0; ix[i] <= an_dir->max_id && p_last_seq == NULL; i++) {
        if (an_dir->data[ix[i]].sequence[0]) {
            p_last_seq = &an_dir->data[ix[i]].sequence[0];
        }
    }
    
    //check that p_last_seq was set
    if (p_last_seq == NULL) {
        printf("validate_transcript_substrings: error - no profiles contained a sequence. aborting...\n");
        abort();
    }
    
    //i was already incremented above
    for (; ix[i] <= an_dir->max_id; i++) { //for every length after min opened
        
        if (strcmp(an_dir->loc[ix[i]], empty_SM2out)) { //if a profile was opened
            
            //check that previous sequence is a substring of the current sequence
            if ((p_substring = strstr(an_dir->data[ix[i]].sequence, p_last_seq)) == NULL) {
                printf("validate_transcript_substrings: error - sequence of transcript %d is not a substring of transcript %d. aborting...\n", ix[i]-1-skipped, ix[i]);
                abort();
            }
            
            //check that the current and previous sequence string start at the same position
            if ((uint64_t)(&p_substring[0]) != (uint64_t)(&an_dir->data[ix[i]].sequence[0])) {
                printf("validate_transcript_substrings: error - sequence of transcripts %d and %d do not start at the same position. aborting...\n", ix[i]-1-skipped, ix[i]);
                abort();
            }
                  
            //check that the current sequence string is the expected number
            //of nucleotides shorter than the previous sequence string
            if (strlen(p_last_seq) != (strlen(an_dir->data[ix[i]].sequence)-1-skipped)) {
                printf("validate_transcript_substrings: error - sequence of transcript %d is not %d nt shorter than sequence of transcript %d. aborting...\n", ix[i]-1-skipped, 1+skipped, ix[i]);
                abort();
            }
                        
            p_last_seq = &an_dir->data[ix[i]].sequence[0]; //set the last seq pointer to the current seq
            skipped = 0;                                   //reset the skipped counter to zero
            
        } else {       //profile was not opened for current transcript length
            skipped++; //increment skipped counter
        }
    }
    
    return;
}

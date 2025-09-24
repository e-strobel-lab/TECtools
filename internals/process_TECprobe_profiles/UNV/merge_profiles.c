//
//  merge_profiles.c
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "../../global/global_defs.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "../../cotrans_preprocessor/cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor/run_script_gen/UNV/config_struct.h"

#include "../../seq_utils/isRNAbase.h"

#include "../process_TECprobe_profiles_defs.h"
#include "../process_TECprobe_profiles_structs.h"
#include "../UNV/calculate_normalization_factor.h"

#include "merge_profiles.h"

/* merge_profiles: combine data from input reactivity profiles to generate merged reactivity profile files */
int merge_profiles(SM2_analysis_directory * an_dir, int dir_count, SM2_analysis_directory * mrg, int min_depth, double max_bkg)
{
    
    int i = 0;     //general purpose index
    int j = 0;     //general purpose index
    int col = 0;   //column index
    int d = 0;     //directory index
    
    char * p_sq = NULL; //pointer to the nucleotide string that will be copied into the merged profile
    
    //allocate memory for mrg profile data
    if ((mrg->data = calloc(an_dir[0].max_id+1, sizeof(*mrg->data))) == NULL) {
        printf("merge_profiles: error - failed to allocate memory for profile data. aborting...\n");
        abort();
    }
    
    //allocate memory for lookup table
    if ((mrg->indx = calloc(an_dir[0].sd_cnt+1, sizeof(*mrg->indx))) == NULL) {
        printf("merge_profiles: error - failed to allocate memory for index table. aborting...\n");
        abort();
    }
    
    //copy index array to mrg. this includes the
    //sentinel that signals the end of the array
    for (i = 0; i <= an_dir[0].sd_cnt; i++) {
        mrg->indx[i] = an_dir[0].indx[i];
    }
    
    //set variables for merged analysis directory to that of the first input analysis
    //directory. these variables are identical for all input analysis directories
    mrg->chnls.mod = an_dir[0].chnls.mod;     //set modified  channel detection flag
    mrg->chnls.unt = an_dir[0].chnls.unt;     //set untreated channel detection flag
    mrg->chnls.den = an_dir[0].chnls.den;     //set denatured channel detection flag
    mrg->outs_cnt = an_dir[0].outs_cnt;       //set SM2 output directories count
    mrg->trgt_start = an_dir[0].trgt_start;   //set target start index
    mrg->min_id = an_dir[0].min_id;           //set minimum target id
    mrg->max_id = an_dir[0].max_id;           //set maximum target id
    mrg->len[MIN] = an_dir[0].len[MIN];       //set minimum transcript length
    mrg->len[MAX] = an_dir[0].len[MAX];       //set maximum transcript length
    mrg->trg_rct_cnt = an_dir[0].trg_rct_cnt; //set target nucleotide reactivity count
       
    int * ix = &an_dir[0].indx[0]; //set pointer to target indices
    
    for (i = 0; ix[i] <= an_dir[0].max_id; i++) {  //for every target id
        
        //check data compatibility
        for (d = 0; d < dir_count; d++) {
            
            //check total nt count
            if (an_dir[d].data[ix[i]].tot_nt_cnt != an_dir[0].data[ix[i]].tot_nt_cnt) {
                printf("merge_profiles: discordant total nucleotide count for target %d. aborting...\n", ix[i]);
                abort();
            }
            
            //check target nt count
            if (an_dir[d].data[ix[i]].trg_nt_cnt != an_dir[0].data[ix[i]].trg_nt_cnt) {
                printf("merge_profiles: discordant target nucleotide count for target %d. aborting...\n", ix[i]);
                abort();
            }
            
            //check target start
            if (an_dir[d].data[ix[i]].trgt_start != an_dir[0].data[ix[i]].trgt_start) {
                printf("merge_profiles: discordant target start for target %d. aborting...\n", ix[i]);
                abort();
            }
        }
        
        //allocate memory for storing the merged profile of the current target
        allocate_SM2_profile_memory(&mrg->data[ix[i]], an_dir[0].data[ix[i]].tot_nt_cnt);

        mrg->data[ix[i]].trgt_start = an_dir[0].data[ix[i]].trgt_start;
        mrg->data[ix[i]].trg_nt_cnt = an_dir[0].data[ix[i]].trg_nt_cnt;
        mrg->data[ix[i]].tot_nt_cnt = an_dir[0].data[ix[i]].tot_nt_cnt;
        
        //set channel tracker flags
        mrg->data[ix[i]].chnls.mod = an_dir[0].chnls.mod;
        mrg->data[ix[i]].chnls.unt = an_dir[0].chnls.unt;
        mrg->data[ix[i]].chnls.den = an_dir[0].chnls.den;
                
        //check if any of the current target id files contain the nucleotide string
        //(The nt seq of empty SM2 out file placeholders is a string of n/N's)
        for (d = 0, p_sq = NULL; d < dir_count && p_sq == NULL; d++) {
            if (isRNAbase(an_dir[d].data[ix[i]].sequence[0])) {
                p_sq = &an_dir[d].data[ix[i]].sequence[0];
            }
        }
         
        //if all SM2 out files are empty, point to the nucleotide string of the first file
        if (p_sq == NULL) {
            p_sq = &an_dir[0].data[ix[i]].sequence[0];
        }
        
        //TODO: check that tot nucleotide count is the same for all inputs
        for (j = 0; j < an_dir[0].data[ix[i]].tot_nt_cnt; j++) { //for every line of the profile file
                        
            mrg->data[ix[i]].nucleotide[j] = an_dir[0].data[ix[i]].nucleotide[j]; //copy nucleotide number
            mrg->data[ix[i]].sequence[j]   = p_sq[j];                             //copy sequence
                        
            mrg->data[ix[i]].std_err[j] = NAN;     //mask std err
            mrg->data[ix[i]].hq_stderr[j] = NAN;   //mask hq std err
            mrg->data[ix[i]].norm_stderr[j] = NAN; //mask normalized std err
                        
            for (col = 0; col < PRFL_CLMNS; col++) { //for each column
                for (d = 0; d < dir_count; d++) {
                    switch (col) {
                        //variables copied outside of loop
                        case NUCLEOTIDE:         break; //do nothing, copied above loop
                        case SEQUENCE:           break; //do nothing, copied above loop
                            
                        //calculations performed outside of loop
                        case MOD_RATE:           break; //do nothing, calculated below loop
                        case UNT_RATE:           break; //do nothing, calculated below loop
                        case DEN_RATE:           break; //do nothing, calculated below loop
                        case REACTIVITY_PROFILE: break; //do nothing, calculated below loop
                        case HQ_PROFILE:         break; //do nothing, set below loop
                        case NORM_PROFILE:       break; //do nothing, generated outside of function
                        case STD_ERR:            break; //do nothing, set above loop
                        case HQ_STDERR:          break; //do nothing, set above loop
                        case NORM_STDERR:        break; //do nothing, set above loop
                            
                        //variables summed inside of loop
                        case MOD_MUTATIONS:
                            mrg->data[ix[i]].mod_mutations[j] += an_dir[d].data[ix[i]].mod_mutations[j];
                            break;
                            
                        case MOD_READ_DEPTH:
                            mrg->data[ix[i]].mod_read_depth[j] += an_dir[d].data[ix[i]].mod_read_depth[j];
                            break;
                            
                        case MOD_EFF_DEPTH:
                            mrg->data[ix[i]].mod_eff_depth[j] += an_dir[d].data[ix[i]].mod_eff_depth[j];
                            break;
                            
                        case MOD_OFF_TARGET_DEPTH:
                            mrg->data[ix[i]].mod_off_target_depth[j] += an_dir[d].data[ix[i]].mod_off_target_depth[j];
                            break;
                            
                        case MOD_LOW_MAPQ_DEPTH:
                            mrg->data[ix[i]].mod_low_mapq_depth[j] += an_dir[d].data[ix[i]].mod_low_mapq_depth[j];
                            break;
                            
                        case MOD_MAPPED_DEPTH:
                            mrg->data[ix[i]].mod_mapped_depth[j] += an_dir[d].data[ix[i]].mod_mapped_depth[j];
                            break;
                            
                        case UNT_MUTATIONS:
                            mrg->data[ix[i]].unt_mutations[j] += an_dir[d].data[ix[i]].unt_mutations[j];
                            break;
                            
                        case UNT_READ_DEPTH:
                            mrg->data[ix[i]].unt_read_depth[j] += an_dir[d].data[ix[i]].unt_read_depth[j];
                            break;
                            
                        case UNT_EFF_DEPTH:
                            mrg->data[ix[i]].unt_eff_depth[j] += an_dir[d].data[ix[i]].unt_eff_depth[j];
                            break;
                            
                        case UNT_OFF_TARGET_DEPTH:
                            mrg->data[ix[i]].unt_off_target_depth[j] += an_dir[d].data[ix[i]].unt_off_target_depth[j];
                            break;
                            
                        case UNT_LOW_MAPQ_DEPTH:
                            mrg->data[ix[i]].unt_low_mapq_depth[j] += an_dir[d].data[ix[i]].unt_low_mapq_depth[j];
                            break;
                            
                        case UNT_MAPPED_DEPTH:
                            mrg->data[ix[i]].unt_mapped_depth[j] += an_dir[d].data[ix[i]].unt_mapped_depth[j];
                            break;
                            
                        case DEN_MUTATIONS:
                            mrg->data[ix[i]].den_mutations[j] += an_dir[d].data[ix[i]].den_mutations[j];
                            break;
                            
                        case DEN_READ_DEPTH:
                            mrg->data[ix[i]].den_read_depth[j] += an_dir[d].data[ix[i]].den_read_depth[j];
                            break;
                            
                        case DEN_EFF_DEPTH:
                            mrg->data[ix[i]].den_eff_depth[j] += an_dir[d].data[ix[i]].den_eff_depth[j];
                            break;
                            
                        case DEN_OFF_TARGET_DEPTH:
                            mrg->data[ix[i]].den_off_target_depth[j] += an_dir[d].data[ix[i]].den_off_target_depth[j];
                            break;
                            
                        case DEN_LOW_MAPQ_DEPTH:
                            mrg->data[ix[i]].den_low_mapq_depth[j] += an_dir[d].data[ix[i]].den_low_mapq_depth[j];
                            break;
                            
                        case DEN_MAPPED_DEPTH:
                            mrg->data[ix[i]].den_mapped_depth[j] += an_dir[d].data[ix[i]].den_mapped_depth[j];
                            break;
                            
                        default:
                            printf("store_SM2_profile: error - exceeded expected number of profile columns. aborting...\n");
                            abort();
                            break;
                    }
                }
                
            }
            
            //calculate mutation rates using summed mutation/eff depth values
            mrg->data[ix[i]].mod_rate[j] = ((double)mrg->data[ix[i]].mod_mutations[j])/((double)mrg->data[ix[i]].mod_eff_depth[j]);
            mrg->data[ix[i]].unt_rate[j] = ((double)mrg->data[ix[i]].unt_mutations[j])/((double)mrg->data[ix[i]].unt_eff_depth[j]);
            mrg->data[ix[i]].den_rate[j] = ((double)mrg->data[ix[i]].den_mutations[j])/((double)mrg->data[ix[i]].den_eff_depth[j]);
            
            //calculate reactivity
            if (mrg->chnls.mod && !mrg->chnls.unt && !mrg->chnls.den) {
                mrg->data[ix[i]].reactivity_profile[j] = mrg->data[ix[i]].mod_rate[j];
            
            } else if (mrg->chnls.mod && mrg->chnls.unt && !mrg->chnls.den) {
                mrg->data[ix[i]].reactivity_profile[j] = mrg->data[ix[i]].mod_rate[j] - mrg->data[ix[i]].unt_rate[j];
                                
            } else if (mrg->chnls.mod && mrg->chnls.unt && mrg->chnls.den) {
                mrg->data[ix[i]].reactivity_profile[j] = (mrg->data[ix[i]].mod_rate[j] - mrg->data[ix[i]].unt_rate[j])/mrg->data[ix[i]].den_rate[j];
                
            } else {
                printf("merge_profiles: error - unexpected channel configuration. aborting...\n");
                abort();
            }
            
            //force reactivity values that are 0 within the
            //value that will be printed to be true zero
            if (fabs(mrg->data[ix[i]].reactivity_profile[j]) < 0.000001) {
                mrg->data[ix[i]].reactivity_profile[j] = 0;
            }
            
            //set hq profile value
            if (isHQnuc(&mrg->data[ix[i]], j, min_depth, max_bkg)) {
                mrg->data[ix[i]].hq_profile[j] = mrg->data[ix[i]].reactivity_profile[j];
                
            } else {
                mrg->data[ix[i]].hq_profile[j] = NAN;
            }   
        }
        
        mrg->data[ix[i]].sequence[j] = '\0'; //terminate sequence string
        
    }
    return 1;
}




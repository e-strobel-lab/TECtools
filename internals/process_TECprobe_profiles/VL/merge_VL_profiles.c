//
//  merge_VL_profiles.c
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
#include "../../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"



#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"
#include "../global/store_SM2_profile.h" //TODO: temporary to allow compiling before revision of the function
#include "../global/calculate_normalization_factor.h"

#include "merge_VL_profiles.h"

/* merge_VL_profiles: combine data from input reactivity profiles to generate merged reactivity profile files */
int merge_VL_profiles(SM2_analysis_directory * an_dir, int dir_count, SM2_analysis_directory * mrg, int min_depth, double max_bkg)
{
    
    int tl  = 0;   //transcript length index
    int i   = 0;   //general purpose index
    int col = 0;   //column index
    int d   = 0;   //directory index
    
    //set variables for merged analysis directory to that of the first input analysis directory.
    //these variables are identical for all input analysis directories
    mrg->trgt_start = an_dir[0].trgt_start;   //set target start index
    mrg->min_tl = an_dir[0].min_tl;           //set minimum transcript length
    mrg->max_tl = an_dir[0].max_tl;           //set maximum transcript length
    mrg->trg_rct_cnt = an_dir[0].trg_rct_cnt; //set target nucleotide reactivity count
    mrg->chnls.mod = an_dir[0].chnls.mod;     //set modified  channel detection flag
    mrg->chnls.unt = an_dir[0].chnls.unt;     //set untreated channel detection flag
    mrg->chnls.den = an_dir[0].chnls.den;     //set denatured channel detection flag
        
    for (tl = an_dir[0].min_tl; tl <= an_dir[0].max_tl; tl++) {  //for every transcript length

        //allocate memory for storing the merged profile of the current transcript length
        allocate_SM2_profile_memory(&mrg->data[tl], an_dir[0].data[tl].tot_nt_cnt);

        //TODO: check that tot nucleotide count is the same for all inputs
        for (i = 0; i < an_dir[0].data[tl].tot_nt_cnt; i++) { //for every line of the profile file
            
            mrg->data[tl].trgt_start  = an_dir[0].data[tl].trgt_start;
            mrg->data[tl].trg_nt_cnt = an_dir[0].data[tl].trg_nt_cnt;
            mrg->data[tl].tot_nt_cnt = an_dir[0].data[tl].tot_nt_cnt;
            
            //set channel tracker flags
            mrg->data[tl].chnls.mod = an_dir[0].chnls.mod;
            mrg->data[tl].chnls.unt = an_dir[0].chnls.unt;
            mrg->data[tl].chnls.den = an_dir[0].chnls.den;
            
            
            mrg->data[tl].nucleotide[i] = an_dir[0].data[tl].nucleotide[i]; //copy nucleotide number
            mrg->data[tl].sequence[i] = an_dir[0].data[tl].sequence[i];     //copy sequence
                        
            mrg->data[tl].std_err[i] = NAN;     //mask std err
            mrg->data[tl].hq_stderr[i] = NAN;   //mask hq std err
            mrg->data[tl].norm_stderr[i] = NAN; //mask normalized std err
                        
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
                            mrg->data[tl].mod_mutations[i] += an_dir[d].data[tl].mod_mutations[i];
                            break;
                            
                        case MOD_READ_DEPTH:
                            mrg->data[tl].mod_read_depth[i] += an_dir[d].data[tl].mod_read_depth[i];
                            break;
                            
                        case MOD_EFF_DEPTH:
                            mrg->data[tl].mod_eff_depth[i] += an_dir[d].data[tl].mod_eff_depth[i];
                            break;
                            
                        case MOD_OFF_TARGET_DEPTH:
                            mrg->data[tl].mod_off_target_depth[i] += an_dir[d].data[tl].mod_off_target_depth[i];
                            break;
                            
                        case MOD_LOW_MAPQ_DEPTH:
                            mrg->data[tl].mod_low_mapq_depth[i] += an_dir[d].data[tl].mod_low_mapq_depth[i];
                            break;
                            
                        case MOD_MAPPED_DEPTH:
                            mrg->data[tl].mod_mapped_depth[i] += an_dir[d].data[tl].mod_mapped_depth[i];
                            break;
                            
                        case UNT_MUTATIONS:
                            mrg->data[tl].unt_mutations[i] += an_dir[d].data[tl].unt_mutations[i];
                            break;
                            
                        case UNT_READ_DEPTH:
                            mrg->data[tl].unt_read_depth[i] += an_dir[d].data[tl].unt_read_depth[i];
                            break;
                            
                        case UNT_EFF_DEPTH:
                            mrg->data[tl].unt_eff_depth[i] += an_dir[d].data[tl].unt_eff_depth[i];
                            break;
                            
                        case UNT_OFF_TARGET_DEPTH:
                            mrg->data[tl].unt_off_target_depth[i] += an_dir[d].data[tl].unt_off_target_depth[i];
                            break;
                            
                        case UNT_LOW_MAPQ_DEPTH:
                            mrg->data[tl].unt_low_mapq_depth[i] += an_dir[d].data[tl].unt_low_mapq_depth[i];
                            break;
                            
                        case UNT_MAPPED_DEPTH:
                            mrg->data[tl].unt_mapped_depth[i] += an_dir[d].data[tl].unt_mapped_depth[i];
                            break;
                            
                        case DEN_MUTATIONS:
                            mrg->data[tl].den_mutations[i] += an_dir[d].data[tl].den_mutations[i];
                            break;
                            
                        case DEN_READ_DEPTH:
                            mrg->data[tl].den_read_depth[i] += an_dir[d].data[tl].den_read_depth[i];
                            break;
                            
                        case DEN_EFF_DEPTH:
                            mrg->data[tl].den_eff_depth[i] += an_dir[d].data[tl].den_eff_depth[i];
                            break;
                            
                        case DEN_OFF_TARGET_DEPTH:
                            mrg->data[tl].den_off_target_depth[i] += an_dir[d].data[tl].den_off_target_depth[i];
                            break;
                            
                        case DEN_LOW_MAPQ_DEPTH:
                            mrg->data[tl].den_low_mapq_depth[i] += an_dir[d].data[tl].den_low_mapq_depth[i];
                            break;
                            
                        case DEN_MAPPED_DEPTH:
                            mrg->data[tl].den_mapped_depth[i] += an_dir[d].data[tl].den_mapped_depth[i];
                            break;
                            
                        default:
                            printf("store_SM2_profile: error - exceeded expected number of profile columns. aborting...\n");
                            abort();
                            break;
                    }
                }
                
            }
            
            //calculate mutation rates using summed mutation/eff depth values
            mrg->data[tl].mod_rate[i] = ((double)mrg->data[tl].mod_mutations[i])/((double)mrg->data[tl].mod_eff_depth[i]);
            mrg->data[tl].unt_rate[i] = ((double)mrg->data[tl].unt_mutations[i])/((double)mrg->data[tl].unt_eff_depth[i]);
            mrg->data[tl].den_rate[i] = ((double)mrg->data[tl].den_mutations[i])/((double)mrg->data[tl].den_eff_depth[i]);
            
            //calculate reactivity
            if (mrg->chnls.mod && !mrg->chnls.unt && !mrg->chnls.den) {
                mrg->data[tl].reactivity_profile[i] = mrg->data[tl].mod_rate[i];
            
            } else if (mrg->chnls.mod && mrg->chnls.unt && !mrg->chnls.den) {
                mrg->data[tl].reactivity_profile[i] = mrg->data[tl].mod_rate[i] - mrg->data[tl].unt_rate[i];
                                
            } else if (mrg->chnls.mod && mrg->chnls.unt && mrg->chnls.den) {
                mrg->data[tl].reactivity_profile[i] = (mrg->data[tl].mod_rate[i] - mrg->data[tl].unt_rate[i])/mrg->data[tl].den_rate[i];
                
            } else {
                printf("merge_VL_profiles: error - unexpected channel configuration. aborting...\n");
                abort();
            }
            
            //force reactivity values that are 0 within the
            //value that will be printed to be true zero
            if (fabs(mrg->data[tl].reactivity_profile[i]) < 0.000001) {
                mrg->data[tl].reactivity_profile[i] = 0;
            }
            
            //set hq profile value
            if (isHQnuc(&mrg->data[tl], i, min_depth, max_bkg)) {
                mrg->data[tl].hq_profile[i] = mrg->data[tl].reactivity_profile[i];
                
            } else {
                mrg->data[tl].hq_profile[i] = NAN;
            }   
        }
        
        mrg->data[tl].sequence[i] = '\0'; //terminate sequence string
        
    }
    return 1;
}




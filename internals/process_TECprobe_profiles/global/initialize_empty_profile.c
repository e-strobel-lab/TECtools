//
//  initialize_empty_profile.c
//  
//
//  Created by Eric Strobel on 3/15/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../../global/global_defs.h"
#include "../../utils/debug.h"
#include "../../utils/io_management.h"

#include "./store_SM2_profile.h"
#include "initialize_empty_profile.h"

/* initialize_empty_profile: initialize values for empty profile structure */
void initialize_empty_profile(SM2_profile * prf, int trg_nt_cnt, int trgt_start_ix)
{
    int i = 0; //general purpose index
    
    prf->trgt_start = trgt_start_ix;              //set target start index
    prf->trg_nt_cnt = trg_nt_cnt;                 //set target nucleotide count
    prf->tot_nt_cnt = trg_nt_cnt + trgt_start_ix; //set total  nucleotide count
    
    //set channels
    prf->chnls.mod = 0;
    prf->chnls.unt = 0;
    prf->chnls.den = 0;
    
    //set profile values
    for (i = 0; i < prf->tot_nt_cnt; i++) {
        prf->nucleotide[i] = i+1;
        prf->sequence[i] = (i < trgt_start_ix) ? 'n' : 'N';
        prf->mod_mutations[i] = 0;
        prf->mod_read_depth[i] = 0;
        prf->mod_eff_depth[i] = 0;
        prf->mod_rate[i] = NAN;
        prf->mod_off_target_depth[i] = 0;
        prf->mod_low_mapq_depth[i] = 0;
        prf->mod_mapped_depth[i] = 0;
        prf->unt_mutations[i] = 0;
        prf->unt_read_depth[i] = 0;
        prf->unt_eff_depth[i] = 0;
        prf->unt_rate[i] = NAN;
        prf->unt_off_target_depth[i] = 0;
        prf->unt_low_mapq_depth[i] = 0;
        prf->unt_mapped_depth[i] = 0;
        prf->den_mutations[i] = 0;
        prf->den_read_depth[i] = 0;
        prf->den_eff_depth[i] = 0;
        prf->den_rate[i] = NAN;
        prf->den_off_target_depth[i] = 0;
        prf->den_low_mapq_depth[i] = 0;
        prf->den_mapped_depth[i] = 0;
        prf->reactivity_profile[i] = NAN;
        prf->std_err[i] = NAN;
        prf->hq_profile[i] = NAN;
        prf->hq_stderr[i] = NAN;
        prf->norm_profile[i] = NAN;
        prf->norm_stderr[i] = NAN;
    }
    
    prf->sequence[i] = '\0';
}

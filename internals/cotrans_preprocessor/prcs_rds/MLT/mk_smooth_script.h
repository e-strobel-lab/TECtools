//
//  mk_smooth_script.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef mk_smooth_script_h
#define mk_smooth_script_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

/* mk_smooth_script: generate shell script for concatenating 3' end split files
 the minimum target is concatenated as: n, n+1.
 the maximum target is concatenated as: n-1, n.
 intermediate targets are concatenated as: n-1, n, n+1
 */

void mk_smooth_script(names * nm, target3p_params trg_prms);	//generate script for transcript length smoothing

#endif /* mk_smooth_script_h */

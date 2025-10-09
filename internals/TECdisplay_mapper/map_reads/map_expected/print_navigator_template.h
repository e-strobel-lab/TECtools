//
//  print_navigator_template.h
//  
//
//  Created by Eric Strobel on 10/9/25.
//

#ifndef print_navigator_template_h
#define print_navigator_template_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../TECdisplay_mapper_structs.h"
#include "../../TECdisplay_mapper_defs.h"

#include "../../../utils/io_management.h"
#include "../../../utils/gen_utils.h"

/* print_navigator_template: make template files for TECdisplay_navigator */
void print_navigator_template(target *refs, TDSPLY_fasta * wt, target_params * trg_prms);

#endif /* print_navigator_template_h */

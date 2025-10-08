//
//  mk_barcoded_target_fastas.h
//  
//
//  Created by Eric Strobel on 10/6/25.
//

#ifndef mk_barcoded_target_fastas_h
#define mk_barcoded_target_fastas_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/io_management.h"

#include "./mk_MUX_trgts.h"

/* mk_barcoded_target_fastas: generate directory that contains a fasta file for each target and a file that contains a list of barcode ids */
void mk_barcoded_target_fastas(TPROBE_names * nm, compact_target * ctrg, target_params * trg_prms);

#endif /* mk_barcoded_target_fastas_h */

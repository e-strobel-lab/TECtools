//
//  read_vbase.h
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#ifndef read_vbase_h
#define read_vbase_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/is_dgnrt_mtch.h"
#include "../seq_utils/basemap.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"

/* read_vbase: parse variable base information from an input string and store in a base_params structure.
 read_vbase performs extensive checks to confirm the validity of the parsed variable base with respect to
 the reference sequence and the variable base reference sequence.
 
 format:
   variable base, no insertion: <pos><seq>          e.g., 19N
   variable base, w/ insertion: <pos>i<ins><seq>    e.g., 19i1N
   constant base:               c<pos><seq>         e.g., c19A
   deletion site:               d<pos><seq>         e.g., d19
 */
char * read_vbase(char * vbase, char delimiter, base_params * prsd_vbase, basemap * bmap, int mode);


/* print_read_vbase_error: print error message failure to match vbase to variable base reference sequence */
void print_read_vbase_error(char * vbase, base_params * prsd_vbase, basemap * bmap, uint64_t bi, int typ_error, int pos_error, int ins_error, int seq_error);

#endif /* read_vbase_h */

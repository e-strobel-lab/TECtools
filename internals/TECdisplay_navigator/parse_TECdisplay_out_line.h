//
//  parse_TECdisplay_out_line.h
//  
//
//  Created by Eric Strobel on 4/4/24.
//

#ifndef parse_TECdisplay_out_line_h
#define parse_TECdisplay_out_line_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "../global/global_defs.h"

#include "../TECdisplay_mapper/TECdisplay_output_column_headers.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"

#define TDSPLY_HDR_LINE  0
#define TDSPLY_DATA_LINE 1

void parse_TECdisplay_out_line(char * line, char ** p_id, char ** p_vals, int * bnd, int * unb, double * frc, int mode, int nonstandard);

void check_header_string(char * haystack, const char * needle);

#endif /* parse_TECdisplay_out_line_h */

//
//  mk_output_files.h
//  
//
//  Created by Eric Strobel on 10/11/22.
//

#ifndef mk_output_files_h
#define mk_output_files_h

#include <stdio.h>
#include <stdlib.h>

#include "../global/global_defs.h"

#include "./cotrans_mtrx.h"
#include "./mkmtrx_defs.h"
#include "./mkmtrx_structs.h"

/* print_csv: print csv file of value matrix */
int print_csv(char * prnt_dir, char * smpl_nm, char * output_nm, char * (*mtrx)[MAX_ROW][MAX_COL], int row_min, int row_max, int col_min, int col_max);

/* print_alignment stats: print summary of alignment states to files and screen.
 per length alignment rates are sent to <smpl_nm>_alignment_rates.txt
 overall aligment totals are setn to <smpl_nm>_alignment_totals.txt */
int print_alignment_stats(char * prnt_dir, char * smpl_nm, alignment_stats * algn, cotrans_matrix * mtrx);

int columnize_mtrx(char * prnt_dir, char * smpl_nm, cotrans_matrix * mtrx);


int mk_linebars(char * prnt_dir, char * smpl_nm, cotrans_matrix * mtrx);

#endif /* mk_output_files_h */

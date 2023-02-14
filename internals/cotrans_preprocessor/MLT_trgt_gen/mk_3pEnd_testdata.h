//
//  mk_3pEnd_testdata.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef mk_3pEnd_testdata_h
#define mk_3pEnd_testdata_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../seq_utils/revcomp.h"

/* mk_3pEnd_testdata: generate individual fasta files for every intermediate
 length of the input sequence, starting at min_len+1
 ***arguments***
 char name[MAX_LINE]: sequence name
 char seq[MAX_LINE]: input sequence
 target3p_genVals * ends: end target generation values
 int randomize_end: flag to mutate 3' end region when making test data
 int mltplir: multiplier value for generating more than the default amount of test data
 int min_len: minimum transcript length-1
 */
int mk_3pEnd_testdata(char name[MAX_LINE], char seq[MAX_LINE], target3p_genVals * ends, int randomize_end, int mltplir, int min_len);


/* mk_rndmzd_bc: generate a randomized channel barcode with variable channel and match settings.
 chnl variable:
 UNT = untreated barcode (YYYRR)
 MOD = modified barcode  (RRRYY)
 ERR = undetermined barcode (<4/5 position match to expected sequence)
 
 mtch variable:
 FULL = 5/5 position match to expected sequence
 PART = 4/5 position match to expected sequence
 
 ***arguments***
 char * bc: array to store randomized barcode
 int chnl: barcode generation channel setting
 int mtch: barcode generation match setting
 */
void mk_rndmzd_bc(char * bc, int chnl, int mtch);


/* rndmz_end: randomize 3' end region  of testdata insert sequence
 ***arguments***
 char * out: array to store randomized output sequence
 char * ipt: input sequence (substring of seq)
 char * seq: full input sequence
 target3p_genVals * ends: end target generation values
 int min_len: minimum transcript length-1
 */
int rndmz_end(char * out, char * ipt, char * seq, target3p_genVals * ends, int min_len);


/* print_fq: construct read sequences and print to fastq file
 ***arguments***
 FILE * out_rd1: file pointer for read 1 output
 FILE * out_rd2: file pointer for read 2 output
 char * insrt2use: insert for test data read generation
 char * chnl_bc: channel barcode for test data read generation
 int end3p: transcript length of sequence that was used to generate insert
 int end_rnd_typ: randomization type
 */
void print_fq(FILE * out_rd1, FILE * out_rd2, char * insrt2use, char * chnl_bc, int end3p, int end_rnd_typ);


/* print_var_map: print map that illustrates 3' end randomization
 this is for debugging purposes and is not typically used. it is
 not accessible unless the code is edited directly.
 ***arguments***
 char * ipt: input sequence
 char * out: randomized output sequence
 int rnd_pos: position that was randomized
 int mode: randomization type
 */
void print_var_map(char * ipt, char * out, int rnd_pos, int mode);

#endif /* mk_3pEnd_testdata_h */

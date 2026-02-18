//
//  mk_SGL_target.c
//  
//
//  Created by Eric Strobel on 1/17/23.
//

#include "mk_SGL_target.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

#include "../../utils/io_management.h"
//mk_MLT_targets is only included for the parse_fasta
//function. should move parse_fasta elsewhere
#include "../MLT_trgt_gen/mk_MLT_trgts.h"


int mk_SGL_target(FILE * fp_fasta)
{
    FILE * out_fp = {NULL};        //FILE pointer for generating output fasta files
    
    int i = 0; //general purpose index

    //out_dir/outname:
    //+9 accounts for extra printed chars to out_dir,
    //+14 for extra printed chars to out_name
    //+3 for transcript length digits
    //+1 for terminating null
    char out_dir[MAX_LINE+9+1] = {0};           //array for output directory name
    char out_name[(MAX_LINE*2)+9+14+3+1] = {0}; //array for output fasta filename
    
    static const char ldr[21]  = "atggccttcgggccaa"; //leader sequence, masked by lowercase

    char fa_name[MAX_LINE] = {0};  //fasta sequence name
    char seq[MAX_LINE] = {0};      //fasta sequence
    int ipt_len = 0;               //length of input sequence
    char target[MAX_LINE] = {0};   //array to store target sequence
    
    if (!parse_fasta(fp_fasta, &fa_name[0], &seq[0])) { //parse fasta file
        printf("mk_SGL_target: error - input fasta file was empty. aborting...\n");
        abort();
    }
    
    ipt_len = strlen(seq); //get sequence length
    
    //copy target sequence to new array while converting all characters to uppercase
    for (i = 0; seq[i]; i++) {
        target[i] = toupper(seq[i]); //ensure target is uppercase
    }
    target[i] = '\0';
    
    sprintf(out_dir, "./%s_target", fa_name);  //assemble target directory name
    mk_out_dir(out_dir);                       //make output directory for target file
    
    sprintf(out_name, "%s/%s_%dnt_target.fa", out_dir, fa_name, ipt_len); //assemble target file name
    if ((out_fp = fopen(out_name,"w")) == NULL) {
        printf("mk_SGL_target: error - could not open target output file. Aborting program...\n");
        abort();
    }
    
    fprintf(out_fp, ">%s\n", fa_name);         //print target name
    fprintf(out_fp, "%s%s\n", ldr, target); //print leader+target sequence
    
    if ((fclose(out_fp)) == EOF) {
        printf("mk_SGL_trgt: error - error occurred when closing target output file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

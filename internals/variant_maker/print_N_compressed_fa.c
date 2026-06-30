//
//  print_N_compressed_fa.c
//  
//
//  Created by Eric Strobel on 6/19/26.
//

#include "print_N_compressed_fa.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../seq_utils/basemap.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"
#include "./print_output.h"

/* print_N_compressed_fa: print fasta file in which variants are compressed by replacing nucleotides that can be any base with 'N' */
void print_N_compressed_fa(names * nm, basemap * bmap, int vTmpCnt, int varCnt, char * out_dir, int append_priming, int append_barcode, char * lnkr, int make_fasta, int lib_type)
{
    extern fasta *vrnts;       //array of variant sequences
    extern uint64_t v_indx;    //index for vrnts array
    
    if (lib_type != TECDISPLAY_LIB || append_barcode) {
        printf("print_N_compressed_fa: error - N compression can only be applied to unbarcoded TECdisplay variant libraries. aborting...\n");
        abort();
    }
    
    FILE * Ncmp_fp = NULL;        //output file pointer for n-compressed fasta
    char Ncmp_nm[MAX_LINE] = {0}; //output file name
    int ret = 0;                  //variable for storing snprintf return value
    
    //construct n-compressed fasta file name
    ret = snprintf(Ncmp_nm, MAX_LINE, "./%s/%s_variants_Ncompressed.fa", out_dir, nm->vTmp);
    if (ret >= MAX_LINE || ret < 0) {
        printf("print_N_compressed_fa: error - error when constructing variant output file name. aborting...\n");
        abort();
    }
    
    //generate output file
    if ((Ncmp_fp = fopen(Ncmp_nm, "w")) == NULL) {
        printf("print_N_compressed_fa: ERROR - could not generate variants output file. aborting...\n");
        abort();
    }
        
    int b = 0; //basemap index
    int v = 0; //variant index
    int n = 0; //n-compressed variant index
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
    
    int len = 0;     //sequence length
    int pr_char = 0; //flag that current nucleotide is paired
    
    int Ncmp_var_cnt = 0; //number of n-compressed variants
    
    fasta * Ncmp_vrnts = NULL; //pointer for allocating fasta structures for n-compressed variants
    
    //allocate memory for N-compressed variants
    if ((Ncmp_vrnts = calloc(varCnt, sizeof(*Ncmp_vrnts))) == NULL) {
        printf("print_N_compressed_fa: error - failed to allocate memory for N compression variants. aborting...\n");
        abort();
    }
    
    for (b = 0, v = 0, n = 0; b < vTmpCnt; b++) {  //for each basemap
        
        //check that the first entry of vrnts is a reference sequence
        if (strstr(vrnts[v].nm, "REF") == NULL) { //if not, abort
            printf("print_N_compressed_fa: error - expected first variant for variant template %d to be a reference sequence. aborting...\n", b);
            abort();
            
        } else if (!strcmp(vrnts[v].nm, bmap[b].rS)) { //if variants and basemap are not aligned, abort
            printf("print_N_compressed_fa: error - basemap and variant list are not aligned. aborting...\n");
            abort();
            
        } else { //increment variant index
            v++;
        }
    
        len = strlen(bmap[b].rS); //set variant length using reference sequence
        
        for (i = 0; i < bmap[b].cnt[FILTERED]; i++, v++, n++) { //for every variant associate with the current basemap
            
            //allocate memory for storing n-compressed variant sequence
            if ((Ncmp_vrnts[n].sq = malloc((len+1) * sizeof(Ncmp_vrnts[n].sq))) == NULL) {
                printf("print_N_compressed_fa: error - failed to allocate memory for N compressed variant sequence. aborting...\n");
                abort();
            }
            
            for (j = 0; j < len; j++) { //for each character in the current variant sequence
                
                for (k = 0, pr_char = 0; k < bmap[b].rP_cnt; k++) { //for each ref pairing constraint associated with the current variant
                    
                    if (bmap[b].rP[k][j] == '(' || bmap[b].rP[k][j] == ')') { //if a pair is indicated
                        pr_char = 1;                                          //set pr_char flag
                    }
                }
                
                if ((bmap[b].rS[j] == 'N' || bmap[b].rS[j] == 'n') && !pr_char) { //if ref seq is an unpaired 'n'
                    Ncmp_vrnts[n].sq[j] = 'N';                                    //substitute variant base with N in N-compressed seq
                } else {
                    Ncmp_vrnts[n].sq[j] = vrnts[v].sq[j];                         //otherwise, do not substitute variant base
                }
            }
            
            Ncmp_vrnts[n].sq[j] = '\0'; //terminate variant sequence string
            Ncmp_var_cnt++;             //increment Ncmp_var_cnt
        }
    }
    
    //generate array for tracking whether a variant or its identical sequence has been printed
    int * prntd = NULL;
    if ((prntd = calloc(varCnt, sizeof(*prntd))) == NULL) {
        printf("print_N_compressed_fa: failed to allocate memory for printed line flags. aborting...\n");
        abort();
    }
    
    int non_rdndnt = 0; //tracks the number of non-redundant variant sequences that have been printed
        
    for (i = 0; i < Ncmp_var_cnt; i++) { //for every n-compressed variant
                
        if (!prntd[i]) {  //if the variant was not yet printed

            non_rdndnt++; //increment non_rdndnt then print variant to file using print_standard_variant function
            
            print_standard_variant(NULL, Ncmp_fp, &Ncmp_vrnts[i], non_rdndnt, append_priming, lib_type, 0, make_fasta);
            
            for (j = i+1; j < Ncmp_var_cnt; j++) {                     //for every variant after the current variant
                if (!prntd[j]) {                                       //if the variant has not yet been printed
                    if (!strcmp(Ncmp_vrnts[j].sq, Ncmp_vrnts[i].sq)) { //check if the variant seq matches that of the current variant
                        prntd[j] = 1;                                  //if so, set flag that the variant was already printed
                    }
                }
            }
        }
    }
    
    //close n-compressed fasta output file
    if (fclose(Ncmp_fp) == EOF) {
        printf("print_N_compressed_fa: error - error occurred when closing N-compressed fasta file. Aborting program...\n");
        abort();
    }
}



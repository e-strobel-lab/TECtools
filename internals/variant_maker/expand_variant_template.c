//
//  expand_variant_template.c
//  
//
//  Created by Eric Strobel on 8/3/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"

#include "../seq_utils/isIUPACbase.h"

#include "./mk_variant_nm.h"
#include "./make_barcodes.h"
#include "./filter_variants.h"

#include "expand_variant_template.h"

/* expand_variant_template: recursively construct targets from a variant template base map */
void expand_variant_template(basemap * bmap, int nxt, char output[MAXLEN+1], struct trgt_vb lcl_bases, int mode)
{
    extern fasta *vrnts;    //array of variant sequences
    extern uint64_t v_indx; //index for vrnts array
    
    extern barcode_bank * bc_crnt; //pointer to current barcode bank
    extern uint64_t barcode_count; //number of barcodes that passed filter
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    struct varbase tmp_vb; //temporary copy of vb_tbl entry, used when branching recursion path
    
    /* increment local bases index if variable base was previously found. this code
     differentiates the initial expand_variant_template call, in which no variable
     bases have been found yet and the local bases index should remain zero, from
     all recursive calls of expand_variant_template, which only occur when a variable
     base has been found.
     */
    if (lcl_bases.bs[0]) { //if variable base was previously found
        lcl_bases.cnt++;   //increment lcl_bases.cnt
    }
        
    //copy sequence to output until variable base is found
    for (i = nxt; bmap->nts[i] != '*' && bmap->nts[i] && i < MAXLEN; i++) {
        
        if (mode == MAKE_BARCODES && bmap->nts[i] == 'P') { //nucleotide is the second base of a simple pair
            
            /* if the basemap was generated using a barcode template, the second base
             simple pairs is indicated by a 'P' in the nts array. the value of these
             nucleotides in the output sequence is the complement of the first base
             of the simple pair.
             */
            
            if (bmap->rP_cnt != 1) { //check that barcode template only contains 1 secondary structure
                printf("expand_variant_template: error - barcode template basemaps must contain exactly one pair constraint. aborting...\n");
                abort();
            }

            //set sequence to complement of the simple pair first mate
            switch (output[bmap->prs[0][i]]) {
                case 'A': output[i] = 'T'; break;
                case 'T': output[i] = 'A'; break;
                case 'G': output[i] = 'C'; break;
                case 'C': output[i] = 'G'; break;
                default: break;
            }
            
        } else { //the character is a simple DNA base or is a spacer
            output[i] = bmap->nts[i]; //set output to corresponding char in nts array
        }
    }
    
    /* set variable base and branch the recursion path */
    if (bmap->nts[i] == '*') { //current nt is a variable base with no dependence on previously determined variable base
        
        /* copy vb_tbl entry to tmp_vb. this generates a 'workng' copy
         of the current vb_tbl entry that can be used to iterate through
         each possible variable base when branching the recursion path
         */
        tmp_vb = *bmap->p_vb[i];
        
        /* set variable base and branch recursion path. in each iteration of the loop, the current
         nucleotides array index is recorded and the next local base is set to the value of the current
         alternate base that is specified by the tmp_vb entry. expand_variant_template is then called
         recursively to proceed expansion of the variant template
         */
        while (tmp_vb.alt[tmp_vb.crnt]) {             //repeat until all vb_tbl possibilities have been applied
            output[i] = tmp_vb.alt[tmp_vb.crnt];      //set output[i] to current vb_tbl entry
            lcl_bases.ix[lcl_bases.cnt] = i;          //record index of variable base
            lcl_bases.bs[lcl_bases.cnt] = output[i];  //record current variable base entry in this recursion path
            
            //recursively call expand_variant_template to generate all variants
            expand_variant_template(bmap, i+1, output, lcl_bases, mode);
            tmp_vb.crnt++;
        }
        
    } else if (bmap->nts[i]) { //unrecognized base error
        printf("expand_variant_template: unexpected character %c. aborting...\n", bmap->nts[i]);
        abort();
    }
    
    char vrnt_nm[256] = {0}; //array for storing variant name
    
    if (!bmap->nts[i]) {           //reached the end of vbase string
        bmap->cnt[EXPANDED]++;     //increment number of variants that resulted from variant template expansion
        
        /* in MAKE_VARIANTS mode, check that the variant meets all
         base pairing constraints that are specified by the basemap.
         variants that pass this filter are stored variants that do
         not pass this filter are discarded.
         */
        if (mode == MAKE_VARIANTS && !filter_variants(output, bmap)) {  //test if variant passes filter
            bmap->cnt[FILTERED]++;                                      //increment number of variants that passed filter
            mk_variant_name(vrnt_nm, lcl_bases, bmap);                  //construct variant name
        
            //allocate memory for variant name and sequence
            vrnts[v_indx].nm = malloc((strlen(vrnt_nm)+1) * sizeof(*(vrnts[v_indx].nm)));
            vrnts[v_indx].sq = malloc((strlen(output)+1) * sizeof(*(vrnts[v_indx].sq)));
            
            if (vrnts[v_indx].nm == NULL || vrnts[v_indx].sq == NULL) { //check memory allocation success
                printf("expand_variant_template: variant memory allocation failed. aborting...");
                abort();
            }
            
            strcpy(vrnts[v_indx].nm, vrnt_nm); //store variant name in vrnts array
            strcpy(vrnts[v_indx].sq, output);  //store variant sequence in vrnts array
            v_indx++;                          //increment variant index
            
        /* in MAKE_BARCODES mode, check that the variant meets all barcode
         specifications. variants that pass this filter are stored. variants
         that do not pass this filter are discarded
         */
        } else if (mode == MAKE_BARCODES && !filter_barcode(output)) {
            
            //allocate memory for barcode sequence
            if ((bc_crnt->fa[bc_crnt->cnt].sq = malloc((strlen(output)+1) * sizeof(*(bc_crnt->fa[bc_crnt->cnt].sq)))) == NULL) {
                printf("expand_variant_template: barcode sequence memory allocation failed. aborting...");
                abort();
            }
            
            strcpy(bc_crnt->fa[bc_crnt->cnt].sq, output); //store barcode sequence
            
            if (++bc_crnt->cnt == BC_BLOCK_SIZE) {        //if current bank is now full
                add_barcode_bank(bc_crnt);                //allocate a new barcode bank
                bc_crnt = bc_crnt->nxt;                   //set bc_crnt to point to the next barcode bank
            }
            
            barcode_count++; //increment total passed filter barcode count
        }
    }
    
    return;
}

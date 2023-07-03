//
//  count_variants.c
//  
//
//  Created by Eric Strobel on 8/4/22.
//

#include <stdio.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"

#include "../utils/gen_utils.h"
#include "../seq_utils/basemap.h"

#include "count_variants.h"

/* count_variants: determine the number of variants specified by the supplied variant templates
 variants are counted by multiplying current_count by the number of possible bases at each variable
 position unless the identity of the variable position is constrained by a previously counted
 variable base (this constraint only occurs in MAKE_BARCODES mode */
int count_variants(basemap * bmap, int tot_templates, int mode)
{
    extern FILE * prcs_ofp;            //output file pointer for processing messages
    extern char prcs_out_nm[MAX_LINE]; //name of processing message output file
    extern char out_msg[MAX_LINE];     //output message
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int total_count = 0;   //total number of variants encoded by all variant templates
    int current_count = 1; //number of variants encoded by the current variant template
                           //initialized to 1 because it is multiplied by the variable
                           //base composition of the sequence

    sprintf(out_msg, "\nvariant templates:\n\n"); //print first line of output
    printf2_scrn_n_fl(prcs_ofp, out_msg);
    
    for (i = 0; i < tot_templates; i++) { //for every variant template...
        
        sprintf(out_msg, "nm  %s_%s\n", bmap[i].wt->nm, bmap[i].nm); //print variant template name
        printf2_scrn_n_fl(prcs_ofp, out_msg);
        
        sprintf(out_msg, "seq %s\n", bmap[i].rS); //print variant template sequence
        printf2_scrn_n_fl(prcs_ofp, out_msg);
        
        if (mode == MAKE_VARIANTS) {              //in MAKE_VARIANTS mode
            for (j = 0; bmap[i].rP[j][0]; j++) {  //print all variant template pair constraints
                sprintf(out_msg, "str %s\n", bmap[i].rP[j]);
                printf2_scrn_n_fl(prcs_ofp, out_msg);
            }
            
        } else if (mode == MAKE_BARCODES) { //in MAKE_BARCODES mode
            
            if (bmap[i].rP_cnt != 1) { //check that barcode template only contains 1 secondary structure
                printf("count_variants: error - barcode template basemap must contain exactly 1 pair constraint. aborting...\n");
                abort();
            }
            
            sprintf(out_msg, "str %s\n", bmap[i].rP[0]); //print the barcode template pair constraint
            printf2_scrn_n_fl(prcs_ofp, out_msg);
        }
        
        current_count = 1;                 //set currrent_count to 1, then multiply based on variable bases below
        for (j = 0; bmap[i].rS[j]; j++) {  //for each base of the variant template sequence...
            
            //in MAKE_VARIANTS mode, all variable positions are counted at this stage.
            //in MAKE_BARCODES mode, only the first member of a variable base pair is
            //counted barcodes can only have on structure and all base pairs within
            //this structure are guaranteed to be WC pairs, so calculating the number
            //of variants that will remain after filtering is simple
            
            if ((mode == MAKE_VARIANTS) ||
                (mode == MAKE_BARCODES && bmap[i].rP[0][j] != ')')) {
                
                switch (toupper(bmap[i].rS[j])) {
                    case 'N': current_count *= 4; break;
                    case 'R': current_count *= 2; break;
                    case 'Y': current_count *= 2; break;
                    case 'M': current_count *= 2; break;
                    case 'K': current_count *= 2; break;
                    case 'S': current_count *= 2; break;
                    case 'W': current_count *= 2; break;
                    case 'B': current_count *= 3; break;
                    case 'D': current_count *= 3; break;
                    case 'H': current_count *= 3; break;
                    case 'V': current_count *= 3; break;
                    case 'A': break;
                    case 'T': break;
                    case 'G': break;
                    case 'C': break;
                    case '.': break;
                    case '-': break;
                    default:
                        printf("count_variants: error - unrecognized base %c. aborting...", bmap[i].rS[j]);
                        abort();
                        break;
                }
            }
        }
        
        bmap[i].cnt[CALCULATED] = current_count; //set calculated number of variants for current variant template
        total_count += current_count;            //add current count to total count
        
        //print calculated variants for current template
        if (mode == MAKE_VARIANTS) {
            sprintf(out_msg, "%d variants (without considering base pairing constraints)\n", bmap[i].cnt[CALCULATED]);
        } else {
            sprintf(out_msg, "%d variants\n\n", bmap[i].cnt[CALCULATED]);
        }
        printf2_scrn_n_fl(prcs_ofp, out_msg);
        print_vbases(&bmap[i]);
        printf("--------------------------------------------------------\n\n");
    }
    
    //print variant total
    if (mode == MAKE_VARIANTS) {
        sprintf(out_msg, "the above variant templates encode %d variants before base pairing constraints are applied.\n\n", total_count);
    } else {
        sprintf(out_msg, "the above variant templates encode %d variants.\n\n", total_count);
    }
    printf2_scrn_n_fl(prcs_ofp, out_msg);
    
    return total_count; //return variant total
}

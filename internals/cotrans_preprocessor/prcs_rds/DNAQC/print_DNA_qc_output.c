//
//  print_DNA_qc_output.c
//  
//
//  Created by Eric Strobel on 9/16/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

#include "../../../seq_utils/mapping_metrics.h"

#include "../../../variant_maker/bc_ind.h"

#include "print_DNA_qc_output.h"

/* print_DNA_qc_output: print DNA QC output files */
void print_DNA_qc_output(char * out_prfx, compact_target * ctrg, int brcd_cnt, int * pairing)
{
    FILE * ofp = NULL; //output file pointer
    
    int i = 0; //general purpose index

    compact_target * crnt_ctrg = NULL; //pointer to current ctrg
    compact_target * lnkd_ctrg = NULL; //pointer to linked ctrg
    
    opt_BC * crnt_BC_val = NULL; //pointer to current ctrg barcode values
    opt_BC * lnkd_BC_val = NULL; //pointer to linked ctrg barcode values
    
    association * assoc = NULL; //pointer to current ctrg association structure
    
    char * mrg = NULL; //pointer for generating merged barcode id
    
    //generate output file
    if ((ofp = fopen("concordance.txt", "w")) == NULL) {
        printf("assess_brcd_concordance: error - failed to generate output file. aborting...\n");
        abort();
    }
    
    //print concordance summary
    fprintf(ofp, "category\tcount\tfraction\n");
    fprintf(ofp, "total\t%d\n", pairing[CONCORDANT] + pairing[DISCORDANT]);
    fprintf(ofp, "concordant\t%d\t%.4f\n", pairing[CONCORDANT],
            ((double)(pairing[CONCORDANT])/((double)(pairing[CONCORDANT]) + (double)(pairing[DISCORDANT]))));
    fprintf(ofp, "discordant\t%d\t%.4f\n", pairing[DISCORDANT],
            ((double)(pairing[DISCORDANT])/((double)(pairing[CONCORDANT]) + (double)(pairing[DISCORDANT]))));
    
    //close output file
    if (fclose(ofp) == EOF) {
        printf("assess_brcd_concordance: error - failed to close output file\n");
    }
    
    //print concordance result to screen
    printf("\n%10d concordant (%6.2f%%)\n", pairing[CONCORDANT],
           100 * ((double)(pairing[CONCORDANT])/((double)(pairing[CONCORDANT]) + (double)(pairing[DISCORDANT]))));
    printf("%10d discordant (%6.2f%%)\n\n", pairing[DISCORDANT],
           100 * ((double)(pairing[DISCORDANT])/((double)(pairing[CONCORDANT]) + (double)(pairing[DISCORDANT]))));
        
    return;
}

/* merge_bc_id: merge the ids of two barcode targets */
void merge_bc_id(char ** mrg, char * str1, char * str2)
{
    extern char bc1_ind[3];
    
    char * id1 = NULL;       //pointer to start of barcode 1 id
    char * p_bc_ind1 = NULL; //pointer to barcode indicator 1
    char * p_bc_ind2 = NULL; //pointer to barcode indicator 2
    char * p_2bc = NULL;     //pointer to second barcode id, starting at '2bc'
    
    int id1_set = 0; //counts number of times id1 was set
    int id2_set = 0; //counts number of times id2 was set
    
    p_bc_ind1 = strstr(str1, bc1_ind); //set pointer to bc indicator in barcode 1 id
    p_bc_ind2 = strstr(str2, bc1_ind); //set pointer to bc indicator in barcode 2 id
    
    if (p_bc_ind1 == NULL || p_bc_ind2 == NULL) { //test whether barcode indicator strings were identified
        printf("merge_bc_id: error - failed to identify barcode indicator string in barcode id. aborting...\n");
        abort();
    }
    
    //determine if str1 is barcode 1 or 2
    if (p_bc_ind1[-1] == '_') { //if barcode indicator is preceded by '_', found barcode 1 id
        id1 = str1;             //point id1 to start of str1
        id1_set++;              //increment number of times id 1 was set
        
    } else if (p_bc_ind1[-1] == '2') { //if barcode indicator is preceded by '2', found barcode 2 id
        p_2bc = &p_bc_ind1[-1];        //set pointer to '2bc' in barcode 2 id
        id2_set++;                     //increment number of times barcode 2 was set
        
    } else { //barcode indicator not detected, throw error and abort
        printf("merge_bc_id: error - unrecognized id format. aborting...\n");
        abort();
    }
    
    //determine if str2 is barcode 1 or 2
    if (p_bc_ind2[-1] == '_') { //if barcode indicator is preceded by '_', found barcode 1 id
        id1 = str2;             //point id1 to start of str2
        id1_set++;              //increment number of times id 1 was set
        
    } else if (p_bc_ind2[-1] == '2') { //if barcode indicator is preceded by '2', found barcode 2 id
        p_2bc = &p_bc_ind2[-1];        //set pointer to '2bc' in barcode 2 id
        id2_set++;                     //increment number of times barcode 2 was set
    
    } else { //barcode indicator not detected, throw error and abort
        printf("merge_bc_id: error - unrecognized id format. aborting...\n");
        abort();
    }
    
    //merge barcode ids
    if (id1_set == 1 && id2_set == 1) { //if single barcode 1 and barcode 2 ids were found
        
        //allocate memory for storing the merged barcode id
        if ((*mrg = malloc((strlen(id1) + strlen(p_2bc) + 2) * sizeof(**mrg))) == NULL) {
            printf("merge_bc_id: error - failed to allocate memory for merged id. aborting...\n");
            abort();
        }
        sprintf(*mrg, "%s_%s", id1, p_2bc); //generate the merged barcode id
            
    } else { //if an incorrect number of barcode 1 and barcode 2 ids were found, throw error and abort
        printf("merge_bc_id: error - expected a pair of barcode 1 and barcode 2 ids to be provided. aborting...\n");
        abort();
    }
    
    return;
}

/* print_DNA_qc_metrics: print mapping metrics for DNA qc sample */
void print_DNA_qc_metrics(mapping_metrics * met)
{
    FILE * ofp = NULL; //output file poitner
    
    //generate mapping metrics file
    if ((ofp = fopen("mapping_matrics.txt", "w")) == NULL) {
        printf("print_DNA_qc_metrics: error - failed to generate mapping metrics file. aborting...\n");
        abort();
    }
    
    //print mapping metrics to output file
    fprintf(ofp, "%10d (%6.2f%%) barcodes mapped to whitelisted targets\n",
           met->mapped, ((double)(met->mapped)/(double)(met->reads_processed*2)) * 100);
    
    fprintf(ofp, "%10d (%6.2f%%) barcodes mapped to blacklisted targets\n",
           met->bl_mapped, ((double)(met->bl_mapped)/(double)(met->reads_processed*2)) * 100);
    
    fprintf(ofp, "%10d (%6.2f%%) barcodes did not map to a target\n",
           met->unmapped, ((double)(met->unmapped)/(double)(met->reads_processed*2)) * 100);
    
    //close mapping metrics file
    if (fclose(ofp) == EOF) {
        printf("print_DNA_qc_metrics: error - failed to close mapping metrics file. aborting...\n");
        abort();
    }
    
    //print mapping metrics to screen
    printf("\n%10d (%6.2f%%) barcodes mapped to whitelisted targets\n",
           met->mapped, ((double)(met->mapped)/(double)(met->reads_processed*2)) * 100);
    
    printf("%10d (%6.2f%%) barcodes mapped to blacklisted targets\n",
           met->bl_mapped, ((double)(met->bl_mapped)/(double)(met->reads_processed*2)) * 100);
    
    printf("%10d (%6.2f%%) barcodes did not map to a target\n\n",
           met->unmapped, ((double)(met->unmapped)/(double)(met->reads_processed*2)) * 100);
    
    return;
}


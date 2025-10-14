//
//  printQC_prcsMLT.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/gen_utils.h"
#include "../../../seq_utils/mapping_metrics.h"

#include "../MLT/testdata3pEnd_analysis.h"
#include "../MUX/testdataMUX_analysis.h"

#include "print_splitting_metrics.h"

/* print_metrics: prints a report of channel and 3' end processing metrics */
void print_splitting_metrics(TPROBE_names * nm, mapping_metrics  * met, fastp_params fastp_prms)
{
    extern struct testdata_3pEnd_vars testdata_3pEnd; //structure containing test data read analysis variables
    extern struct testdata_MUX_vars testdata_MUX;     //structure containing test data read analysis variables
    
    //(MAX_LINE*3)+512 should exceed the max string printed to out_str
    char out_str[(MAX_LINE*3)+512] = {0};
    
    FILE * log_fp = NULL; //output file pointer for processing log
    
    //open output files
    if ((log_fp = fopen("processing.txt", "w")) == NULL) {
        printf("print_metrics: ERROR - could not generate processing log file. Aborting program...\n");
        abort();
    }
    
    printf("\n");
    
    //print input file paths
    sprintf(out_str, "input files\nread1: %s\nread2: %s\ntrgts: %s\n", nm->file[READ1], nm->file[READ2], nm->trgts);
    printf2_scrn_n_fl(log_fp, out_str);
    
    sprintf(out_str, "\nsample names\nread1: %s\nread2: %s\n", nm->smpl[READ1], nm->smpl[READ2]);
    printf2_scrn_n_fl(log_fp, out_str);
    
    //print number of reads processed
    sprintf(out_str, "\n%d reads were processed\n\n", met->reads_processed);
    printf2_scrn_n_fl(log_fp, out_str);
    
    //printf channel splitting metrics
    sprintf(out_str,"channels:\n  total assessed: %d reads\n  untreated: %6.2f%% (%9d reads)\n  modified:  %6.2f%% (%9d reads)\n  unmapped:  %6.2f%% (%9d reads)\n",
            met->reads_processed,
            ((float)met->chan_count[UNT]/(float)(met->reads_processed)) *100, met->chan_count[UNT],
            ((float)met->chan_count[MOD]/(float)(met->reads_processed)) *100, met->chan_count[MOD],
            ((float)met->chan_count[ERR]/(float)(met->reads_processed)) *100, met->chan_count[ERR]);
    printf2_scrn_n_fl(log_fp, out_str);
    
    //print untreated channel barcode mapping metrics
    sprintf(out_str, "    of %d untreated reads,\n     %9d (%6.2f%%) were identified with a full barcode match and\n     %9d (%6.2f%%) with a partial (%d/%d) barcode match\n",
            met->chan_count[UNT],
            met->full_match[UNT], ((float)(met->full_match[UNT])/(float)(met->chan_count[UNT]))*100,
            met->part_match[UNT], ((float)(met->part_match[UNT])/(float)(met->chan_count[UNT]))*100,
            VL_MIN_MATCH, VL_MAX_MATCH);
    printf2_scrn_n_fl(log_fp, out_str);
    
    //print treated channel barcode mapping metrics
    sprintf(out_str, "    of %d modified reads,\n     %9d (%6.2f%%) were identified with a full barcode match and\n     %9d (%6.2f%%) with a partial (%d/%d) barcode match\n\n",
            met->chan_count[MOD],
            met->full_match[MOD], ((float)(met->full_match[MOD])/(float)(met->chan_count[MOD]))*100,
            met->part_match[MOD], ((float)(met->part_match[MOD])/(float)(met->chan_count[MOD]))*100,
            VL_MIN_MATCH, VL_MAX_MATCH);
    printf2_scrn_n_fl(log_fp, out_str);
    
    
    if (fastp_prms.mode == MULTI) {
        //3' end mapping metrics for multi-length cotranscriptional structure probing experiments
        
        //3' end mapping efficiency metrics
        sprintf(out_str, "3' end:\n  total assessed: %d reads\n  mapped:   %6.2f%% (%9d reads)\n  unmapped: %6.2f%% (%9d reads)\n",
                met->reads_processed,
                ((float)(met->mapped)   / (float)(met->reads_processed))*100, met->mapped,
                ((float)(met->unmapped) / (float)(met->reads_processed))*100, met->unmapped);
        printf2_scrn_n_fl(log_fp, out_str);
        
        //test data analysis-specific error indicating that <100% of test data reads were mappable
        //this is also reported from the print_3pEnd_testdata_analysis function, however this failure is
        //so severe that it is worth printing more than once
        if (testdata_3pEnd.run) {
            if ((met->reads_processed - met->mapped) != 0) {
                sprintf(out_str, "TEST DATA ANALYSIS: error - expected 100%% of reads to map correctly.\n");
            }
            printf2_scrn_n_fl(log_fp, out_str);
        }
        
        //3' end mapping quality metrics
        sprintf(out_str, "    of %d mapped reads,\n     %9d ends (%6.2f%%) mapped without substitutions or indels,\n     %9d ends (%6.2f%%) mapped without indels, and\n     %9d ends (%6.2f%%) match the expected sequence\n\n",
                met->mapped,
                met->nat_cnt, ((float)(met->nat_cnt)/(float)(met->mapped))*100,
                met->nat_cnt+met->sub_cnt, ((float)(met->nat_cnt+met->sub_cnt)/(float)(met->mapped))*100,
                met->hits, ((float)(met->matches)/(float)(met->hits))*100);
        printf2_scrn_n_fl(log_fp, out_str);
        
    } else if (fastp_prms.mode == MULTIPLEX) {
        //barcode mapping metrics for multi-length cotranscriptional structure probing experiments
        
        //barcode mapping efficiency metrics
        sprintf(out_str, "variant barcodes:\n  total assessed: %d reads\n  mapped:   %6.2f%% (%9d reads)\n  unmapped: %6.2f%% (%9d reads)\n",
                met->reads_processed,
                ((float)(met->mapped) / (float)(met->reads_processed))*100, met->mapped,
                ((float)(met->unmapped) / (float)(met->reads_processed))*100, met->unmapped);
        printf2_scrn_n_fl(log_fp, out_str);
        
        //test data analysis-specific error indicating that <100% of test data reads were mappable
        //this is also reported from the print_MUX_testdata_analysis function, however this failure is
        //so severe that it is worth printing more than once
        if (testdata_MUX.run) {
            if ((met->reads_processed - met->mapped) != 0) {
                sprintf(out_str, "TEST DATA ANALYSIS: error - expected 100%% of reads to map correctly.\n");
                printf2_scrn_n_fl(log_fp, out_str);
            }
        }
        
        //barcode mapping quality metrics
        sprintf(out_str, "    of %d mapped reads,\n     %9d ends (%6.2f%%) mapped without substitutions or indels,\n     %9d ends (%6.2f%%) mapped without indels, and\n     %9d ends (%6.2f%%) match the expected sequence\n\n",
                met->mapped,
                met->nat_cnt, ((float)(met->nat_cnt)/(float)(met->mapped))*100,
                met->nat_cnt+met->sub_cnt, ((float)(met->nat_cnt+met->sub_cnt)/(float)(met->mapped))*100,
                met->hits, ((float)(met->matches)/(float)(met->hits))*100);
        printf2_scrn_n_fl(log_fp, out_str);
    }
    
    
    int paranoid = 1;
    if (paranoid) {
        //sanity check that read 1 sequences are unchanged
        sprintf(out_str, "%d of %d read1 sequences (%.2f%%) match the expected sequence after processing\n",
                met->read_matches[READ1], met->reads_processed,
                ((float)(met->read_matches[READ1])/(float)(met->reads_processed))*100);
        printf2_scrn_n_fl(log_fp, out_str);
        
        //sanity check that read 2 sequences are unchanged
        sprintf(out_str, "%d of %d read2 sequences (%.2f%%) match the expected sequence after processing\n\n",
                met->read_matches[READ2], met->reads_processed,
                ((float)(met->read_matches[READ2])/(float)(met->reads_processed))*100);
        printf2_scrn_n_fl(log_fp, out_str);
        
        if ((fclose(log_fp)) == EOF) {
            printf("print_metrics: error - error occurred when closing metrics log file. Aborting program...\n");
            abort();
        }
    }
}



/* print_len_dist: print the observed 3' end distribution */
void print_len_dist(mapping_metrics *met, target3p_params trg_prms)
{
    extern struct testdata_3pEnd_vars testdata_3pEnd;    //structure containing test data read analysis variables
    
    int i = 0;
    
    char out_str[8192]; //array for storing output lines
    
    //open length distribution file
    FILE * dist_fp = NULL;
    if ((dist_fp = fopen("length_distribution.txt", "w")) == NULL) {
        printf("make_end_files: ERROR - could not generate length distribution file. Aborting program...\n");
        abort();
    }
    
    printf("\n");
    //print headings
    sprintf(out_str, "transcript length distribution:\n");
    printf2_scrn_n_fl(dist_fp, out_str);
    
    //if running test data read analysis, print length distribution column
    //that reports the number of reads in which the mapped 3' end value matched
    //the expected 3' end value exactly for each transcript length. this section
    //of the code prints a header for this additional 'distance zero' (dist0) column
    if (testdata_3pEnd.run) {
        sprintf(out_str, "length\tmapped\tdist0\n");
    } else {
        sprintf(out_str, "length\tmapped\n");
    }
    printf2_scrn_n_fl(dist_fp, out_str);
    
    
    //print length distribution
    for (i = trg_prms.min; i <= trg_prms.max; i++) {
        if (testdata_3pEnd.run) {    //if running test dat analysis, print distance zero column
            sprintf(out_str, "%d\t%d\t%d\n", i, met->len_dist[i], testdata_3pEnd.perfect[i]);
        } else {
            sprintf(out_str, "%d\t%d\n", i, met->len_dist[i]);
        }
        printf2_scrn_n_fl(dist_fp, out_str);
    }
    
    if ((fclose(dist_fp)) == EOF) {
        printf("print_len_dist: error - error occurred when closing length distribution file. Aborting program...\n");
        abort();
    }
}

//
//  mapping_metrics.h
//  
//
//  Created by Eric Strobel on 10/14/25.
//

#ifndef mapping_metrics_h
#define mapping_metrics_h

#include <stdio.h>
#include <stdlib.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

/* metrics: structure containing seq read processing metrics */
typedef struct mapping_metrics {
    int srcTrgs;                 //number of source barcode targets
    int targets;                 //number of barcode targets, including mutated barcodes
    int blacklisted;             //number of blacklisted targets
    int reads_processed;         //number of processed reads
    int * chan_count;            //count of reads from each channel
    int * full_match;            //counts number of full channel matches
    int * part_match;            //counts number of partial channel matches
    int hits;                    //tracks number of reads that were mapped by hash table
    int matches;                 //tracks number of hits with expected sequence (should be 100%)
    int nat_cnt;                 //tracks number of reads that match native targets
    int sub_cnt;                 //tracks number of reads that match substitution targets
    int ins_cnt;                 //tracks number of reads that match insertion targets
    int del_cnt;                 //tracks number of reads that match deletion targets
    int mapped;                  //mapped reads count (in doubly-barcoded libraries, tracks BC1 mapping)
    int bc2_mapped;              //barcode 2 mapped reads count
    int dblBC_mapped;            //number of doubly-barcoded reads in which both barcodes mapped
    int noBC_mapped;             //neither barcode mapped
    int bl_mapped;               //blacklisted mapped reads count
    int unmapped;                //unmapped reads count
    int concordant;              //tracks concordance when a second barcode is present
    int read_matches[READ_MAX];  //count of post-processing verified reads
    int len_dist[MAX_LINE];      //array for tracking the length distribution of UNT, MOD, and ERR reads
} mapping_metrics;

/* init_chnl_mtrcs_mem: initialize channel tracking metrics memory */
void init_chnl_mtrcs_mem(mapping_metrics * met, int channels);

#endif /* mapping_metrics_h */

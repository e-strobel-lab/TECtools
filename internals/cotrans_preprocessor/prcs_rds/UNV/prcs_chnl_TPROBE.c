//
//  prcs_chnl_TPROBE.c
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

#include "../../../seq_utils/mapping_metrics.h"

#include "../../../variant_maker/make_barcodes.h"

#include "prcs_chnl_TPROBE.h"

/* prcs_chnl_TPROBE: identify channel code from read 1 ID.
 these barcodes are used to split reads into modified and untreated channels
 barcodes:
 modified  = RRRYY
 untreated = YYYRR
 */
int prcs_chnl_TPROBE(char * read1_ID, mapping_metrics  * met, int mode)
{
    extern int debug;    //flag to turn on debug mode
    
    int i = 0;
    int j = 0;
    
    int UNT_mtch = 0;    //number of matches to UNT channel barcode
    int MOD_mtch = 0;    //number of matches to MOD channel barcode
    
    char barcode[VL_CHNL_BC_LEN+1] = {0}; //array to start channel barcode sequence
    
    int offset; //number of characters from the first space in the read id to the first base of the channel barcode
    
    if (mode == MULTI) {            //in PROCESS_MULTI mode, use TECprobe-VL UMI length
        offset = VL_UMI_LENGTH;
        
    } else if (mode == MULTIPLEX) { //in PROCESS_MULTIPLEX mode, use maximum barcode length
        offset = MAX_BARCODE_LEN;
        
    } else {
        printf("prcs_chnl_TPROBE: error -unrecognized mode. aborting...\n");
        abort();
    }
    
    //for compatibility with cotranscriptional structure probing, we pull the channel barcode as a
    //9 nt or 16 nt UMI even though the channel barcode is only 5 nts. The remaining 4/11 nt are constant sequence
    //and are ignored in the processing below
    //
    //TECprobe-VL format:               vvvvv-channel barcode
    //[illumina read id info]_NNNNNNNNN_NNNNNATGG [illumina read id info]
    //                        |R1 head| |R2 head|
    //
    //TECprobe-MUX format:                      vvvvv-channel barcode
    //[illumina read id info]_NNNNNNNNNNNNNNNN_NNNNNATGGCCTTCGG [illumina read id info]
    //                        |---R1 head----| |----R2 head---|
    
    if (debug) {
        printf(">channel determination\n");
        printf("R1 id:\t%s\n", read1_ID);
    }
    
    //iterate to first space in read 1 name, which indicates end of channel barcode
    for (i = 0; read1_ID[i] != ' ' && read1_ID[i]; i++) { ;}
    
    //check that read1_ID contained space and that expected UMI location is within array bounds
    if (!read1_ID[i]) { //check that loop did not exit on null character
        printf("prcs_chnl_TPROBE: error - unexpected read1 id line format. aborting...\n");
        abort();
    } else if (i <= offset) {//check that 1-offset will not be negative
        printf("prcs_chnl_TPROBE: error - unexpected short read1 id line. aborting...\n");
        abort();
    }
    
    i -= offset; //start of channel barcode is <offset> chars before first space in read 1 name
    
    //copy channel barcode to barcode array
    for (j = 0; j < VL_CHNL_BC_LEN; j++, i++) {
        barcode[j] = read1_ID[i];
    }
    barcode[j] = '\0';
    
    if (debug) {printf("barc:\t%s\nbarc:\t", barcode);}
    
    //determine whether channel barcode bases match MOD or UNT channels
    //UNT = YYYRR
    //MOD = RRRYY
    for (i = 0; i < VL_CHNL_BC_LEN && barcode[i]; i++) {
        switch (barcode[i]) {
            case 'A':
                (i <= 2) ? MOD_mtch++ : UNT_mtch++;
                if (debug) {printf("R");}
                break;
            case 'G':
                (i <= 2) ? MOD_mtch++ : UNT_mtch++;
                if (debug) {printf("R");}
                break;
            case 'T':
                (i <= 2) ? UNT_mtch++ : MOD_mtch++;
                if (debug) {printf("Y");}
                break;
            case 'C':
                (i <= 2) ? UNT_mtch++ : MOD_mtch++;
                if (debug) {printf("Y");}
                break;
            case 'N':
                if (debug) {printf("N");}
                break;
            default:
                printf("prcs_chnl_TPROBE: error - unexpected character %c(%d) in channel barcode. aborting...\n", barcode[i], barcode[i]);
                abort();
                break;
        }
    }
    if (debug) {printf("\nuntBC:\t%d\nmodBC:\t%d\n",UNT_mtch, MOD_mtch);}
    
    //check if the observed barcode is a match to an expected barcode
    //VL_MIN_MATCH is currently set to 4 (4/5 bases must match expected barcode)
    if (UNT_mtch >= VL_MIN_MATCH) {
        (UNT_mtch == VL_MAX_MATCH) ? met->full_match[UNT]++ : met->part_match[UNT]++;
        if (debug) {printf("chan:\tuntreated\n\n");}
        return UNT;
    } else if (MOD_mtch >= VL_MIN_MATCH) {
        (MOD_mtch == VL_MAX_MATCH) ? met->full_match[MOD]++ : met->part_match[MOD]++;
        if (debug) {printf("chan:\tmodified\n\n");}
        return MOD;
    } else {
        if (debug) {printf("chan:\tundetermined\n\n");}
        return ERR;
    }
}

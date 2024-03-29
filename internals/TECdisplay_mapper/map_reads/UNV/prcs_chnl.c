//
//  prcs_chnl.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"

#include "prcs_chnl.h"

/* prcs_chnl: identify channel code from read 1 ID.
 these barcodes are used to split reads into modified and untreated channels
 barcodes:
 bound  = RYYY
 unbound = YRRR
 
 in TECdisplay experiments the channel barcode is the first four nucleotides of read2
 */
int prcs_chnl(char * read_ID, metrics  * met, int * chnl_match_type)
{
    extern int debug;    //flag to turn on debug mode
    
    int i = 0;           //general purpose index
    int j = 0;           //general purpose index
    
    int BND_mtch = 0;    //number of matches to BND channel barcode
    int UNB_mtch = 0;    //number of matches to UNB channel barcode
    
    char barcode[CHANNEL_BC_LENGTH+1] = {0}; //array to start channel barcode sequence
    
    if (debug) {
        printf(">channel determination\n");
        printf("R1 id:\t%s\n", read_ID);
    }
    
    //iterate to first space in the merged read id, which indicates the end of the channel barcode
    for (i = 0; read_ID[i] != ' ' && read_ID[i]; i++) { ;}
    
    //check that the read_ID contained a space and that expected UMI location is within array bounds
    if (!read_ID[i]) { //check that loop did not exit on null character
        printf("prcs_chnl: error - unexpected read id line format. aborting...\n");
        abort();
    } else if (i <= UMI_LENGTH) {//check that 1-UMI_Length will not be negative
        printf("prcs_chnl: error - unexpected short read id line. aborting...\n");
        abort();
    }
    
    i -= UMI_LENGTH; //start of channel barcode is UMI_LENGTH chars before first space in read 1 name
    
    //copy channel barcode to barcode array
    for (j = 0; j < CHANNEL_BC_LENGTH; j++, i++) {
        barcode[j] = read_ID[i];
    }
    barcode[j] = '\0';
    
    if (debug) {printf("barc:\t%s\nbarc:\t", barcode);}
    
    //determine whether channel barcode bases match MOD or UNT channels
    //BND = RYYY
    //UNB = YRRR
    for (i = 0; i < CHANNEL_BC_LENGTH && barcode[i]; i++) {
        switch (barcode[i]) {
            case 'A':
                (i < 1) ? BND_mtch++ : UNB_mtch++;
                if (debug) {printf("R");}
                break;
            case 'G':
                (i < 1) ? BND_mtch++ : UNB_mtch++;
                if (debug) {printf("R");}
                break;
            case 'T':
                (i < 1) ? UNB_mtch++ : BND_mtch++;
                if (debug) {printf("Y");}
                break;
            case 'C':
                (i < 1) ? UNB_mtch++ : BND_mtch++;
                if (debug) {printf("Y");}
                break;
            case 'N':
                if (debug) {printf("N");}
                break;
            default:
                printf("prcs_chnl: error - unexpected character %c(%d) in channel barcode. aborting...\n", barcode[i], barcode[i]);
                abort();
                break;
        }
    }
    if (debug) {printf("\nbndBC:\t%d\nunbBC:\t%d\n",BND_mtch, UNB_mtch);}
    
    //check if the observed barcode is a match to an expected barcode
    //MIN_MATCH is currently set to 3 (3/4 bases must match expected barcode)
    if (UNB_mtch >= MIN_MATCH) {          //read is an unbound channel match
        if (UNB_mtch == MAX_MATCH) {      //read is a full match
            met->full_match[UNB]++;       //increment unbound full match counter
            *chnl_match_type = FULL;      //set channel match type to FULL
        } else {                          //read is a part match
            met->part_match[UNB]++;       //increment unbound part match counter
            *chnl_match_type = PART;      //set channel match type to PART
        }
        
        if (debug) {printf("chan:\tunbound\n\n---------------\n\n");}
        return UNB;
        
    } else if (BND_mtch >= MIN_MATCH) {   //read is a bound channel match
        if (BND_mtch == MAX_MATCH) {      //read is a full match
            met->full_match[BND]++;       //increment bound full match counter
            *chnl_match_type = FULL;      //set channel match type to FULL
        } else {                          //read is a part match
            met->part_match[BND]++;       //increment bound part match counter
            *chnl_match_type = PART;      //set channel match type to PART
        }
        
        if (debug) {printf("chan:\tbound\n\n---------------\n\n");}
        return BND;
        
    } else {                              //read did not match a channel
        *chnl_match_type = NO_CHNL_MATCH; //set channel match type to NO_CHNL_MATCH
        if (debug) {printf("chan:\tundetermined\n\n---------------\n\n");}
        return ERR;
    }
}

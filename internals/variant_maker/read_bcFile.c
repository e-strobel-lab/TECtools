//
//  read_bcFile.c
//  
//
//  Created by Eric Strobel on 8/28/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../utils/io_management.h"
#include "../utils/gen_utils.h"
#include "../seq_utils/isDNAbase.h"

//#include "./variant_maker_defs.h"
//#include "./variant_maker_structs.h"

#include "read_bcFile.h"

/* read_bcFile: parse barcode file header and barcode lines */
int read_bcFile(FILE * fp_brcd, int mode, char * str, int array_len)
{
    int i = 0;      //general purpose index
    int j = 0;      //general purpose index
    int bc_cnt = 0; //barcode count

    char line[MAX_LINE+1] = {0};    //line storage
    char tmp_str[MAX_LINE+1] = {0}; //temporary string storage
    
    const char lnkr_tag[8] = "linker=";   //tag indicating line is linker line
    const char strct_tag[9] = "bStruct="; //tag indicating line is barcode secondary structure line
    
    char * p_str = NULL; //pointer to a string
    
    if (mode == HEADER) { //reading header lines
        
        if (!get_line(line, fp_brcd)) { //get first header line
            return EOF;
        }
                                      
        for (i = 0; isdigit(line[i]) && i < MAX_LINE; i++) {  //store bc count string
            tmp_str[i] = line[i];
        }
        tmp_str[i] = '\0';
        
        if (!tmp_str[0]) { //throw error if barcode file doesn't start with a digit
            printf("read_bcFile: ERROR - expected file to start with the number of barcodes. aborting...\n");
            abort();
        }
        bc_cnt = atoi(tmp_str); //set barcode count
        
        
        get_line(line, fp_brcd);                           //get second header line
        if (!memcmp(line, strct_tag, strlen(strct_tag))) { //check line for linker tag
            p_str = &line[strlen(strct_tag)];
            strcpy(str, p_str);
        } else {
            printf("read_bcFile: ERROR - unexpected format for barcodes file. Could not find barcode structure line. aborting...\n");
            abort();
        }
        
        return bc_cnt; //return barcode count
        
    } else if (mode == BC_LINE) { //reading barcode lines
        
        if (!get_line(line, fp_brcd)) { //get barcode line
            return EOF;                 //no more barcodes to get
            
        } else {
            for (i = 0; isdigit(line[i]) && i < MAX_LINE; i++) { } //iterate past bc index
            if (line[i] == '\t') {  //check that next char is a tab
                line[i] = '\0';     //terminate bc index string
                p_str = &line[i+1]; //set pointer to barcode sequence
                //TODO: need array length check to ensure str is long enough
                strcpy(str, p_str); //store barcode sequence in str
            } else {
                printf("read_bcFile: ERROR - unexpected format for barcodes file. expected barcode index and sequence to be tab-separated. aborting...\n");
                abort();
            }
            return atoi(line); //return barcode index
        }
        
    } else {
        printf("read_bcFile: ERROR - unexpected mode. aborting...\n");
        abort();
    }
}

//
//  parse_brcd_id_list.c
//  
//
//  Created by Eric Strobel on 9/24/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/io_management.h"

#include "../UNV/config_struct.h"

#include "parse_brcd_id_list.h"

/* parse_brcd_id_list: parse barcode ids file and store barcode ids */
int parse_brcd_id_list(tprobe_configuration * config, char *** brcd_id)
{
    const char brcd_id_sffx[17] = "_barcode_ids.txt"; //barcode id file suffix
    
    char brcd_id_loc[MAX_LINE+1] = {0}; //barcode id file location
    int  brcd_id_cnt = 0;               //number of barcode ids
    
    int i = 0;   //general purpose index
    int j = 0;   //general purpose index
    int ret = 0; //return value for snprintf
    
    char line[MAX_LINE+1] = {0};    //storage for barcode id file lines
    char tmp_str[MAX_LINE+1] = {0}; //temp storage for barcode id strings
    
    //construct barcode id file location
    ret = snprintf(brcd_id_loc, MAX_LINE, "%s%s%s", config->trg_files_loc, config->trg_files_prfx, brcd_id_sffx);
    if (ret >= MAX_LINE || ret < 0) {
        printf("parse_brcd_id_list: error - error when constructing barcode id file location. aborting...\n");
        abort();
    }
    
    FILE * fp_bid = NULL; //barcode id file pointer
    
    //determine number of barcode ids and allocate memory for barcode strings
    if ((fp_bid = fopen(brcd_id_loc, "r")) == NULL) { //open barcode id file
        printf("parse_brcd_id_list: error - failed to open barcode id list. aborting...\n");
        abort();
    }
    
    for (i = 0; get_line(line, fp_bid); i++) {;} //count number of barcode ids
    brcd_id_cnt = i;                             //store number of barcode ids
    
    if (((*brcd_id) = calloc(brcd_id_cnt, sizeof(**brcd_id))) == NULL) { //allocate memory for barcode ids
        printf("parse_brcd_id_list: error - failed to allocate memory for barcode ids. aborting...\n");
        abort();
    }
    
    if (fclose(fp_bid) == EOF) { //close barcode id file
        printf("parse_brcd_id_list: error - failed to close barcode id list. aborting...\n");
    }
    
    //open barcode ids file again. this time, allocate memory for and store barcode ids
    if ((fp_bid = fopen(brcd_id_loc, "r")) == NULL) {
        printf("parse_brcd_id_list: error - failed to open barcode id list. aborting...\n");
        abort();
    }
    
    //read, allocate memory for, and store barcode ids
    for (i = 0; i < brcd_id_cnt; i++) { //for every barcode id
        get_line(line, fp_bid);         //get barcode id line
        
        //read barcode id line and check that all chars in barcode are digits
        for (j = 0; line[j]; j++) {     //for every char in barcode id
            if (isdigit(line[j])) {     //if the char is a digit
                tmp_str[j] = line[j];   //store the digit in tmp_str
            } else {                    //otherwise, throw error and abort
                printf("parse_brcd_id_list: error - non-digit character detected in barcode id. aborting...\n");
                abort();
            }
        }
        tmp_str[j] = '\0'; //terminate string
        
        //allocate memory for barcode id
        if (((*brcd_id)[i] = malloc(strlen(tmp_str+1) * sizeof(*((*brcd_id)[i])))) == NULL) {
            printf("parse_brcd_id_list: error - failed to allocate memory for barcode ids. aborting...\n");
            abort();
        }
        strcpy((*brcd_id)[i], tmp_str); //store barcode id
    }
    
    if (fclose(fp_bid) == EOF) { //close barcode id file
        printf("parse_brcd_id_list: error - failed to close barcode id list. aborting...\n");
    }
    
    return brcd_id_cnt;
}

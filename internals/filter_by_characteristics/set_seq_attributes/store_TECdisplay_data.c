//
//  store_TECdisplay_data.c
//  
//
//  Created by Eric Strobel on 3/5/26.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"
#include "../../utils/gen_utils.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"

#include "../compare_variants/parse_comparison_file.h"

#include "store_TECdisplay_data.h"

/* store_TECdisplay_data: store input TECdisplay data in the sequence attributes structure */
void store_TECdisplay_data(sequence_attributes * sq_att, char ** tecd_data_hdr, char *** prsd_tecd_data_hdr, int seq_cnt, char * tecd_path, comparison_values * cmp)
{
    FILE * p_tecd = NULL; //pointer to TECdisplay data file
    
    int i = 0;                   //general purpose index
    int j = 0;                   //general purpose index
    char line[MAX_LINE+1] = {0}; //array to store lines from TECdisplay data file
    char * p_hdr = NULL;         //pointer to header string
    char * p_nm = NULL;          //pointer to sequence name string
    char * p_data = NULL;        //pointer to data string
    int reached_eof = 0;         //flag that end of file was reached
    
    int fields = 0;              //number of fields in data header
    
    get_file(&p_tecd, tecd_path); //open TECdisplay data file
    
    //get column headear string. the first column header (for variant ids)
    //is skipped since this is printed to the output file separately
    get_line(line, p_tecd);               //read header line
    for (i = 0; line[i] != '\t'; i++) {;} //iterate to first tab
    line[i] = '\0';                       //truncate string at first tab
    p_hdr = &line[i+1];                   //set pointer to char after the first tab
    if ((*tecd_data_hdr = malloc((strlen(p_hdr)+1) * sizeof(**tecd_data_hdr))) == NULL ) { //allocate memory for TECdisplay data header
        printf("store_TECdisplay_data: failed to allocate memory for TECdisplay data header storage. aborting...");
        abort();
    }
    strcpy(*tecd_data_hdr, p_hdr); //store TECdisplay data header
    
    //clear trailing space characters from TECdisplay data file header
    clear_trailing_spaces(*tecd_data_hdr);
    
    fields = parse_tecd_data_hdr(prsd_tecd_data_hdr, *tecd_data_hdr); //parse data header and store field names
    find_cmp_index(*prsd_tecd_data_hdr, fields, cmp);                 //store the column index to be used for comparison
    
    for (i = 0, reached_eof = 0; i < seq_cnt && !reached_eof; i++) { //process all data lines
        if (get_line(line, p_tecd)) {                                //if a data line was found
            for (j = 0; line[j] != '\t'; j++) {;}                    //iterate to the first tab
            line[j] = '\0';                                          //replace the first tab with a null character
            p_nm = &line[0];                                         //set the sequence name pointer to index zero
            p_data = &line[j+1];                                     //set the data pointer to one character after the first tab
            clear_trailing_spaces(p_data);                           //clear trailing spaces from the data string
            
            if (!strcmp(p_nm, sq_att[i].nm)) { //check that the variant id from the TECdisplay data line matches the current variant id
                
                //allocate memory for storing the TECdisplay data string
                if ((sq_att[i].tecd = malloc((strlen(p_data)+1) * sizeof(*(sq_att[i].tecd)))) == NULL ) {
                    printf("store_TECdisplay_data: failed to allocate memory for TECdisplay data storage. aborting...");
                    abort();
                }
                strcpy(sq_att[i].tecd, p_data); //store the TECdisplay data string
                                
            } else { //error - msa and TECdisplay data file are discordant
                printf("store_TECdisplay_data: variant ids in msa and TECdisplay data inputs are discordant. aborting...\n");
                abort();
            }

            parse_tecd_data_line(&sq_att[i].td_vals, sq_att[i].tecd, fields); //parse data line and store values
            
        } else {
            reached_eof = 1; //set flag that eof was reached
        }
    }
    
    if (i == seq_cnt && get_line(line, p_tecd)) { //check that there are no more data lines in the TECdisplay data file
        printf("store_TECdisplay_data: error - more data lines in TECdisplay data file than variants in multiple sequence alignment. aborting...\n");
        abort();
        
    } else if (i < seq_cnt) { //check that there was a TECdisplay data line for every variant in the msa
        printf("store_TECdisplay_data: error - too few data lines in TECdisplay data file. aborting...\n");
        abort();
    }
            
    return;
}

/* parse_tecd_data_hdr: parse TECdisplay data header string and store individual fields */
int parse_tecd_data_hdr(char *** prsd_tecd_data_hdr, char * tecd_data_hdr)
{
    int i = 0;                 //general purpose index
    int j = 0;                 //general purpose index
    int k = 0;                 //general purpose index
    int fields = 0;            //number of fields
    int fnd_end = 0;           //flag that end of string was found
    char tmp[MAX_LINE] = {0};  //array for storing field names temporarily
    
    //count number of fields
    for (i = 0, fields = 0, fnd_end = 0; !fnd_end; i++) {
        if (tecd_data_hdr[i] == '\t' || tecd_data_hdr[i] == '\0') { //if a tab or null char was found
            fields++;                                               //increment field count
        }
        
        if (tecd_data_hdr[i] == '\0') { //if a null char was found
            fnd_end = 1;                //set flag to terminate the loop
        }
    }

    //allocate pointers for allocating memory to store field name strings
    if ((*prsd_tecd_data_hdr = malloc(fields * sizeof(**prsd_tecd_data_hdr))) == NULL) {
        printf("parse_tecd_data_hdr: error - failed to allocate memory for storing parsed TECdisplay data header.\n");
        abort();
    }
    
    //store field name strings
    for (i = 0, j = 0; i < fields; i++, j++) { //for every field
        
        //store field name in tmp string
        for (k = 0; tecd_data_hdr[j] != '\t' && tecd_data_hdr[j] != '\0'; j++, k++) {
            tmp[k] = tecd_data_hdr[j];
        }
        tmp[k] = '\0';
        
        //allocate memory for storing current field name
        if (((*prsd_tecd_data_hdr)[i] = malloc((strlen(tmp)+1) * sizeof(***prsd_tecd_data_hdr))) == NULL) {
            printf("parse_tecd_data_hdr: error - failed to allocate memory for storing parsed TECdisplay data header.\n");
            abort();
        }
        strcpy((*prsd_tecd_data_hdr)[i], tmp); //store field name
    }
    
    return fields; //return number of parsed fields
}

/* parse_tecd_data_line: parse TECdisplay data line and store values as doubles */
void parse_tecd_data_line(double ** vals, char * p_data, int fields)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
    
    char tmp[MAX_LINE+1] = {0}; //array for storing field value temporarily
    
    //allocate memory for storing field values
    if ((*vals = malloc(fields * sizeof(**vals))) == NULL) {
        printf("parse_tecd_data_line: error - failed to allocate memory for storing parsed TECdisplay data values.\n");
        abort();
    }
    
    //parse field substrings and store as doubles
    for (i = 0, j = 0; i < fields; i++, j++) { //for every field
        
        //parse substring and store in tmp
        for (k = 0; p_data[j] != '\t' && p_data[j] != '\0'; j++, k++) {
            tmp[k] = p_data[j];
        }
        tmp[k] = '\0';
        
        check_float_str(tmp, ABORT_FAILURE); //check that all chars in tmp are valid float characters
        (*vals)[i] = strtod(tmp, NULL);      //convert tmp to double and store value
    }
    
    return;
}

/* find_cmp_index: identify index of comparison column specified by input comparison file */
void find_cmp_index(char ** prsd_tecd_data_hdr, int fields, comparison_values * cmp)
{
    int i = 0;         //general purpose index
    int fnd_match = 0; //flag that match was found
    
    //find index of comparison column
    for (i = 0, fnd_match = 0; i < fields  && !fnd_match; i++) { //for every field
        if (!strcmp(prsd_tecd_data_hdr[i], cmp->nm)) {           //check if field name matches comparison column name
            cmp->ix = i;                                         //if so, set comparison index and fnd_match flag
            fnd_match = 1;
        }
    }
    
    if (!fnd_match) { //if no match was found, throw error and abort
        printf("find_cmp_index: error - failed to find match to comparison column name in TECdisplay data header. aborting...\n");
        abort();
    }
        
    return;
}

//
//  parse_comparison_file.c
//  
//
//  Created by Eric Strobel on 6/29/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"

#include "parse_comparison_file.h"

/* parse_comparison file: parse and store comparison details from input comparison values file */
int parse_comparison_file(char * cmp_path, comparison_values ** cmp)
{
    FILE * cfp = NULL; //pointer for comparison file
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int cmp_cnt;                          //number of comparison entries in comparison input file
    cmp_cnt = count_entries(cmp_path, 1); //count number of comparison entries
    
    char line[MAX_LINE+1] = {0}; //line storage
    
    char * p_nm = NULL;  //pointer to comparison entry name
    char * p_typ = NULL; //pointer to comparison entry type
    
    int fnd_tab = 0; //flag that tab was found
    
    char dif_str[4] = {"dif"}; //comparison type value: different
    char sim_str[4] = {"sim"}; //comparison type value: similar
    
    //allocate memory for descriptors
    if ((*cmp = calloc(cmp_cnt, sizeof(**cmp))) == NULL) {
        printf("parse_comparison_file: error - failed to allocate memory for comparison values. aborting...");
        abort();
    }
    
    //open comparison values file for parsing
    if ((cfp = fopen(cmp_path, "r")) == NULL) {
        printf("parse_comparison_file: error - failed to open comparison values file. aborting...\n");
        abort();
    }
    
    /* process comparison entries*/
    for (i = 0; get_line(line, cfp); i++) {
        
        //throw error and abort if leading whitespace is present
        if (isspace(line[0])) {
            printf("parse_comparison_file: error - unexpected leading whitespace in comparison entry. aborting...\n");
            abort();
        }
        
        p_nm = &line[0]; //set pointer to beginning of line
        p_typ = NULL;    //set p_typ to NULL
        
        //split line into name and type
        for (j = 0, fnd_tab = 0; line[j]; j++) { //for each character of the line
            if (line[j] == '\t') {       //if the character is a tab..
                if (fnd_tab) {           //if a tab was found previously, throw error and abort
                    printf("parse_comparison_file: error - detected multiple tabs in comparison entry line. aborting...\n");
                    abort();
                } else {                 //if this is the first tab
                    line[j] = '\0';      //change tab to a null terminator
                    p_typ = &line[j+1];  //set p_typ pointer to the next character in the array
                                         //TODO: consider checking if next character is alpha/digit
                    fnd_tab = 1;         //set flag that tab was found
                }
            }
        }
        
        if (fnd_tab && line[j] == '\0') { //if a tab was found and the whole line was tested
            
            //allocate memory for storing name in comparison_values struct
            if (((*cmp)[i].nm = malloc(strlen(p_nm) * sizeof(*((*cmp)[i].nm)))) == NULL) {
                printf("parse_comparison_file: error - failed to allocate memory for comparison values name\n");
                abort();
            }
            strcpy((*cmp)[i].nm, p_nm); //copy name to comparison values truct
            
            if (!strcmp(p_typ, sim_str)) { //if comparison type matches sim_str
                (*cmp)[i].typ = SIM;       //set type to similar
                
            } else if (!strcmp(p_typ, dif_str)) { //if comparison type matches dif_str
                (*cmp)[i].typ = DIF;              //set type to different
            
            } else { //unexpected type, throw error and abort
                printf("parse_comparison_file: error - unexpected comparison entry '%s'. aborting...\n", p_typ);
                abort();
            }
            
        } else { //unexpected line format, throw error and abort
            printf("parse_comparison_file: error - unexpected format for comparison entry line. aborting...\n");
            abort();
        }
    }
    
    //close comparison file
    if (fclose(cfp) == EOF) {
        printf("parse_comparison_file: error - failed to close comparison file. aborting...\n");
        abort();
    }
    
    return cmp_cnt;
}


/* count_entries: count number of lines in input comparison file */
int count_entries(char * path, int lpe)
{
    FILE * ipt = NULL; //input file pointer
    
    char line[MAX_LINE+1] = {0}; //line storage
    
    int l = 0;     //general purpose index
    int mod = 0;   //modulus
    
    //open comparison values file
    if ((ipt = fopen(path, "r")) == NULL) {
        printf("count_entries: error - failed to open input file. aborting...\n");
        abort();
    }
    
    //count number of lines
    for (l = 0; get_line(line, ipt); l++) {;}
    
    //close comparison values file
    if (fclose(ipt) == EOF) {
        printf("count_comparisons: error - failed to close comparison values file. aborting...\n");
        abort();
    }
    
    if ((mod = l % lpe)) { //if modulus is not 0, incomplete entry is present
        printf("count_entries: error - detected incomplete entry (%d/%d lines). aborting...\n", mod, lpe);
        abort();
    } else {
        return l/lpe; //return number of entries
    }
}

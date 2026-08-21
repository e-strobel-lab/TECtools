//
//  parse_comparison_file.c
//  
//
//  Created by Eric Strobel on 6/29/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <float.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"
#include "../../utils/gen_utils.h"

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
    
    char * p_des = NULL;   //pointer to descriptor name string
    char * p_val = NULL;   //pointer to comparison value string
    char * p_oth = NULL;   //pointer to other values string
    char * p_minZ = NULL;  //pointer to minimum z-score string
    char * p_maxdG = NULL; //pointer to maximum deltaG string
    
    int fnd_tab = 0; //flag that tab was found
    int fnd_end = 0; //flag that terminating null was found
    
    const char min_z_dif_str[11]  = {"min_z_dif="};  //leader string ahead of minimum z-score difference value
    const char max_dG_dif_str[12] = {"max_dG_dif="}; //leader string ahead of maximum delta G difference value
    
    char * p_lmt = NULL;     //pointer to limit value
    double * lmt2set = NULL; //pointer to limit value that will be set
 
    //allocate memory for descriptors
    if ((*cmp = calloc(cmp_cnt, sizeof(**cmp))) == NULL) {
        printf("parse_comparison_file: error - failed to allocate memory for comparison values. aborting...");
        abort();
    }
    (*cmp)->minZ  = LIMIT_INIT; //initialize minZ to min float value
    (*cmp)->maxdG = LIMIT_INIT; //initialize maxdG to min float value
    
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
        
        p_des = &line[0]; //point p_des to beginning of line
        p_val = NULL;     //set p_val to NULL
        p_oth = NULL;     //set p_oth to NULL
        p_minZ = NULL;    //set p_minZ to NULL
        p_maxdG = NULL;   //set p_maxdG to NULL

        //set pointer to value string
        for (j = 0; line[j] != '\t' && line[j]; j++) {;}
        if (line[j] == '\t') {
            line[j] = '\0';
            p_val = &line[j+1];
        } else {
            printf("parse_comparison_file: error - unexpected format for comparison file line. aborting...\n");
            abort();
        }
        
        //check if there are additional fields
        for (j++; line[j] != '\t' && line[j]; j++) {;}
        if (line[j] == '\t') {  //if there is another field
            line[j] = '\0';     //terminate the previous field string
            p_oth = &line[j+1]; //point p_oth to the start of the next field
            
            //parse remaining string
            fnd_end = 0; //initialize fnd_end
            while (!fnd_end) { //until the end of the string is found
                
                if (!memcmp(p_oth, min_z_dif_str, strlen(min_z_dif_str))) { //test if value is min_z_dif
                    lmt2set = &(*cmp)->minZ;                                //set minZ as limit to set
                    p_lmt = &p_oth[strlen(min_z_dif_str)];                  //point p_lmt to value
                    
                } else if (!memcmp(p_oth, max_dG_dif_str, strlen(max_dG_dif_str))) { //test if value is max_dG_dif
                    lmt2set = &(*cmp)->maxdG;                                        //set maxdG as limit to set
                    p_lmt = &p_oth[strlen(max_dG_dif_str)];                          //pointe p_lmt to value
                    
                } else {
                    printf("parse_comparison_file: error - unrecognized limit specifier. aborting...\n");
                    abort();
                }
                
                fnd_tab = 0; //initialize fnd_tab
                for (j = 0; p_lmt[j] != '\t' && p_lmt[j]; j++) {;} //iterate to tab or end of string
                if (p_lmt[j] == '\t') {                            //if tab was found
                    fnd_tab = 1;                                   //set fnd_tab to 1
                } else if (!p_lmt[j]) {                            //if end of string was found
                    fnd_end = 1;                                   //set fnd_end to 1
                }
                
                p_lmt[j] = '\0';                       //terminate limit value string
                check_float_str(p_lmt, ABORT_FAILURE); //check that string is comprised of valid float characters
                *lmt2set = strtod(p_lmt, NULL);        //store limit value
                
                if (fnd_tab) {           //if a tab was found
                    p_oth = &p_lmt[j+1]; //point p_oth to the start of the next field
                }
            }
        }
        
        //allocate memory for storing name in comparison_values struct
        if (((*cmp)[i].des = malloc(strlen(p_des) * sizeof(*((*cmp)[i].des)))) == NULL) {
            printf("parse_comparison_file: error - failed to allocate memory for comparison descriptor name\n");
            abort();
        }
        strcpy((*cmp)[i].des, p_des); //copy name to comparison values truct
        
        //allocate memory for storing name in comparison_values struct
        if (((*cmp)[i].val = malloc(strlen(p_val) * sizeof(*((*cmp)[i].val)))) == NULL) {
            printf("parse_comparison_file: error - failed to allocate memory for comparison values name\n");
            abort();
        }
        strcpy((*cmp)[i].val, p_val); //copy name to comparison values truct
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

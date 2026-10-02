//
//  parse_id_link.c
//  
//
//  Created by Eric Strobel on 9/25/26.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

#include "../../TECdisplay_mapper/TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper/TECdisplay_mapper_structs.h"

#include "../../utils/io_management.h"

#include "parse_id_link.h"

/* parse_id_link: parse id link file and store linked concise/verbose ids and barcode sequences in var_id_link structure */
void parse_id_link(var_id_link * id_link, int cnt, char * f_nm)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
    
    char line[MAX_LINE+1] = {0}; //line storage
    char * crnt_fld = NULL;      //pointer to start of current field within a line
    char ** fld2set = NULL;      //field to set in current iteration of line parsing
    int L0_flds = 0;             //number of fields in line zero, for comparison against other lines
    
    int fnd_end = 0; //flag that line end was found
    
    FILE * ifp = NULL; //input file pointer
    
    //open id_link file
    if ((ifp = fopen(f_nm, "r")) == NULL) {
        printf("parse_id_link: error - failed to open id link file. aborting...\n");
        abort();
    }
    
    for (i = 0, L0_flds = 0; get_line(line, ifp); i++) { //for every line in the file
        
        //parse each tab-delimited field of the line
        for (j = 0, crnt_fld = line, fnd_end = 0; !fnd_end && j < MAX_ID_LINK_FIELDS; j++) {
            
            //set the field that will be stored in the current iteration
            switch (j) {
                case 0: fld2set = &id_link[i].v;   break;
                case 1: fld2set = &id_link[i].c;   break;
                case 2: fld2set = &id_link[i].bc1; break;
                case 3: fld2set = &id_link[i].bc2; break;
                default:
                    printf("parse_id_link: error - >%d fields in id link file line. aborting...", MAX_ID_LINK_FIELDS);
                    abort();
                    break;
            }
            
            //iterate to the next tab or the end of the line
            for (k = 0; crnt_fld[k] != '\t' && crnt_fld[k] != '\n' && crnt_fld[k]; k++) {;}
            
            if (crnt_fld[k] == '\n' || crnt_fld[k] == '\0') { //if the end of the line was found
                fnd_end = 1;     //set flag that the end of the line was found
                if (i == 0) {    //if parsing first line of the file
                    L0_flds = j; //set number of fields
                    
                } else {                //otherwise
                    if (j != L0_flds) { //check that the number of fields matches that of the first line
                        printf("parse_id_link: error - number of fields in id link file lines is not uniform. aborting...\n");
                        abort();
                    }
                }
            }
            
            crnt_fld[k] = '\0'; //terminate the field string
            
            //allocate memory for storing the current field string
            if ((*fld2set = malloc((strlen(crnt_fld)+1) * sizeof(**fld2set))) == NULL) {
                printf("parse_id_link: error  - failed to allocate memory for id link string. aborting...\n");
                abort();
            }
            strcpy(*fld2set, crnt_fld); //store the field string
            crnt_fld = &crnt_fld[k+1];  //set crnt_fld to the next field
        }
    }
    
    if (i != cnt) { //if the number of lines in the id link file does not match the number of targets, abort
        printf("parse_id_link: error - number of entries in id link file does not match the number of targets in the fasta file. aborting...\n");
        abort();
    }
    
    if (fclose(ifp) == EOF) { //close id link file
        printf("parse_id_link: error - failed to close id link file. aborting...\n");
        abort();
    }
    
    return;
}

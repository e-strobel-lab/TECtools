//
//  get_constraint_metadata.c
//  
//
//  Created by Eric Strobel on 11/7/23.
//

#include <stdio.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"

#include "../seq_utils/basemap.h"
#include "../seq_utils/is_dgnrt_mtch.h"

#include "../TECdisplay_navigator/parse_reference.h"
#include "../TECdisplay_navigator/parse_constraints.h"
#include "../TECdisplay_navigator/read_vbase.h"
#include "../TECdisplay_navigator/search_4_vbase_match.h"

#include "./TECdisplay_Hnav_defs.h"
#include "./TECdisplay_Hnav_structs.h"
#include "./TECdisplay_Hnav_global_vars.h"

#include "get_constraint_metadata.h"

/* get_constraint_metadata: read constraint files to get information about individual constraints within
 the file. This process also validates each constraint before it is used by TECdisplay_navigator */
void get_constraint_metadata(char * ipt_fn, constraint_metadata * cons_meta, char type, int layr_cnt)
{
    char tmp_sn[MAX_LINE+1] = {0}; //temporary storage for sample name
    char * sn_ptr = NULL;          //pointer to sample name start
    
    get_file(&(cons_meta->fp), ipt_fn); //set file pointer to constraints file
    if (strlen(ipt_fn) < MAX_NAME) {
        strcpy(cons_meta->fn, ipt_fn);       //store file name
        get_sample_name(ipt_fn, &tmp_sn[1]); //get sample name from constraints file name, store at index 1
        if (type == 'c') {                   //if keeping matches,
            sn_ptr = &tmp_sn[1];             //set pointer to index 1 of tmp_sn
        } else if (type == 'x') {            //if excluding matches,
            tmp_sn[0] = 'x';                 //set index 0 of tmp_sn to 'x' to indicate exclusion
            sn_ptr = &tmp_sn[0];             //set pointer to index 0 of tmp_sn
        } else {
            printf("get_constraint_metadata: error - unrecognized constraint type. aborting...");
            abort();
        }
        strcpy(cons_meta->sn, sn_ptr); //copy sample name to cons_meta structure
        
    } else {
        printf("get_constraint_metadata: error - input file name exceeds the maximum length (%d chars). aborting...\n", MAX_NAME);
        abort();
    }
    cons_meta->typ = type;  //set constraint type
    
    FILE * tmp_ifp; //temporary file pointer for getting constraint names
    
    int i = 0; //general purpose index
    
    get_file(&(tmp_ifp), ipt_fn); //set temp file pointer to constraints file
    
    //parse constraints file for reference sequence header and individual constraints
    parse_reference(tmp_ifp, &cons_meta->bmap, &cons_meta->wt, &cons_meta->cnstnt_indels, !layr_cnt);
    cons_meta->c_cnt = parse_constraints(tmp_ifp, &cons_meta->con[0], &cons_meta->bmap, cons_meta->cnstnt_indels);
    
    //set pointers to individual constraint names
    for (i = 0; i < cons_meta->c_cnt; i++) {
        cons_meta->cn[i] = &cons_meta->con[i].nm[0];
    }
    
    /* close output file */
    if (fclose(tmp_ifp) == EOF) {
        printf("filter_values: error - error occurred when closing output file %s.txt. Aborting program...\n", ipt_fn);
        abort();
    }
}

/* valid8_constraint_compatiblity: confirm all constraints files are compatible */
void valid8_constraint_compatiblity(int layr_cnt, constraint_metadata * cons_meta)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int error = 0; //flag that an error was detected
    
    if (layr_cnt <= 1) { //one layer, no validation necessary
        return;
        
    } else {
        
        //check that the wt source seq, ref seq, and constant indels string
        //for each constraint file match those of the first constraint file
        
        printf("\nvalidating constraint compatibility...\n");
        
        for (i = 1, error = 0; i < layr_cnt; i++) {
            
            if (strcmp(cons_meta[i].wt.sq, cons_meta[0].wt.sq)) { //check wt source sequence
                error = 1;
                printf(">error - wt source for constraint file '%s' does not match the first constraint file (%s)\n", cons_meta[i].fn, cons_meta[0].fn);
            }
            
            if (strcmp(cons_meta[i].bmap.rS, cons_meta[0].bmap.rS)) { //check reference sequence
                error = 1;
                printf(">error - variable base reference sequence for constraint file '%s' does not match the first constraint file (%s)\n", cons_meta[i].fn, cons_meta[0].fn);
            }
            
            if (strcmp(cons_meta[i].cnstnt_indels, cons_meta[0].cnstnt_indels)) { //check constant indels string
                error = 1;
                printf(">error - constant indels string for constraint file '%s' does not match the for first constraint file (%s)\n", cons_meta[i].fn, cons_meta[0].fn);
            }
        }
    }
    
    if (error) { //if an error was detected, abort
        printf("aborting...\n");
        abort();
    } else {
        printf("...all checks passed\n\n");
    }
    
    return;
}

/* print_constraint_metadata: print metadata for each input constraint file to output file */
void print_constraint_metadata(int layr_cnt, constraint_metadata * cons_meta) {
    
    FILE * constraint_md_out = NULL; //pointer for constraint metadata output file
    if ((constraint_md_out = fopen("constraint_metadata.txt", "w")) == NULL) {
        printf("print_constraint_metadata: error - could not open constraint metadata output file. Aborting program...\n");
        abort();
    }
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    for (i = 0; i < layr_cnt; i++) { //print constraint metadata for each layer
        
        fprintf(constraint_md_out, "%s\n", cons_meta[i].fn);                  //print constraint filename
        fprintf(constraint_md_out, ">name:  %s\n", cons_meta[i].sn);          //print constraint sample name
        fprintf(constraint_md_out, ">type:  %c\n", cons_meta[i].typ);         //print constraint type
        fprintf(constraint_md_out, ">count: %d\n", cons_meta[i].c_cnt);       //print constraint count
        for (j = 0; j < cons_meta[i].c_cnt; j++) {                            //print constraint names
            fprintf(constraint_md_out, "  %2d: %s\n", j, cons_meta[i].cn[j]);
        }
        fprintf(constraint_md_out, "\n");
    }
    
    
    if (fclose(constraint_md_out) == EOF) { //close constraint metadata file
        printf("print_constraint_metadata: error - error occurred when closing constraint info output file. Aborting program...\n");
        abort();
    }
}

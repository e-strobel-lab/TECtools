//
//  process_output_files.c
//  
//
//  Created by Eric Strobel on 11/7/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/stat.h>

#include "../utils/io_management.h"

#include "./TECdisplay_Hnav_defs.h"
#include "./TECdisplay_Hnav_structs.h"
#include "./TECdisplay_Hnav_global_vars.h"

#include "process_output_files.h"

/* clean_up_output: clean up output directory after TECdisplay_navigator analysis is complete */
void clean_up_output(int layr_cnt, constraint_metadata * cons_meta)
{
    int i = 0; //general purpose index
    
    char rm_mrgd_out_command[MAX_LINE+1] = {0}; //array to store command for removing merged_out files
    char mv_outFiles_command[MAX_LINE+1] = {0}; //array to store command for moving TECdisplay_navigator output file
    char rm_out_dirs_command[MAX_LINE+1] = {0}; //array to store command for removing TECdisplay_navigator output directories
    
    char crnt_layr[16] = {0}; //array to store current layer name
    
    system("rm merged_out.txt"); //remove top level merged_out file
    
    for (i = 0; i < layr_cnt; i++) { //for each layer
        
        sprintf(crnt_layr, "layer%d", i+1); //generate the current layer name
        chdir(crnt_layr);                   //change to the current layer directory
        
        sprintf(rm_mrgd_out_command, "rm -r */*merged_out.txt"); //generate command to remove merged out files
        system(rm_mrgd_out_command);                             //run remove merged out command
        
        if (i == 0) {                                                         //first layer specific:
            sprintf(mv_outFiles_command, "mv %s/*.txt .", cons_meta[i].sn);   //generate command to move TECd_nav output files
            sprintf(rm_out_dirs_command, "rm -r %s", cons_meta[i].sn);        //generate command to remove TECd_nav output dir
                  
        } else {
            sprintf(mv_outFiles_command, "mv *_%s/*.txt .", cons_meta[i].sn); //non-first layer specific
            sprintf(rm_out_dirs_command, "rm -r *_%s", cons_meta[i].sn);      //generate command to remove TECd_nav output dir
        }
        
        system(mv_outFiles_command); //run move output files command
        system(rm_out_dirs_command); //run remove output directories command
        
        chdir(".."); //return to the parent directory
    }
    
    return;
}

/* aggregate_output: aggregate all fracBound (or user-specified) columns for a given layer into an output file */
void aggregate_output(int vals_cnt, values_input * vals, int layr_cnt, constraint_metadata * cons_meta, char * out_prefix, char * col_id)
{
    extern output_file_names out_fns[MAX_LAYERS]; //storage for output file names of each layer
    
    FILE ** TECDnav_out = NULL;     //pointers to TECdisplay_navigator output files
    FILE ** ofp = NULL;             //pointers to aggregated output files
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    char layr_dir[MAX_NAME+1] = {0};   //array for generating layer directory name
    char out_suffix[MAX_LINE+1] = {0}; //array for storing output file name, used to build output names across loop iterations
    char agg_out_nm[MAX_LINE+1] = {0}; //array for generating aggregated output file names for each values file
    char tmp_str[MAX_LINE+1] = {0};    //array for temporary string storage
    
    int EOF_tot = 0;                 //total EOFs reached
    int * EOF_rchd = NULL;           //array for tracking which EOFs have been reached
    
    int cols2incld[MAX_VALS] = {0};  //array for storing the included column numbers
    
    for (i = 0; i < layr_cnt; i++) { //for each layer
        
        EOF_tot = 0; //zero EOF tot for current layer
         
        sprintf(layr_dir, "layer%d", i+1); //generate current layer directory name
        chdir(layr_dir);                   //change to current layer directory
        
        //allocate file pointers for TECdisplay_navigator output files
        if ((TECDnav_out = calloc(out_fns[i].f_cnt, sizeof(*(TECDnav_out)))) == NULL) {
            printf("aggreggate_output: error - memory allocation for file pointers failed. aborting...\n");
            abort();
        }
        
        //open TECdisplay_navigator output files
        for (j = 0; j < out_fns[i].f_cnt; j++) {
            get_file(&TECDnav_out[j], out_fns[i].fn[j]);
        }
        
        //allocate memory for EOF_rchd array
        if ((EOF_rchd = calloc(out_fns[i].f_cnt, sizeof(*(EOF_rchd)))) == NULL) {
            printf("aggreggate_output: error - memory allocation for EOF detection failed. aborting...\n");
            abort();
        }
        
        //allocate memory for aggregated data output files
        if ((ofp = calloc(vals_cnt, sizeof(*(ofp)))) == NULL) {
            printf("aggreggate_output: error - memory allocation for output file pointers failed. aborting...\n");
            abort();
        }
        
        //generate aggregated data output file names
        if (!i) { //if processing layer zero, copy the constraint sample name to the output suffix
            strcpy(out_suffix, cons_meta[i].sn);
            
        } else {  //if processing a non-zero layer, append the constraint sample name to the output suffix
            strcpy(tmp_str, out_suffix);
            if ((snprintf(out_suffix, MAX_LINE, "%s_%s", tmp_str, cons_meta[i].sn)) >= MAX_LINE) {
                printf("aggregate_output: error - output suffix exceeded buffer. aborting...\n");
                abort();
            }
        }
        
        //generate an aggregated output file for each input values file
        for (j = 0; j < vals_cnt; j++) {
            
            //generate file name by appending out_suffix to the values file name
            if ((snprintf(agg_out_nm, MAX_LINE, "%s_%s.txt", vals[j].nm, out_suffix)) >= MAX_LINE) {
                printf("aggregate_output: error - aggregate output file name exceeded buffer. aborting...\n");
                abort();
            }
            
            //open the aggregated output file
            if ((ofp[j] = fopen(agg_out_nm, "w")) == NULL) {
                printf("aggregate_output: error - could not open aggregate output file. Aborting program...\n");
                abort();
            }
        }
        
        //read TECdisplay_navigator output file headers and print the headers of included colums to each output file
        read_output_hdrs(vals_cnt, col_id, i, TECDnav_out, &cols2incld[0], ofp);
        
        //read TECdisplay_navigator output file data lines and print included columns to the relevant output file
        while (EOF_tot < out_fns[i].f_cnt) {
            read_data_line(vals_cnt, i, TECDnav_out, &cols2incld[0], &EOF_rchd[0], &EOF_tot, ofp);
        }
        
        //close the TECdisplay_navigator output files that were aggregated
        for (j = 0; j < out_fns[i].f_cnt; j++) {
            
            /* close TECDnav_out */
            if (fclose(TECDnav_out[j]) == EOF) {
                printf("aggregate_output: error - error occurred when closing file. Aborting program...\n");
                abort();
            }
        }
        
        /* close aggreggate output files */
        for (j = 0; j < vals_cnt; j++) {
            if (fclose(ofp[j]) == EOF) {
                printf("aggregate_output: error - error occurred when closing file. Aborting program...\n");
                abort();
            }
        }
        
        free(TECDnav_out); //free TECdisplay_navigator output file points
        free(EOF_rchd);    //free the EOF tracking array
        chdir("..");       //move up one directory
    }
}

/* read_output_hdrs: read the header line of TECdisplay_navigator output files
 and print the headers of included columns to the aggregated data output file*/
void read_output_hdrs(int vals_cnt, char * col_id, int crnt_layr, FILE ** TECDnav_out, int * cols2incld, FILE ** ofp)
{
    
    extern output_file_names out_fns[MAX_LAYERS]; //storage for output file names of each layer
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
    
    int v = 0; //values index
    
    int col = 0;         //current column number
    int incldd_cols = 0; //number of included cols
    int line_end = 0;    //flag that the end of the line was reached
    
    char line[MAX_LINE+1] = {0};    //array to store line
    char tmp_val[MAX_LINE+1] = {0}; //array to store value strings
    
    char * p_col_id = NULL; //pointer to column id in value string
    
    //process header line for each TECdisplay_navigator output file of the current layer
    for (i = 0; i < out_fns[crnt_layr].f_cnt; i++) {
        
        v = 0;   //initialize the values index to zero
        col = 0; //initialize the column number to zero
        
        if (get_line(line, TECDnav_out[i])) { //while there is still a line to get
            
            tmp_val[0] = '\0'; //intialize the first member of the value string array to zero
            
            for (line_end = 0, j = 0, k = 0; !line_end; j++) { //until the end of the line is reached
                
                if (line[j] != '\t' && line[j]) { //if the char is not a tab or a null char
                    tmp_val[k++] = line[j];       //copy the char to the value string array
                    
                } else { //if the char is a tab or a null char
                    
                    if (!line[j]) {   //if the char is a null char
                        line_end = 1; //set line_end to true
                    }
                    
                    tmp_val[k] ='\0'; //terminate the value string
                    k = 0;            //reset k to zero for next iteration
                    
                    if ((p_col_id = strstr(tmp_val, col_id)) != NULL) { //if the column header contains the target column id
                        
                        //TODO: consider being more flexible about requirements for non standard col ids?
                        if (p_col_id[-1] == '_' && !p_col_id[strlen(p_col_id)]) { //check that the col id is at the header end
                            
                            if (v >= vals_cnt) { //check that number of included columns will not exceed number of vals files
                                printf("read_output_hdrs: error - number of columns with match to target column id exceeds the number of input values files. aborting...\n");
                                abort();
                            }

                            if (!i) {                //if reading the first file of the current layer
                                cols2incld[v] = col; //store the number of the included column at the current vals file index
                                incldd_cols++;       //increment the number of included columns
                                fprintf(ofp[v++], "%s", tmp_val); //print the header to the current vals file and
                                                                  //increment the values file index
                                
                            } else if (col == cols2incld[v]) { //if reading non-first file, check cols2incld is same as 1st file
                                fprintf(ofp[v++], "\t%s", tmp_val); //print the header to the curent vals file and
                                                                    //increment the values file index
                                
                            } else { //included column numbers do not match across files
                                printf("read_output_hdrs: error - columns to include are not the same in all files. aborting...\n");
                                abort();
                            }
                        }
                    }
                    
                    col++; //increment column index
                }
            }
        } else {
            //TODO: add file name to message
            printf("read_output_hdrs: error - TECdisplay_navigator output file %s did not contain any headers. aborting...\n", out_fns[crnt_layr].fn[i]);
            abort();
        }
        
        if (v != vals_cnt) { //check that number of included columns matches values file count
            printf("read_output_hdrs: error - number of included columns (%d) does not match the number of values files (%d). aborting...\n", v, vals_cnt);
            abort();
        }
    }
    
    //print a newline to all aggregated output files
    for (v = 0; v < vals_cnt; v++) {
        fprintf(ofp[v], "\n");
    }
}

/* read_data_line: read data lines from TECdisplay_navigator output files and
 print the included columns to the relevant aggregated data output file */
void read_data_line(int vals_cnt, int crnt_layr, FILE ** TECDnav_out, int * cols2incld, int * EOF_rchd, int * EOF_tot, FILE ** ofp)
{
    extern output_file_names out_fns[MAX_LAYERS]; //storage for output file names of each layer
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
    
    int v = 0; //values index
    
    int col = 0;      //current column number
    int line_end = 0; //flag that the end of the line was reached
        
    char line[MAX_LINE+1] = {0};    //array to store line
    char tmp_val[MAX_LINE+1] = {0}; //array to store value strings
    
    char * p_col_id = NULL; //pointer to column id in value string
    
    //process data line for each TECdisplay_navigator output file of the current layer
    for (i = 0; i < out_fns[crnt_layr].f_cnt; i++) {
        
        v = 0;   //zero values index
        col = 0; //zero column number
        
        if (!EOF_rchd[i]) { //if the EOF was not reached for the current TECdisplay_navigator output file
            
            if (get_line(line, TECDnav_out[i])) { //get the next line
                
                tmp_val[0] = '\0'; //intialize the first member of the value string array to zero
                
                for (line_end = 0, j = 0, k = 0; !line_end; j++) {
                    
                    if (line[j] != '\t' && line[j]) { //if the char is not a tab or a null char
                        tmp_val[k++] = line[j];       //copy the char to the value string array
                        
                    } else { //if the char is a tab or a null char
                        
                        if (!line[j]) {   //if the char is a null char
                            line_end = 1; //set line_end to true
                        }
                        
                        tmp_val[k] ='\0'; //terminate the value string
                        k = 0;            //reset k to zero for next iteration
                        
                        if (col == cols2incld[v]) { //if the current column is to be included. note it is confirmed that the
                                                    //number of cols 2 include does not exceed the number of values files above
                                                    //in the read_output_hdrs function
                            
                            if (!i) {                               //if reading first file line
                                fprintf(ofp[v++], "%s", tmp_val);   //print value w/o leading tab
                                
                            } else {                                //if reading non-first file line
                                fprintf(ofp[v++], "\t%s", tmp_val); //print value with leading tab
                            }
                        }
                        
                        col++; //increment column index
                    }
                }
            } else { //EOF of the current file was reached on this loop iteration
                
                EOF_rchd[i] = 1; //set EOF_rchd flag for current file to true
                (*EOF_tot)++;    //increment count of total EOF reached
                
                for (v = 0; v < vals_cnt; v++) { //for each aggregated data output file
                    if (i) {                     //if reading a non-first file line
                        fprintf(ofp[v], "\t");   //print a tab
                    }
                }
            }
            
        } else { //EOF of the current file was reached on a previous loop iteration
            
            for (v = 0; v < vals_cnt; v++) { //for each aggregated data output file
                if (i) {                     //if reading a non-first file line
                    fprintf(ofp[v], "\t");   //print a tab
                }
            }
        }
        
        if (v != vals_cnt) { //check that number of included columns matches values file count
            printf("read_output_hdrs: error - number of included columns (%d) does not match the number of values files (%d). aborting...\n", v, vals_cnt);
            abort();
        }
    }
    
    //print a newline to all aggregated output files
    for (v = 0; v < vals_cnt; v++) {
        fprintf(ofp[v], "\n");
    }
}

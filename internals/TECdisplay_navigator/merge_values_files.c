//
//  merge_values_files.c
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"

#include "../TECdisplay_mapper/TECdisplay_output_column_headers.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"
#include "./parse_TECdisplay_out_line.h"

#include "merge_values_files.h"

/* merge_values_files: merge multiple values files into a single file */
//TODO: change vals_cnt variable to something more representative of what it is
int merge_values_files(values_input * vals_ipt, int vals_cnt, char * merged_out_nm, int nonstandard)
{
    extern const char TECdsply_clmn_hdrs[4][32]; //column headers from TECdisplay_output_column_headers.c
    
    int v = 0;  //index of values input file
    int i = 0;  //general purpose index variable
    
    int got_line[MAX_VALS] = {0};                    //array to flag the success of get_line for each values file
    int proceed = 1;                                 //flag to continue processing
    
    int line_cnt = 0;                                //tracks number of lines in values files
    char line[MAX_VALS][MAX_LINE] = {{0}};           //arrays to store lines from each values file
    char *p_id[MAX_VALS] = {NULL};                   //pointers to variant id
    char *p_vals[MAX_VALS] = {NULL};                 //pointers to variant values
    
    int bnd[MAX_VALS] = {0};                         //bound read count values
    int unb[MAX_VALS] = {0};                         //unbound read count values
    double frc[MAX_VALS] = {0};                      //fraction bound read count values
    
    int field_count = 0;                             //variable for counting number of fields in values file line
    int found_term = 0;                              //flag that terminating null was found
    
    /* generate merged output file */
    FILE * vals_merged = NULL; //pointer for merged output file
    if ((vals_merged = fopen(merged_out_nm, "w")) == NULL) {
        printf("merge_values_files: error - could not open merged output file. Aborting program...\n");
        abort();
    }
    fprintf(vals_merged, "variant_%s", TECdsply_clmn_hdrs[TDSPLY_VID_HDR]); //print id column header
    
    for (line_cnt = 0; proceed; line_cnt++) {   //iterate through each line of the input files
        for (v = 0; v < vals_cnt; v++) {        //process each input file
            
            /* get line from current values file, test that success status is the same
             for all input files. if success status is not the same for all files, the
             input files are not the same length and cannot be processed together. */
            got_line[v] = (get_line(line[v], vals_ipt[v].fp)) ? 1 : 0;  //if successful, got_line=1, else got_line=0
            if (got_line[v] != got_line[0]) { //test if got_line success status equals that of the first input
                printf("merge_values_files: error - input values files are not the same length. aborting...\n");
                abort();
            }
            
            if (got_line[v]) {

                //in each iteration, the number of tabs that precede a non-null character
                //are counted to determine the number of fields in the data line.
                if (line_cnt == 0) { //reading column header line (first line of each file)
                    
                    parse_TECdisplay_out_line(&line[v][0], &p_id[v], &p_vals[v], NULL, NULL, NULL, TDSPLY_HDR_LINE, nonstandard);
                    
                    if (!nonstandard) { //input data uses standard TECdisplay format, print standard headers
                        
                        /* print column headers for the current values file to the merged output file */
                        fprintf(vals_merged, "\t%s_%s\t%s_%s\t%s_%s",
                                vals_ipt[v].nm, TECdsply_clmn_hdrs[TDSPLY_BND_HDR],
                                vals_ipt[v].nm, TECdsply_clmn_hdrs[TDSPLY_UNB_HDR],
                                vals_ipt[v].nm, TECdsply_clmn_hdrs[TDSPLY_FRC_HDR]);
                        
                    } else { //input data files use non-standard format
                                                
                        fprintf(vals_merged, "\t%s", p_vals[v]); //print non-standard headers
                    }
                    
                } else { //reading data line
                    
                    parse_TECdisplay_out_line(&line[v][0], &p_id[v], &p_vals[v], &bnd[v], &unb[v], &frc[v], TDSPLY_DATA_LINE, nonstandard);
                    
                    /* data lines are validated by comparing the variant id from each input file
                     to the variant id of the first input file. all variant ids must be identical
                     in order for the values files to be merged. */
                    if (strcmp(p_id[v], p_id[0])) {
                        printf("%s\n%s\n", p_id[v], p_id[0]);
                        printf("merge_values_files: error - variant ids are not aligned. aborting...\n");
                        abort();
                    }
                    
                    //all data checks passed, print output
                    if (v == 0) {                            //if reading line from first values file
                        fprintf(vals_merged, "%s", p_id[0]); //print variant id to merged output file
                    }
                    fprintf(vals_merged, "\t%s", p_vals[v]); //print values string to merged output file
                }
            }
        }
        
        if (!got_line[0]) { //reached last line. to reach this part of the code, the end
            proceed = 0;    //of every values file must have been reached simultaneously
        } else {
            fprintf(vals_merged, "\n"); //printed data line, so print newline to start next line
        }
        
    }
    
    /* close merged values file */
    if (fclose(vals_merged) == EOF) {
        printf("merge_values_files: error - error occurred when closing merged input file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

//
//  merge_TECdisplay_replicates.c
//  
//
//  Created by Eric Strobel on 6/25/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <sys/stat.h>

#include "../utils/io_management.h"
#include "../utils/gen_utils.h"

#include "../TECdisplay_mapper/TECdisplay_output_column_headers.h"
#include "../TECdisplay_mapper/map_reads/map_expected/mk_output_files.h"

#include "../TECdisplay_navigator/TECdisplay_navigator_defs.h"
#include "../TECdisplay_navigator/TECdisplay_navigator_structs.h"


void merge_replicates(values_input * vals_ipt, int vals_cnt, char * out_nm);

int main(int argc, char *argv[])
{
    extern int debug; //flag to run debug mode
    
    values_input vals[MAX_VALS] = {0}; //storage for sample name, file name, and file pointer of values files
    
    int vals_cnt = 0;                  //number of values files provided
    int out_nm_provided = 0;           //number of output file names provided
    
    char out_nm[MAX_NAME] = {0};       //array to store output file name
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"values",        required_argument,  0,  'v'},  //values file input
            {"debug",         required_argument,  0,  'd'},  //turn on debug mode
            {"out-name",      required_argument,  0,  'o'},  //set output file name
            
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "v:c:d:o:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
                
           
            case 'v': //get values file
                if (vals_cnt < MAX_VALS) {                              //check that MAX_VALS will not be exceeded
                    strcpy(vals[vals_cnt].fn, argv[optind-1]);          //store file name
                    get_sample_name(argv[optind-1], vals[vals_cnt].nm); //get sample name from file name
                    get_file(&(vals[vals_cnt].fp), argv[optind-1]);     //set file pointer to values file
                    vals_cnt++;                                         //count values files provided
                } else {
                    printf("merge_TECdisplay_replicates: error - too many values files were provided. the maximum is %d. aborting...\n", MAX_VALS);
                    abort();
                }
                break;
                
            case 'o': //set output_directory name
                if (!out_nm_provided) {                        //check that output directory was not previously provided
                    if (strlen(argv[optind-1]) < (MAX_NAME)) { //check output directory name length
                        strcpy(out_nm, argv[optind-1]);        //store output file name
                    } else {
                        printf("merge_TECdisplay_replicates: error - output directory name is longer than the maximum length (%d). setting output directory to default name 'out'.\n", MAX_NAME);
                    }
                    out_nm_provided++; //increment output name counter
                    
                } else {
                    printf("merge_TECdisplay_replicates: error - more than one output directory name was provided. aborting\n");
                    abort();
                }
                break;
            
            case 'd': //turn on debug mode
                debug = 1;
                break;
                
            default:
                printf("error: unrecognized option. Aborting program...\n");
                abort();
        }
    }
    
    if (optind < argc) {
        printf("\nnon-option ARGV-elements:\n");
        while (optind < argc)
            printf("\n%s \n", argv[optind++]);
        putchar('\n');
        printf("Aborting program.\n\n");
        abort();
    }
    /*********** end of option parsing ***********/
    
    if (!out_nm_provided) {
        printf("merge_TECdisplay_replicates: error - no output name was provided. aborting...\n");
        abort();
    }
    
    //TODO: check for existing files with same output name
    //TODO: check whether output name was supplied with a suffix?
    
    merge_replicates(&vals[0], vals_cnt, out_nm);  //merge replicate values
}

/* merge_replicates: generate output file that combines the read counts for replicate bound and unbound channels */
void merge_replicates(values_input * vals_ipt, int vals_cnt, char * out_nm)
{
    extern const char TECdsply_clmn_hdrs[4][32]; //column headers from TECdisplay_output_column_headers.c
    
    int v = 0;  //index of values input file
    int i = 0;  //general purpose index variable
    int j = 0;  //general purpose index variable
    int k = 0;  //general purpose index variable
    
    int got_line[MAX_VALS] = {0};                    //array to flag the success of get_line for each values file
    int proceed = 1;                                 //flag to continue processing
    
    int line_cnt = 0;                                //tracks number of lines in values files
    char line[MAX_VALS][MAX_LINE] = {{0}};           //arrays to store lines from each values file
    char col_nm[XPCTD_FIELDS][MAX_COL_NM+1] = {{0}}; //array to store input column names for error-checking
    char *p_id[MAX_VALS] = {NULL};                   //pointers to variant id
    char *p_vals[MAX_VALS] = {NULL};                 //pointers to variant values
    int field_count = 0;                             //variable for counting number of fields in values file line
    int found_term = 0;                              //flag that terminating null was found
    
    char *crrnt_val = NULL;                          //pointer to current value field string
    int crrnt_bnd_cnt = 0;                           //bound read total for current variant
    int crrnt_unb_cnt = 0;                           //unbound read total for current variant
    
    /* generate merged output file */
    FILE * vals_merged = NULL;                       //pointer for merged output file
    char output_file_name[MAX_LINE] = {0};           //array to store output file name
    
    sprintf(output_file_name, "%s.txt", out_nm);     //generate output file name by appending '.txt' suffix

    if ((vals_merged = fopen(output_file_name, "w")) == NULL) { //open output file
        printf("merge_replicates: error - could not open merged output file. Aborting program...\n");
        abort();
    }
    
    print_output_header(vals_merged, out_nm);   //print output file headers

    for (line_cnt = 0; proceed; line_cnt++) {   //iterate through each line of the input files
        crrnt_bnd_cnt = crrnt_unb_cnt = 0;      //zero read count totals
        for (v = 0; v < vals_cnt; v++) {        //process each input file
            
            /* get line from current values file, test that success status is the same
             for all input files. if success status is not the same for all files, the
             input files are not the same length and cannot be processed together. */
            got_line[v] = (get_line(line[v], vals_ipt[v].fp)) ? 1 : 0;  //if successful, got_line=1, else got_line=0
            if (got_line[v] != got_line[0]) { //test if got_line success status equals that of the first input
                printf("merge_replicates: error - input values files are not the same length. aborting...\n");
                abort();
            }
            
            if (got_line[v]) {

                //in each iteration, the number of tabs that precede a non-null character
                //are counted to determine the number of fields in the data line.

                if (line[v][0] == '\t') { //if id field is missing entry, abort
                    printf("merge_replicates: error - values file lines cannot begin with a tab character. aborting...\n");
                    abort();
                    
                } else if (line[v][0]) {  //line contains at least one field
                    field_count = 1;      //initialize field count to 1
                    
                } else { //it should not be possible to have an empty line here
                    printf("merge_replicates: error - empty line, this should not be possible. aborting...\n");
                    abort();
                }
                
                //split line into id and values strings
                for (i = 0; line[v][i] != '\t' && line[v][i] && i < MAX_LINE; i++) {;} //iterate to first tab
                
                if (line[v][i] == '\t') {       //if the tab separating the id field from the data fields was found
                    line[v][i] = '\0';          //set first tab to null char to split input line
                    p_id[v] = &line[v][0];      //set pointer to variant id string
                    p_vals[v] = &line[v][i+1];  //set pointer to values string
                    
                } else { //unrecognized format error
                    printf("%s\n", line[v]);
                    printf("merge_replicates: error - unrecognized line format. aborting...\n");
                    abort();
                }
                
                if (line_cnt == 0) { //reading column header line (first line of each file)
                    if (v == 0) { //reading column header line of first values file
                        
                        /* when reading the column headers of the first values file, parse
                         the headers and validate against the expected headers (variant_id,
                         bound, unbound, fracBound). all subsequent header lines can then
                         be compared directly against those from the first values file*/
                                                
                        i = 0;                        //initialize i (column name index) to zero
                        strcpy(col_nm[i++], p_id[0]); //copy the variant id column name to col_nm[0], increment i
                        
                        //parse values headers and store in col_nm array
                        for (j = 0; p_vals[0][j] && i < XPCTD_FIELDS; i++, field_count++) {
                            
                            //copy column header into col_nm array
                            for (k = 0; p_vals[0][j] != '\t' && p_vals[0][j] && k < (MAX_COL_NM); j++, k++) {
                                col_nm[i][k] = p_vals[0][j];
                            }
                            col_nm[i][k] = '\0';
                                                        
                            //test that loop exited on a tab or null character.
                            //if test fails, column name is too long.
                            if (p_vals[0][j] == '\t') {        //ended on tab, more headers
                                j++;                           //increment j
                                
                            } else if (p_vals[0][j] != '\0') { //if loop did not exit on tab or null, name is too long
                                printf("merge_replicates: error - unexpected long column name in values file. aborting...\n");
                                abort();
                            }
                        }
                        
                        //test that loop exited on a null character and that
                        //the expected number of value fields were identified
                        if (p_vals[0][j] || field_count != XPCTD_FIELDS) {
                            printf("merge_replicates: error - unexpected headers for values files. aborting...\n");
                            abort();
                        }
                        
                        //test that column headers match expected strings
                        //strstr is used for variant id header because header can be variable but always
                        //contains the substring "id". other headers are not variable, so test is for
                        //an exact match
                        //TODO: make id test look only at last chars of col name?
                        if (strstr(col_nm[TDSPLY_VID_HDR], TECdsply_clmn_hdrs[TDSPLY_VID_HDR]) == NULL ||
                            strcmp(col_nm[TDSPLY_BND_HDR], TECdsply_clmn_hdrs[TDSPLY_BND_HDR])         ||
                            strcmp(col_nm[TDSPLY_UNB_HDR], TECdsply_clmn_hdrs[TDSPLY_UNB_HDR])         ||
                            strcmp(col_nm[TDSPLY_FRC_HDR], TECdsply_clmn_hdrs[TDSPLY_FRC_HDR])) {
                            printf("merge_replicates: error - unexpected headers for values files. aborting...\n");
                            abort();
                        }

                    /* compare all other values file headers to the headers from
                     the first values file, which has already been validated */
                    //TODO: make id test look only at last chars of col name?
                    } else if (strstr(p_id[v], TECdsply_clmn_hdrs[TDSPLY_VID_HDR]) == NULL ||
                               strcmp(p_vals[v], p_vals[0])) {
                        printf("merge_replicates: error - supplied values files contain discordant headers. aborting...\n");
                        abort();
                    }
                    
                } else { //reading data line
                    
                    /* data lines are validated by comparing the variant id from each input file
                     to the variant id of the first input file. all variant ids must be identical
                     in order for the values files to be merged. */
                    if (strcmp(p_id[v], p_id[0])) {
                        printf("%s\n%s\n", p_id[v], p_id[0]);
                        printf("merge_replicates: error - variant ids are not aligned. aborting...\n");
                        abort();
                        
                    } else {
                        
                        /* read p_vals string. confirm expected number of fields (count started
                         above) and to check that all characters are of expected types. total
                         read counts for each channel of the libray (bound and unbound) are
                         determined for the current variant. */
                        
                        for (i = 0, crrnt_val = &p_vals[v][0], found_term = 0; !found_term && i < MAX_LINE; i++) {
                            if (p_vals[v][i] == '\t' || !p_vals[v][i]) { //if a tab or term null was reached
                                
                                if (!p_vals[v][i]) {                            //found terminating null
                                    found_term = 1;                             //set found terminating null flag
                                }
                                
                                if (isdigit(p_vals[v][i-1])) {                  //if the preceding character was a digit
                                    p_vals[v][i] = '\0';                        //set the delimiter to a term null
                                    
                                    if (field_count == TDSPLY_BND_HDR) {        //if reading bound reads field
                                        crrnt_bnd_cnt += atoi(crrnt_val);       //add reads to bound reads total
                                    
                                    } else if (field_count == TDSPLY_UNB_HDR) { //if reading unbound reads field
                                        crrnt_unb_cnt += atoi(crrnt_val);       //add reads to unbound reads total
                                    }
                                    
                                    field_count++;                              //increment field count
                                }
                                
                                crrnt_val = &p_vals[v][i+1];                    //set pointer to start of next value field
                                
                            } else if (!isdigit(p_vals[v][i]) && //char is not a digit
                                       p_vals[v][i] != '.'    && //or a '.'
                                       p_vals[v][i] != '-') {    //or a '-'
                                
                                printf("merge_replicates: error - unexpected character %c (ASCII: %d) in values string. aborting...\n", p_vals[v][i], p_vals[v][i]);
                                abort();
                            }
                        }
                        
                        if (!found_term) { //if the terminating null was not found, throw error and abort
                            printf("merge_replicates: unexpected long data values line. aborting...\n");
                            abort();
                        }
                        
                        if (field_count != XPCTD_FIELDS) { //check if correct number of fields was read
                            printf("%d\n", field_count);
                            printf("merge_replicates: error - unexpected number of fields in data line. aborting...\n");
                            abort();
                        }
                    }
                }
            }
        }
        
        if (!got_line[0]) { //reached last line. to reach this part of the code, the end
            proceed = 0;    //of every values file must have been reached simultaneously
        }
        
        if (line_cnt && proceed) { //read data line, print values to output file
            for (i = 0; i < TDSPLY_HDR_CNT; i++) {
                switch (i) {
                    case TDSPLY_VID_HDR: //printing variant id column
                        fprintf(vals_merged, "%s", p_id[0]);
                        break;
                        
                    case TDSPLY_BND_HDR: //printing bount read count
                        fprintf(vals_merged, "\t%d", crrnt_bnd_cnt);
                        break;
                        
                    case TDSPLY_UNB_HDR: //printing unbound read count
                        fprintf(vals_merged, "\t%d", crrnt_unb_cnt);
                        break;
                        
                    case TDSPLY_FRC_HDR: //printing fraction bound
                        fprintf(vals_merged, "\t%.4f",
                                (double)(crrnt_bnd_cnt)/(double)(crrnt_bnd_cnt + crrnt_unb_cnt));
                        break;
                        
                    default: //this should not be reachable
                        printf("print_output: error - header index exceeded the maximum. aborting...\n");
                        abort();
                        break;
                }
            }
            fprintf(vals_merged, "\n"); //print newline at the end of data line
        }
    }
    
    /* close merged values file */
    if (fclose(vals_merged) == EOF) {
        printf("merge_replicates: error - error occurred when closing merged input file. Aborting program...\n");
        abort();
    }

    return;
}

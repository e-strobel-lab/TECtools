//
//  calc_FracBound_difference.c
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

#include "../TECdisplay_navigator/parse_TECdisplay_out_line.h"

#define tMIN 0
#define tSUB 1
#define MAX_TERMS 2

/* calc_FracBound_difference: generate output file that combines the read counts for replicate bound and unbound channels */
void calc_FracBound_difference(values_input * ipt, char * out_nm);

int main(int argc, char *argv[])
{
    extern int debug;              //flag to run debug mode
    
    values_input ipt[MAX_TERMS] = {0}; //storage for sample name, file name, and file pointer of input files
    
    int minuend_provided = 0;      //flag that minuend file was provided
    int subtrahend_provided = 0;   //flag that subtrahend file was provided
    int out_nm_provided = 0;       //number of output file names provided
    
    char out_nm[MAX_NAME] = {0};   //array to store output file name
    
    int i = 0;                     //general purpose index
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"minuend",       required_argument,  0,  'm'},  //minuend file input
            {"subtrahend",    required_argument,  0,  's'},  //subtrahend file input
            {"debug",         required_argument,  0,  'd'},  //turn on debug mode
            {"out-name",      required_argument,  0,  'o'},  //set output file name
            
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "m:s:d:o:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
                           
            case 'm': //get minuend file
                if (!minuend_provided) { //check that minuend file was not yet provided
                    strcpy(ipt[tMIN].fn, argv[optind-1]);          //store file name
                    get_sample_name(argv[optind-1], ipt[tMIN].nm); //get sample name from file name
                    get_file(&(ipt[tMIN].fp), argv[optind-1]);     //set file pointer
                    minuend_provided++;                      //count minuend files provided
                } else {
                    printf("calc_FracBound_difference: error - more than one minuend file was provided. aborting...\n");
                    abort();
                }
                break;
                
            case 's': //get subtrahend file
                if (!subtrahend_provided) { //check that subtrahend file was not yet provided
                    strcpy(ipt[tSUB].fn, argv[optind-1]);          //store file name
                    get_sample_name(argv[optind-1], ipt[tSUB].nm); //get sample name from file name
                    get_file(&(ipt[tSUB].fp), argv[optind-1]);     //set file pointer
                    subtrahend_provided++;                         //count subtrahend files provided
                } else {
                    printf("calc_FracBound_difference: error - more than one subtrahend file was provided. aborting...\n");
                    abort();
                }
                break;
                
            case 'o': //set output file name
                if (!out_nm_provided) {  //check that output file was not previously provided
                    
                    if (strlen(argv[optind-1]) < MAX_NAME) { //check output file name length
                        strcpy(out_nm, argv[optind-1]);      //store output file name
                        
                        //check that output file name contains only alphanumeric and '_' characters
                        for (i = 0; out_nm[i]; i++) {
                            if (!isalnum(out_nm[i]) && out_nm[i] != '_') {
                                printf("calc_FracBound_difference: error output file name contains the invalid character %c (%d). aborting...\n", out_nm[i], out_nm[i]);
                                abort();
                            }
                        }
                    } else {
                        printf("calc_FracBound_difference: error - output file name is longer than the maximum length (%d). aborting...\n", MAX_NAME);
                        abort();
                    }
                    out_nm_provided++; //increment output name counter
                    
                } else {
                    printf("calc_FracBound_difference: error - more than one output file name was provided. aborting\n");
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
        printf("calc_FracBound_difference: error - no output name was provided. aborting...\n");
        abort();
    }
    
    /* generate input record */
    FILE * input_record = NULL;             //pointer for input record file
    char output_file_name[MAX_LINE] = {0};  //array to store input record file name
    
    sprintf(output_file_name, "%s_input_record.txt", out_nm); //generate input record file name by appending '_input_record.txt' suffix

    if ((input_record = fopen(output_file_name, "w")) == NULL) { //open input record file
        printf("calc_FracBound_difference: error - could not open input record file. Aborting program...\n");
        abort();
    }
    
    //print the input and output file names
    fprintf(input_record, "minuend:    %s\n", ipt[tMIN].fn);
    fprintf(input_record, "subtrahend: %s\n", ipt[tSUB].fn);
    fprintf(input_record, "output:     %s.txt\n", out_nm);
    
    /* close merged values file */
    if (fclose(input_record) == EOF) {
        printf("calc_FracBound_difference: error - error occurred when closing merged input file. Aborting program...\n");
        abort();
    }
        
    calc_FracBound_difference(&ipt[0], out_nm);  //calculate difference
}

/* calc_FracBound_difference: generate output file that combines the read counts for replicate bound and unbound channels */
void calc_FracBound_difference(values_input * ipt, char * out_nm)
{
    extern const char TECdsply_clmn_hdrs[4][32]; //column headers from TECdisplay_output_column_headers.c
    
    int v = 0;  //index of values input file
    int i = 0;  //general purpose index variable
    
    int got_line[MAX_TERMS] = {0};                     //array to flag the success of get_line for each values file
    int proceed = 1;                                   //flag to continue processing
    
    int line_cnt = 0;                                  //tracks number of lines in values files
    char line[MAX_TERMS][MAX_LINE] = {{0}};            //arrays to store lines from each values file
    char *p_id[MAX_TERMS] = {NULL};                    //pointers to variant id
    char *p_vals[MAX_TERMS] = {NULL};                  //pointers to variant values
    int field_count = 0;                               //variable for counting number of fields in values file line
    int found_term = 0;                                //flag that terminating null was found
    
    int bnd = 0;     //bound read count values
    int unb = 0;     //unbound read count values
    double frc = 0;  //fraction bound read count values
    
    int crrnt_bnd[MAX_TERMS] = {0};                  //bound read count for current variant
    int crrnt_unb[MAX_TERMS] = {0};                  //unbound read count for current variant
    double crrnt_frc[MAX_TERMS] = {0};            //frac bound for current variant
    
    /* generate merged output file */
    FILE * vals_merged = NULL;                       //pointer for merged output file
    char output_file_name[MAX_LINE] = {0};           //array to store output file name
    
    sprintf(output_file_name, "%s.txt", out_nm);     //generate output file name by appending '.txt' suffix

    if ((vals_merged = fopen(output_file_name, "w")) == NULL) { //open output file
        printf("calc_FracBound_difference: error - could not open merged output file. Aborting program...\n");
        abort();
    }
        
    fprintf(vals_merged, "%s", TECdsply_clmn_hdrs[TDSPLY_VID_HDR]);
    for (v = 0; v < MAX_TERMS; v++) {
        for (i = 1; i < TDSPLY_HDR_CNT; i++) {
            fprintf(vals_merged, "\t%s_%s", ipt[v].nm, TECdsply_clmn_hdrs[i]);  //print "<tab><oustatements"
        }
    }
    fprintf(vals_merged, "\tdiff");
    //fprintf(vals_merged, "\ttot1\ttot2\ttot_all");
    fprintf(vals_merged, "\n");

    for (line_cnt = 0; proceed; line_cnt++) {   //iterate through each line of the input files
        crrnt_bnd[tMIN] = crrnt_bnd[tSUB] = 0;
        crrnt_unb[tMIN] = crrnt_unb[tSUB] = 0;
        crrnt_frc[tMIN] = crrnt_frc[tSUB] = 0;
        
        for (v = 0; v < MAX_TERMS; v++) {        //process each input file
            
            /* get line from current values file, test that success status is the same
             for all input files. if success status is not the same for all files, the
             input files are not the same length and cannot be processed together. */
            got_line[v] = (get_line(line[v], ipt[v].fp)) ? 1 : 0;  //if successful, got_line=1, else got_line=0
            if (got_line[v] != got_line[0]) { //test if got_line success status equals that of the first input
                printf("calc_FracBound_difference: error - input values files are not the same length. aborting...\n");
                abort();
            }
            
            if (got_line[v]) {
                
                if (line_cnt == 0) { //reading column header line (first line of each file)
                                            
                    parse_TECdisplay_out_line(&line[v][0], &p_id[v], &p_vals[v], NULL, NULL, NULL, TDSPLY_HDR_LINE, 0);
                    
                } else { //reading data line
                    
                    crrnt_bnd[v] = 0; //initialize bound read count to zero
                    crrnt_unb[v] = 0; //initialize unbound read count to zero
                    crrnt_frc[v] = 0; //initialize fraction bound to zero
                    
                    parse_TECdisplay_out_line(&line[v][0], &p_id[v], &p_vals[v], &crrnt_bnd[v], &crrnt_unb[v], &crrnt_frc[v], TDSPLY_DATA_LINE, 0);
                    
                    /* data lines are validated by comparing the variant id from each input file
                     to the variant id of the first input file. all variant ids must be identical
                     in order for the values files to be merged. */
                    if (strcmp(p_id[v], p_id[0])) {
                        printf("%s\n%s\n", p_id[v], p_id[0]);
                        printf("calc_FracBound_difference: error - variant ids are not aligned. aborting...\n");
                        abort();
                    }
                }
            }
        }
        
        if (!got_line[0]) { //reached last line. to reach this part of the code, the end
            proceed = 0;    //of every values file must have been reached simultaneously
        }
        
        if (line_cnt && proceed) { //read data line, print values to output file
            
            fprintf(vals_merged, "%s", p_id[0]);
            
            for (v = 0; v < MAX_TERMS; v++) {
                fprintf(vals_merged, "\t%s", p_vals[v]);
            }
            
            fprintf(vals_merged, "\t%.4f", crrnt_frc[tMIN] - crrnt_frc[tSUB]);
            //fprintf(vals_merged, "\t%d\t%d\t%d", crrnt_bnd[tMIN] + crrnt_unb[tMIN], crrnt_bnd[tSUB] + crrnt_unb[tSUB], crrnt_bnd[tMIN] + crrnt_unb[tMIN] + crrnt_bnd[tSUB] + crrnt_unb[tSUB]);
            fprintf(vals_merged, "\n");
            
        }
    }
    
    /* close merged values file */
    if (fclose(vals_merged) == EOF) {
        printf("calc_FracBound_difference: error - error occurred when closing merged input file. Aborting program...\n");
        abort();
    }

    return;
}



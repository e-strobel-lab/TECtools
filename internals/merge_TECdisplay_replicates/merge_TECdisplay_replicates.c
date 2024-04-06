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

#include "../TECdisplay_navigator/parse_TECdisplay_out_line.h"

/* merge_replicates: generate output file that combines the read counts for replicate bound and unbound channels */
void merge_replicates(values_input * vals_ipt, int vals_cnt, char * out_nm);

int main(int argc, char *argv[])
{
    extern int debug;                  //flag to run debug mode
    
    values_input vals[MAX_VALS] = {0}; //storage for sample name, file name, and file pointer of values files
    
    int vals_cnt = 0;                  //number of values files provided
    int out_nm_provided = 0;           //number of output file names provided
    
    char out_nm[MAX_NAME] = {0};       //array to store output file name
    
    int i = 0;                         //general purpose index
    
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
        
        c = getopt_long(argc, argv, "v:d:o:", long_options, &option_index);
        
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
                
            case 'o': //set output file name
                if (!out_nm_provided) {  //check that output file was not previously provided
                    
                    if (strlen(argv[optind-1]) < MAX_NAME) { //check output file name length
                        strcpy(out_nm, argv[optind-1]);      //store output file name
                        
                        //check that output file name contains only alphanumeric and '_' characters
                        for (i = 0; out_nm[i]; i++) {
                            if (!isalnum(out_nm[i]) && out_nm[i] != '_') {
                                printf("merge_TECdisplay_replicates: error output file name contains the invalid character %c (%d). aborting...\n", out_nm[i], out_nm[i]);
                                abort();
                            }
                        }
                    } else {
                        printf("merge_TECdisplay_replicates: error - output file name is longer than the maximum length (%d). aborting...\n", MAX_NAME);
                        abort();
                    }
                    out_nm_provided++; //increment output name counter
                    
                } else {
                    printf("merge_TECdisplay_replicates: error - more than one output file name was provided. aborting\n");
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
    
    /* generate merge record */
    FILE * merge_record = NULL;             //pointer for merge record file
    char output_file_name[MAX_LINE] = {0};  //array to store merge record file name
    
    sprintf(output_file_name, "%s_replicate_merge_record.txt", out_nm); //generate merge record file name by appending '_replicate_merge_record.txt' suffix

    if ((merge_record = fopen(output_file_name, "w")) == NULL) { //open merge record file
        printf("merge_TECdisplay_replicates: error - could not open merge record file. Aborting program...\n");
        abort();
    }
    
    //print the input and output file names
    for (i = 0; i < vals_cnt; i++) {
        fprintf(merge_record, "input%d: %s\n", i+1, vals[i].fn);
    }
    fprintf(merge_record, "output: %s.txt\n", out_nm);
    
    /* close merged values file */
    if (fclose(merge_record) == EOF) {
        printf("merge_TECdisplay_replicates: error - error occurred when closing merged input file. Aborting program...\n");
        abort();
    }
        
    merge_replicates(&vals[0], vals_cnt, out_nm);  //merge replicate values
}

/* merge_replicates: generate output file that combines the read counts for replicate bound and unbound channels */
void merge_replicates(values_input * vals_ipt, int vals_cnt, char * out_nm)
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
    int field_count = 0;                             //variable for counting number of fields in values file line
    int found_term = 0;                              //flag that terminating null was found
        
    int bnd = 0;     //bound read count values
    int unb = 0;     //unbound read count values
    double frc = 0;  //fraction bound read count values
    
    int crrnt_bnd_cnt = 0;  //bound read total for current variant
    int crrnt_unb_cnt = 0;  //unbound read total for current variant
    
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
                
                if (line_cnt == 0) { //reading column header line (first line of each file)
                    
                    parse_TECdisplay_out_line(&line[v][0], &p_id[v], &p_vals[v], NULL, NULL, NULL, TDSPLY_HDR_LINE, 0);
                    
                } else { //reading data line
                    
                    bnd = 0; //initialize bound read count to zero
                    unb = 0; //initialize unbound read count to zero
                    frc = 0; //initialize fraction bound to zero
                    
                    parse_TECdisplay_out_line(&line[v][0], &p_id[v], &p_vals[v], &bnd, &unb, &frc, TDSPLY_DATA_LINE, 0);
                    
                    /* data lines are validated by comparing the variant id from each input file
                     to the variant id of the first input file. all variant ids must be identical
                     in order for the values files to be merged. */
                    if (strcmp(p_id[v], p_id[0])) {
                        printf("%s\n%s\n", p_id[v], p_id[0]);
                        printf("merge_replicates: error - variant ids are not aligned. aborting...\n");
                        abort();
                        
                    } else {
                        crrnt_bnd_cnt += bnd; //add reads to bound reads total
                        crrnt_unb_cnt += unb; //add reads to unbound reads total
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

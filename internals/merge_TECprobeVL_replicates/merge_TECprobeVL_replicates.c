//
//  merge_TECprobeVL_replicates.c
//  
//
//  Created by Eric Strobel on 1/19/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <getopt.h>
#include <dirent.h>

#include "../global/global_defs.h"
#include "../utils/io_management.h"

#include "../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"

#include "../mkmtrx/cotrans_mtrx.h"
#include "../mkmtrx/mkmtrx_defs.h"

#include "./merge_TECprobeVL_replicates_defs.h"
#include "./merge_TECprobeVL_replicates_structs.h"

#include "./read_analysis_directories.h"
#include "./make_output_directories.h"
#include "./merge_reactivity_profiles.h"
#include "./generate_sample_name.h"

void print_merge_record(sample_names * sn, output_files * outfiles, int * profiles_opened, SM2_analysis_directory * an_dir);

int main(int argc, char *argv[])
{
    char * prnt_dir_nm[MAX_RUNS] = {NULL}; //list of parent directories
        
    int dir_count = 0; //flag that input directory was supplied
    
    sample_names sn = {{{0}}}; //structure for storing input sample names and constructing merged name
    output_files outfiles;   //structure for storing output file pointers and names
    
    char out_dir_sffx[20] = {"merged_SM2_profiles"};
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"input",        required_argument,  0,  'i'}, //shapemapper analysis directory input
            {"out_dir_name", required_argument,  0,  'o'}, //output directory name
            {"sample-name",  required_argument,  0,  'n'}, //user-supplied sample name
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "i:o:n:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
                
            case 'i': //input directory
                if (dir_count < MAX_RUNS) { //if the number of input directories has not exceed the max
                    
                    //allocate memory for parent directory name storage
                    if ((prnt_dir_nm[dir_count] = calloc(strlen(argv[optind-1])+1, sizeof(*(prnt_dir_nm[dir_count])))) == NULL) {
                        printf("merge_replicate_SM2_out: error - memory allocation for directory name failed. aborting...\n");
                        abort();
                    }
                    strcpy(prnt_dir_nm[dir_count], argv[optind-1]); //store input directory name
                    dir_count++; //increment directory count
                    
                } else { //max input directories exceeded
                    printf("merge_replicate_SM2_out: error - too many input directories were provided. The maximum number allowed is %d. aborting...\n", MAX_RUNS);
                    abort();
                }
                
                break;
                
            case 'o': //output directory name
                
                //store output directory name
                if (snprintf(outfiles.out_dir, MAX_NAME, "%s_%s", argv[optind-1], out_dir_sffx) >= MAX_NAME) {
                    printf("merge_replicate_SM2_out: error - output directory name exceeds buffer. aborting...\n");
                    abort();
                }
                
                break;
                
            case 'n': //user-supplied output file name
                
                //store user-supplied output file name
                if (snprintf(sn.usr, MAX_NAME, "%s", argv[optind-1]) >= MAX_NAME) {
                    printf("merge_replicate_SM2_out: error - sample name exceeds buffer. aborting...\n");
                    abort();
                }
                
                break;
            
            default: printf("error: unrecognized option. Aborting program...\n"); abort();
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
    
    //check that input directory was provided
    if (dir_count < 2) {
        printf("merge_replicate_SM2_out: error - too few input directories were provided. aborting...\n");
        abort();
    }
    
    //if no output directory name as provided, set the output directory name
    if (!outfiles.out_dir[0]) {
        strcpy(outfiles.out_dir, out_dir_sffx);
    }
    
    //set the output sample name that will be used
    if (sn.usr[0]) {            //if a user-supplied output sample name was provided
        sn.sn2use = &sn.usr[0]; //use the user-supplied output sample name
    } else {
        sn.sn2use = &sn.mrg[0]; //otherwise, use the auto-generated sample name
    }
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int profiles_opened[MAX_RUNS]; //array to store the number of profiles opened for each input directory
    
    SM2_analysis_directory * an_dir = NULL; //pointer for SM2 analysis directory memory allocation
    
    //allocate memory for storing SM2 analysis directory information
    if ((an_dir = calloc(dir_count, sizeof(*an_dir))) == NULL) {
        printf("merge_replicate_SM2_out: error - memory allocation for directory/file pointers failed. aborting...\n");
        abort();
    }
    
    //initialize SM2 analysis directory structs
    for (i = 0; i < dir_count; i++) {                  //for each input directory
        strcpy(an_dir[i].prnt_dir_nm, prnt_dir_nm[i]); //store the parent directory name
        an_dir[i].min_tl = MIN_INIT;                   //initialize min_tl to MIN_INIT (512 (MAX_ROW) + 1)
        an_dir[i].max_tl = MAX_INIT;                   //initialize max_tl to MAX_INIT (0)
    }
        
    //read each SM2 analysis directory and confirm that the number of
    //files and transcript length range is the same for all samples
    for (i = 0; i < dir_count; i++) {
        
        //open all relevant directories and files in the parent directory
        profiles_opened[i] = read_prnt_directory(&an_dir[i], i, &sn);
        
        if (i) { //if processing a non-first parent directory
            
            //check that the number of profiles opened matches that of the first parent directory
            if (profiles_opened[i] != profiles_opened[0]) {
                printf("merge_replicate_SM2_out: error - the number of reactivity profiles opened is not the same for all input directories. aborting...\n");
                abort();
            }
            
            //check that the minimum transcript length matches that of the first parent directory
            if (an_dir[i].min_tl != an_dir[0].min_tl) {
                printf("merge_replicate_SM2_out: error - the minimum transcript length is not the same for all input directories. aborting...\n");
                abort();
            }
            
            //check that the maximum transcript length matches that of the first parent directory
            if (an_dir[i].max_tl != an_dir[0].max_tl) {
                printf("merge_replicate_SM2_out: error - the maximum transcript length is not the same for all input directories. aborting...\n");
                abort();
            }
        }
        
        //check the transcrip length contiguity of the opened profiles
        for (j = 0; j < MAX_ROW; j++) {
            
            //check that all transript lengths are within the min_tl/max_tl bounds
            if (an_dir[i].prf[j] != NULL && (j < an_dir[0].min_tl || j > an_dir[0].max_tl)) {
                printf("merge_replicate_SM2_out: error - out of bounds transcript length in %s data set. aborting...\n", an_dir[i].prnt_dir_nm);
            }
            
            //check that all transcript lengths within the min_tl/max_tl bounds have an associated profile
            if (an_dir[i].prf[j] == NULL && (j >= an_dir[0].min_tl && j <= an_dir[0].max_tl)) {
                printf("merge_replicate_SM2_out: error - missing transcript length %d  in %s data set. aborting...\n", j, an_dir[i].prnt_dir_nm);
            }
        }
    }
    
    
    //parse input sample names and generate merged sample name
    //an auto-generated sample name is always constructed because
    //sample metadata is validated during this process
    generate_sample_name(&sn);
    
    //print sample names to screen
    printf("\nuser-specified sample name: %s\n", (sn.usr[0]) ? sn.usr : "not provided");
    printf("auto-generated sample name: %s\n\n", sn.mrg);

    make_output_directories(&an_dir[0], &outfiles, &sn); //make output directories/open output files
    print_merge_record(&sn, &outfiles, &profiles_opened[0], &an_dir[0]); //generate merge record
    merge_profiles(&an_dir[0], dir_count, &outfiles);    //generate merged SM2 files
    
    //close output files
    for (i = an_dir[0].min_tl; i <= an_dir[0].max_tl; i++) {
        if (fclose(outfiles.ofp[i]) == EOF) {
            printf("read_SM2out_directory: error - failed to close output file. Aborting program...\n");
            abort();
        }
    }
    
    //close input files
    for (i = 0; i < dir_count; i++) {
        for (j = an_dir[0].min_tl; j <= an_dir[0].max_tl; j++) {
            if (fclose(an_dir[i].prf[j]) == EOF) {
                printf("read_SM2out_directory: error - failed to close input file. Aborting program...\n");
                abort();
            }
        }
    }
}

/* print_merge_record: print record of input files, merged sample
 name,and transcrip length inventory to screen and file*/
void print_merge_record(sample_names * sn, output_files * outfiles, int * profiles_opened, SM2_analysis_directory * an_dir)
{
    int i = 0; //general purpose index
    
    FILE * p_mrg_rcrd = NULL;         //merge record file pointer
    char mrg_rcrd_nm[MAX_LINE] = {0}; //merge record file name
    
    //generate merge record file name
    if (snprintf(mrg_rcrd_nm, MAX_LINE, "./%s/000_merge_record.txt", outfiles->out_dir) >= MAX_LINE) {
        printf("print_merge_record: error - merge record file name exceeded buffer. aborting...\n");
        abort();
    }
    
    //open merge record file
    if ((p_mrg_rcrd = fopen(mrg_rcrd_nm, "w")) == NULL) {
        printf("print_merge_record: error - could not open merge record output file. Aborting program...\n");
        abort();
    }
    
    //print input directory information
    for (i = 0; i < sn->cnt; i++) {
        printf("input %d: %s\n", i+1, sn->ipt[i]);
        fprintf(p_mrg_rcrd, "input %d: %s\n", i+1, sn->ipt[i]);
    }
    printf("merged:  %s\n\n", sn->sn2use);
    fprintf(p_mrg_rcrd, "merged:  %s\n\n", sn->sn2use);
    
    //print transcript length inventory
    printf("%d profiles opened per input directory\n", profiles_opened[0]);
    printf("minimum transcript length: %d\n", an_dir[0].min_tl);
    printf("maximum transcript length: %d\n\n", an_dir[0].max_tl);
    
    fprintf(p_mrg_rcrd, "%d profiles opened per input directory\n", profiles_opened[0]);
    fprintf(p_mrg_rcrd, "minimum transcript length: %d\n", an_dir[0].min_tl);
    fprintf(p_mrg_rcrd, "maximum transcript length: %d\n", an_dir[0].max_tl);
    
    //close merge record file
    if (fclose(p_mrg_rcrd) == EOF) {
        printf("print_merge_record: error - failed to close merge record file. Aborting program...\n");
        abort();
    }
}







//
//  mtrx2cols.c
//  
//
//  Created by Eric Strobel on 8/26/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <getopt.h>

#include "../global/global_defs.h"
#include "../utils/io_management.h"
#include "../mkmtrx/cotrans_mtrx.h"
#include "../seq_utils/parse_fasta.h"

#define DMS 0
#define SHAPE 1

/* structure declarations */


void make_fasta_file(cotrans_matrix * mtrx, char * seq, char * output_name, char * tl_out_nm, int len, int start_nt, /*int end_nt,*/ char * fasta_ipt);
void make_reactivity_file(cotrans_matrix * mtrx, char * output_name, char * tl_out_nm, int len, int start_nt, /*int end_nt,*/ char * reactivity_ipt);


int main(int argc, char *argv[])
{
    FILE * fp_fasta = NULL; //sequence fasta file pointer
    FILE * fp_mtrx = NULL;  //matrix file pointer
    
    char fa_name[MAX_LINE] = {0}; //name from fasta file
    char seq[MAX_LINE] = {0};     //input sequence
    int fasta_supplied = 0;       //flag that FASTA file was supplied
    
    cotrans_matrix mtrx = {{0}};  //cotrans_matrix struct
    int mtrx_cnt = 0;             //counts number of supplied matrices
    
    char RNAstructure_path[MAX_LINE] = {0}; //path to RNAstructure
    int path_supplied = 0;                  //counts number of supplied paths
    
    int probe = 0;          //probe
    int probe_supplied = 0; //counts number of supplied probe options
    
    char output_name[MAX_NAME] = {0}; //ouput file name
    
    int i = 0; //general purpose index
    int ret = 0; //variable for storing snprintf return value
    
    int start_nt = 0;
    //int end_nt = 0;
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"rna-structure", required_argument,  0,  'r'},  //RNAstructure path input
            {"matrix",        required_argument,  0,  'm'},  //matrix file input
            {"sequence",      required_argument,  0,  's'},  //sequence fasta file input
            {"probe",         required_argument,  0,  'p'},  //probe type
            {"output",        required_argument,  0,  'o'},  //output file name
            {"start-nt",      required_argument,  0,  'a'},  //nt to start at
            //{"end-nt",      required_argument,  0,  'b'},  //nt to end at
            {0, 0, 0, 0}
        };
        
        //if enabling end-nt option, short option code here
        c = getopt_long(argc, argv, "r:m:s:p:o:a:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
                
            case 'r': //path to RNAstructure supplied
                ret = snprintf(RNAstructure_path, MAX_LINE, "%s", argv[optind-1]);
                if (ret >= MAX_LINE || ret < 0) {
                    printf("average_TECprobe_matrices: error - error when storing RNAstructure path. aborting...\n");
                    abort();
                }
                path_supplied++;
                break;
                
            case 'm': //matrix supplied
                
                if (mtrx_cnt > 0) { //check if other matrices were supplied previously
                    printf("draw_intermediates: error - more than one reactivity matrix file was supplied. aborting...\n");
                    abort();
                }
                
                get_file(&fp_mtrx, argv[optind-1]);            //set file pointer to matrix file
                get_sample_name(argv[optind-1], &mtrx.fn[0]);  //set input file name
                strcpy(mtrx.nm, mtrx.fn);                      //set input file name as default matrix name
                store_mtrx(fp_mtrx, &mtrx);                     //store matrix in cotrans_mtrx struct
                
                mtrx_cnt++; //increment matrix counter
                break;
                
            case 's': //sequence supplied

                if (fasta_supplied) {
                    printf("draw_intermediates: error - more than one fasta file was supplied. aborting...\n");
                    abort();
                }
                
                get_file(&fp_fasta, argv[optind-1]); //set file pointer to sequence fasta file
                parse_fasta(fp_fasta, &fa_name[0], &seq[0]);
                
                fasta_supplied++;
                break;
                
            
            case 'p':
                if (probe_supplied) {
                    printf("draw_intermediates: error - too many probe options supplied. aborting...\n");
                    abort();
                    
                } else {
                    
                    if (!strcmp(argv[optind-1], "DMS")) {
                        probe = DMS;
                        
                    } else if (!strcmp(argv[optind-1], "SHAPE")) {
                        probe = SHAPE;
                    
                    } else {
                        printf("draw_intermediates: error - unrecognized probe. use 'DMS' or 'SHAPE'.aborting...\n");
                        abort();
                    }
                    
                    probe_supplied++;
                }
                break;
                
            case 'o': //output folder name supplied
                ret = snprintf(output_name, MAX_NAME, "%s", argv[optind-1]);
                if (ret >= MAX_NAME || ret < 0) {
                    printf("average_TECprobe_matrices: error - error when storing output name. aborting...\n");
                    abort();
                }
                break;
                
            case 'a': //start nt
                start_nt = strtol(argv[optind-1], NULL, 10)-1;
                break;
                
            /*case 'b': //end nt
                end_nt = strtol(argv[optind-1], NULL, 10);
                break;
             */
            
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
        
    if (path_supplied > 1) { //no more than one RNAstructure path can be provided
        printf("main: error - %d paths to RNAstructure were supplied. expected a maximum of 1. aborting...\n", path_supplied);
        abort();
    }
    
    if (fasta_supplied != 1) { //1 and only 1 fasta file must be provided
        printf("main: error - %d fasta files were supplied. expected 1. aborting...\n", fasta_supplied);
        abort();
    }
    
    if (mtrx_cnt != 1) { //at least 1 matrix file must be provided
        printf("main: error - %d reactivity matrix were supplied. expected 1. aborting...\n", mtrx_cnt);
        abort();
    }
    
    if (probe_supplied != 1) { //1 and only 1 probe type must be provided
        printf("main: error - %d probe types were supplied. expected 1. aborting...\n", probe_supplied);
        abort();
    }
    
    if (!output_name[0]) { //no output folder name was provided, set to default "out"
        strcpy(output_name, "out");
    }
    
    if (mtrx.col_cnt-1 != strlen(seq)) {
        printf("main: error - the longest reactivity string in the reactivity matrix (%d values) is not equal to the sequence length (%lu). aborting...", mtrx.col_cnt-1, strlen(seq));
    }
    
    if (start_nt) {
        if (start_nt < 0 || start_nt > mtrx.tl[MIN]) {
            printf("draw_intermediates: error - start nucleotide (%d) must be greater than 0 and less than the minimum transcript length of the current matrix (%d nt). aborting...", start_nt, mtrx.tl[MIN]);
            abort();
        }
    }
    
    /*if (end_nt) {
        if (end_nt < 0 || end_nt > mtrx.tl[MAX]) {
            printf("draw_intermediates: error - end nucleotide (%d) must be greater than 0 and less than the maximum transcript length of the current matrix (%d nt). aborting...", end_nt, mtrx.tl[MIN]);
            abort();
        }
    }
    
    if (start_nt && end_nt) {
        if (start_nt >= end_nt) {
            printf("draw_intermediates: error - start nucleotide (%d) must be less than end nucleotide (%d). aborting...", start_nt, end_nt);
            abort();
        }
    }*/
    
    mk_out_dir(output_name); //make parent output directory
        
    
    char tl_out_dir[MAX_LINE] = {0};
    char tl_out_nm[MAX_LINE] = {0};
    char command[MAX_LINE] = {0};
    char fasta_ipt[MAX_LINE] = {0};
    char reactivity_ipt[MAX_LINE] = {0};
    
    char * probe_str = NULL;
    char dms_str[6] = {"--DMS"};
    char shape_str[8] = {"--SHAPE"};
    
    if (probe == DMS) {
        probe_str = &dms_str[0];
    } else if (probe == SHAPE) {
        probe_str = &shape_str[0];
    } else {
        printf("draw_intermediates: error - unexpected probe type ahead of command construction. this error should not be reachable. aborting...\n");
        abort();
    }
    
    for (i = mtrx.tl[MIN]; i <= mtrx.tl[MAX]; i++) {
        
        //make transcript output directory
        ret = snprintf(tl_out_nm, MAX_LINE, "%s_%dnt", mtrx.nm, i-start_nt); //TODO add error check
        ret = snprintf(tl_out_dir, MAX_LINE, "%s/%s", output_name, tl_out_nm); //TODO add error check
        mk_out_dir(tl_out_dir);
        
        make_fasta_file(&mtrx, seq, output_name, tl_out_nm, i, start_nt, /*end_nt,*/ fasta_ipt);
        make_reactivity_file(&mtrx, output_name, tl_out_nm, i, start_nt, /*end_nt,*/ reactivity_ipt);
        
        //TODO add error check
        ret = snprintf(command, MAX_LINE, "%s/exe/Fold %s %s/%s.ct %s %s",
                       RNAstructure_path,
                       fasta_ipt,
                       tl_out_dir,
                       tl_out_nm,
                       probe_str,
                       reactivity_ipt);
        
        system(command);
        
        ret = snprintf(command, MAX_LINE, "%s/exe/draw %s/%s.ct %s/%s.ps -s %s",
                       RNAstructure_path,
                       tl_out_dir,
                       tl_out_nm,
                       tl_out_dir,
                       tl_out_nm,
                       reactivity_ipt);
        
        system(command);
        
    }
}
  

void make_fasta_file(cotrans_matrix * mtrx, char * seq, char * output_name, char * tl_out_nm, int len,  int start_nt, /*int end_nt,*/ char * fasta_ipt)
{
    int i = 0;            //general purpose index
    int ret = 0;          //variable for storing snprintf return value
    FILE * out_fp = NULL;
    
    //make fasta
    //open output file
    ret = snprintf(fasta_ipt, MAX_LINE, "%s/%s/%s.fa", output_name, tl_out_nm, tl_out_nm);
    if (ret >= MAX_LINE || ret < 0) {
        printf("make_fasta_file: error - error when constructing fasta input file name. aborting...\n");
        abort();
    }
    
    if ((out_fp = fopen(fasta_ipt, "w")) == NULL) {
        fprintf(out_fp, "draw_intermediates: error - failed to open output file. Aborting program...\n");
        abort();
    }
    
    fprintf(out_fp, ">%s_%dnt\n", mtrx->nm, len);
    
    
    
    for (i = start_nt; i < len; i++) {
        if (i < len-11) {
            fprintf(out_fp, "%c", seq[i]);
        } else {
            fprintf(out_fp, "%c", tolower(seq[i]));
        }
    }
    fprintf(out_fp, "\n");
    
    //close output file
    if (fclose(out_fp) == EOF) {
        fprintf(out_fp, "draw_intermediates: error - failed to close output file. Aborting program...\n");
        abort();
    }
    
    
}

void make_reactivity_file(cotrans_matrix * mtrx, char * output_name, char * tl_out_nm, int len,  int start_nt, /*int end_nt,*/ char * reactivity_ipt)
{
    int i = 0;            //general purpose index
    int ret = 0;          //variable for storing snprintf return value
    
    FILE * out_fp = NULL;
    
    //make fasta
    //open output file
    ret = snprintf(reactivity_ipt, MAX_LINE, "%s/%s/%s_reactivity.txt", output_name, tl_out_nm, tl_out_nm);
    if (ret >= MAX_LINE || ret < 0) {
        printf("make_reactivity_file: error - error when constructing reactivity input file name. aborting...\n");
        abort();
    }
    
    if ((out_fp = fopen(reactivity_ipt, "w")) == NULL) {
        fprintf(out_fp, "draw_intermediates: error - failed to open output file. Aborting program...\n");
        abort();
    }
    
    
    for (i = start_nt + 1; i <= len; i++) {
        fprintf(out_fp, "%d\t%s\n", i-start_nt, mtrx->vals[len-mtrx->tl[MIN]+1][i]);
    }
    fprintf(out_fp, "\n");
    
    //close output file
    if (fclose(out_fp) == EOF) {
        fprintf(out_fp, "draw_intermediates: error - failed to close output file. Aborting program...\n");
        abort();
    }
}

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

void clean_up_output(int layr_cnt, constraint_metadata * cons_meta)
{
    int i = 0;
    
    char rm_mrgd_out_command[2048] = {0};
    char mv_outFiles_command[2048] = {0};
    char rm_out_dirs_command[2048] = {0};
    
    char crnt_layr[2048] = {0};
    
    system("rm merged_out.txt");
    
    for (i = 0; i < layr_cnt; i++) {
        sprintf(crnt_layr, "layer%d", i+1);
        chdir(crnt_layr);
        
        sprintf(rm_mrgd_out_command, "rm -r */*merged_out.txt");
        system(rm_mrgd_out_command);
        
        if (i == 0) {
            sprintf(mv_outFiles_command, "mv %s/*.txt .", cons_meta[i].sn);
            sprintf(rm_out_dirs_command, "rm -r %s", cons_meta[i].sn);
        } else {
            sprintf(mv_outFiles_command, "mv *_%s/*.txt .", cons_meta[i].sn);
            sprintf(rm_out_dirs_command, "rm -r *_%s", cons_meta[i].sn);
        }
        system(mv_outFiles_command);
        system(rm_out_dirs_command);
        
        chdir("..");
        
        
    }
}





void aggregate_output(int vals_cnt, values_input * vals, int layr_cnt, constraint_metadata * cons_meta, char * out_prefix, char * col_id)
{
    extern output_file_names out_fns[MAX_LAYERS];
    
    FILE ** TECDnav_out = NULL;
    
    FILE ** ofp = NULL;
    char out_nm[2048] = {0};
    char tmp_out_nm[2048] = {0};
    
    char tmp_val[128] = {0};
    
    int EOF_tot = 0;
    int * EOF_rchd = NULL;
    
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    
    int v = 0;
    
    int cols2incld[MAX_VALS] = {0};
    
    char layr_dir[8] = {0};
    
    //strcpy(out_nm, out_prefix);
    
    for (i = 0; i < layr_cnt; i++) {
        
        EOF_tot = 0;
        
        //go to output directory
        sprintf(layr_dir, "layer%d", i+1);
        printf("%s\n", layr_dir);
        chdir(layr_dir);
        system("pwd");
        
        //open files to be aggregated
        //allocate file pointers
        if ((TECDnav_out = calloc(out_fns[i].f_cnt, sizeof(*(TECDnav_out)))) == NULL) {
            printf("aggreggate_output: error - memory allocation for file pointers failed. aborting...\n");
            abort();
        }
        
        if ((EOF_rchd = calloc(out_fns[i].f_cnt, sizeof(*(EOF_rchd)))) == NULL) {
            printf("aggreggate_output: error - memory allocation for EOF detection failed. aborting...\n");
            abort();
        }
        
        for (j = 0; j < out_fns[i].f_cnt; j++) {
            get_file(&TECDnav_out[j], out_fns[i].fn[j]);
        }
        
        //generate aggregate output file name and open aggregate output file
        if ((ofp = calloc(vals_cnt, sizeof(*(ofp)))) == NULL) {
            printf("aggreggate_output: error - memory allocation for output file pointers failed. aborting...\n");
            abort();
        }
        
        if (!i) {
            strcpy(out_nm, cons_meta[i].sn);
        } else {
            sprintf(out_nm, "%s_%s", out_nm, cons_meta[i].sn);
        }
        
        for (j = 0; j < vals_cnt; j++) {
            
            sprintf(tmp_out_nm, "%s_%s.txt", vals[j].nm, out_nm);
            
            if ((ofp[j] = fopen(tmp_out_nm, "w")) == NULL) {
                printf("aggregate_output: error - could not open output file. Aborting program...\n");
                abort();
            }
        }
        
        read_output_hdrs(vals_cnt, col_id, i, TECDnav_out, &cols2incld[0], ofp);
        
        /*for (v = 0; v < vals_cnt; v++) {
            printf("%d\n", cols2incld[v]);
            
        }*/
        
        while (EOF_tot < out_fns[i].f_cnt) {
            read_data_line(vals_cnt, i, TECDnav_out, &cols2incld[0], &EOF_rchd[0], &EOF_tot, ofp);
            //printf("EOF_tot = %d of %d\n", EOF_tot, out_fns[i].f_cnt);
        }
        
        //printf("done processing layer %d\n", i);
        
        //close files that were aggregated
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
        
        free(TECDnav_out);
        free(EOF_rchd);
        chdir("..");
        
    }
}

void read_output_hdrs(int vals_cnt, char * col_id, int crnt_layr, FILE ** TECDnav_out, int * cols2incld, FILE ** ofp)
{
    
    extern output_file_names out_fns[MAX_LAYERS];
    
    int i = 0;
    int j = 0;
    int k = 0;
    
    int v = 0;
    int o = 0;
    
    int col = 0;
    int line_end = 0;
    
    char line[MAX_LINE] = {0};
    char tmp_val[MAX_LINE] = {0};
    
    char * p_col_id = NULL;
    
    
    //add check that number of cols to include matches values file input counts
    for (i = 0; i < out_fns[crnt_layr].f_cnt; i++) {
        
        v = 0;
        col = 0;
        
        if (get_line(line, TECDnav_out[i])) {
            
            tmp_val[0] = '\0';
            for (line_end = 0, j = 0, k = 0; !line_end; j++) {
                
                if (line[j] != '\t' && line[j]) {
                    tmp_val[k++] = line[j];
                } else {
                    if (!line[j]) {
                        line_end = 1;
                    }
                    
                    tmp_val[k] ='\0';
                    k = 0;
                    
                    if ((p_col_id = strstr(tmp_val, col_id)) != NULL) {
                        
                        if (p_col_id[-1] == '_' && !p_col_id[strlen(p_col_id)]) {

                            if (!i) {
                                cols2incld[v] = col;
                                fprintf(ofp[v++], "%s", tmp_val);
                            } else if (cols2incld[v] == col) {
                                fprintf(ofp[v++], "\t%s", tmp_val);
                            } else {
                                printf("read_output_hdrs: error - columns to include are not the same in all files. aborting...\n");
                                abort();
                            }
                        }
                    }
                    
                    col++;
                }
            }
        } else {
            //TODO: add file name to message
            printf("read_output_hdrs: error - output file did not contain any headers. aborting...\n");
        }
        
    }
    
    for (v = 0; v < vals_cnt; v++) {
        fprintf(ofp[v], "\n");
    }
    
}


void read_data_line(int vals_cnt, int crnt_layr, FILE ** TECDnav_out, int * cols2incld, int * EOF_rchd, int * EOF_tot, FILE ** ofp)
{
    extern output_file_names out_fns[MAX_LAYERS];
    
    int i = 0;
    int j = 0;
    int k = 0;
    
    int v = 0;
    
    int col = 0;
    int line_end = 0;
    
    int printed_data = 0;
    
    char line[MAX_LINE] = {0};
    char tmp_val[MAX_LINE] = {0};
    
    char * p_col_id = NULL;
    
    
    
    //add check that number of cols to include matches values file input counts
    for (i = 0; i < out_fns[crnt_layr].f_cnt; i++) {
        
        v = 0;
        col = 0;
        
        if (!EOF_rchd[i]) {
            if (get_line(line, TECDnav_out[i])) {
                
                tmp_val[0] = '\0';
                for (line_end = 0, j = 0, k = 0; !line_end; j++) {
                    
                    if (line[j] != '\t' && line[j]) {
                        tmp_val[k++] = line[j];
                    } else {
                        if (!line[j]) {
                            line_end = 1;
                        }
                        
                        tmp_val[k] ='\0';
                        k = 0;
                        
                        if (col == cols2incld[v]) {
                            if (!i) {
                                fprintf(ofp[v++], "%s", tmp_val);
                            } else {
                                fprintf(ofp[v++], "\t%s", tmp_val);
                            }
                            printed_data = 1;
                        }
                        
                        /*
                        
                        //printf("%s\n", tmp_val);
                        
                        if ((p_col_id = strstr(tmp_val, col_id)) != NULL) {
                            if (p_col_id[-1] == '_' && !p_col_id[strlen(p_col_id)]) {
                                //printf("%d %s\n", v, tmp_val);
                                
                                if (!i) {
                                    cols2incld[v++] = col;
                                } else if (cols2incld[v++] != col) {
                                    printf("read_output_hdrs: error - columns to include are not the same in all files. aborting...\n");
                                    abort();
                                }
                            }
                        }*/
                        
                        col++;
                    }
                }
            } else {
                //abort();
                EOF_rchd[i] = 1;
                (*EOF_tot)++;
                
                for (v = 0; v < vals_cnt; v++) {
                    if (i) {
                        fprintf(ofp[v], "\t");
                    }
                }
            }
        } else {
            for (v = 0; v < vals_cnt; v++) {
                if (i) {
                    fprintf(ofp[v], "\t");
                }
            }
        }
    }
    
    /*if (!printed_data) {
        printf("deleting line (%d)\n", *EOF_tot);
        for (v = 0; v < vals_cnt; v++) {
            for (j = 0; j < *EOF_tot; j++) {
                printf("printing backspace\n");
                fprintf(ofp[v], "\b\0\b");
            }
            fflush(ofp[v]);
        }
    }*/
    
    for (v = 0; v < vals_cnt; v++) {
        fprintf(ofp[v], "\n");
    }

}

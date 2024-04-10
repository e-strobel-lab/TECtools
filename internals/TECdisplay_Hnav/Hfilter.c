//
//  Hfilter.c
//  
//
//  Created by Eric Strobel on 11/7/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "../utils/io_management.h"

#include "../seq_utils/basemap.h"

#include "./TECdisplay_Hnav_defs.h"
#include "./TECdisplay_Hnav_structs.h"
#include "./TECdisplay_Hnav_global_vars.h"

#include "Hfilter.h"

/* Hfilter: hierarchically filter input data using a series of constraints files*/
void Hfilter(constraint_metadata * cons_meta, char * prev_out_prfx, int cl, int layr_cnt, char * prev_TECDnav_out, char * layr_list[MAX_LAYERS])
{
    extern char TDHN_TECDnav_path[MAX_LINE+1]; //TECdisplay_navigator path
    extern char TDHN_merged_out_nm[15];      //merged_output filename
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    char tmp_str[MAX_LINE+1] = {0}; //array for temporary string storage
    
    char cl_dir[MAX_NAME+1] = {0};  //current layer directory
    
    char TECDnav_command[MAX_LINE+1] = {0};  //TECdisplay_navigator command
    char TECDnav_out_dir[MAX_LINE+1] = {0};  //TECdisplay_navigator output directory name
    char TECDnav_ipt_data[MAX_LINE+1] = {0}; //TECdisplay_navigator input data file name
    char crnt_out_prfx[MAX_LINE+1] = {0};    //output file prefix for current TECdisplay_navigator run
    char *prfx2use = NULL;                   //pointer to output file prefix to be used for TECdisplay_navigator run
    char crnt_TECDnav_out[MAX_LINE+1] = {0}; //path to current TECdisplay_navigator out directory, passed to next Hfilter call

    int prev_out_files = 0; //number of output files generated by the previous layer
    
    char out_msg_dir[MAX_LINE+1] = {0};           //directory for storing output message files
    char out_msg_dir_name[13] = {"out_messages"}; //ouput message directory name
    struct stat st = {0};                         //storage for stat output
    
    int ret = 0; //variable for storing snprintf return value
    
    sprintf(cl_dir, "layer%d", cl+1); //generate directory name for the current layer
    
    if (stat(cl_dir, &st) == -1  && cl < layr_cnt) { //if output directory for current layer has not yet been made
        mk_out_dir(cl_dir);                                      //make output directory for current layer
        sprintf(out_msg_dir, "%s/%s", cl_dir, out_msg_dir_name); //construct output message directory name
        mk_out_dir(out_msg_dir);                                 //make output message directory
    }
    
    //set number of prevous output files
    if (cl == 0) {          //if at layer 0
        prev_out_files = 1; //set number of previous output files to 1
        
    } else if (cl > 0 && cons_meta[cl-1].typ == 'x') {  //if at layer >0 and previous constraint excluded matches
        prev_out_files = 1;                             //set number of previous output files to 1
        
    } else if (cons_meta[cl-1].typ == 'c') {    //if previous constraint included matches
        prev_out_files = cons_meta[cl-1].c_cnt; //set number of previous output files to number of previous constraints
        
    } else {
        printf("Hfilter: error - this should be unreachable. aborting...\n");
        abort();
    }
    
    if (cl < layr_cnt) { //if the bottom layer has not yet been processed
        
        for (i = 0; i < prev_out_files; i++) { //for all previous output files
            
            chdir(cl_dir); //change to output directory for current layer
            
            if (cl == 0) { //if processing the first layer
                strcpy(TECDnav_out_dir, cons_meta[cl].sn);    //set the TECDnav out dir name to the first constraint sample name
                strcpy(TECDnav_ipt_data, TDHN_merged_out_nm); //set the input data file to the merged input file
                prfx2use = prev_out_prfx;                     //set the output file prefix to the user-supplied prefix
                
            } else { //non-first layer
                //set input file information for the previous and current layer in layr_list
                //for the previous layer:
                //  if the previous constraint included matches, use the constraint name of the next input data file
                //  if the previous constraint excluded matches, use the sample name of previous constraint
                //for the current layer, use the sample name of the constraint file
                layr_list[cl-1] = (cons_meta[cl-1].typ == 'c') ? &cons_meta[cl-1].cn[i][0] : &cons_meta[cl-1].sn[0];
                layr_list[cl] = &cons_meta[cl].sn[0];
                
                //construct the TECdisplay_navigator output directory name
                strcpy(TECDnav_out_dir, layr_list[0]); //copy the first layer file info to TECDnav_out_dir
                for (j = 1; j <= cl ; j++) {           //append the remaining layer(s') file info to TECDnav_out_dir
                    
                    //append layer name to TECDnav_out_dir
                    strcpy(tmp_str, TECDnav_out_dir);
                    ret = snprintf(TECDnav_out_dir, MAX_LINE, "%s_%s",  tmp_str, layr_list[j]);
                    if (ret >= MAX_LINE || ret < 0) {
                        printf("Hfilter: error - error when constructing TECdisplay_navigator output directory name. aborting...\n");
                        abort();
                    }
                }
                
                //append the info for the next input data file to the previous output file prefix to generate current out prefix
                ret = snprintf(crnt_out_prfx, MAX_LINE, "%s_%s", prev_out_prfx, layr_list[cl-1]);
                if (ret >= MAX_LINE || ret < 0) {
                    printf("Hfilter: error - error when constructing TECdisplay_navigator output prefix. aborting...\n");
                    abort();
                }
                prfx2use = crnt_out_prfx; //set crnt_out_prfx as the prefix to use in the TECdisplay_navigator command
                
                //construct input data file name
                ret = snprintf(TECDnav_ipt_data, MAX_LINE, "%s/%s.txt", prev_TECDnav_out, crnt_out_prfx);
                if (ret >= MAX_LINE || ret < 0) {
                    printf("Hfilter: error - error when constructing input data file name. aborting...\n");
                    abort();
                }
            }
            
            //generate a TECdisplay navigator command
            ret = snprintf(TECDnav_command, MAX_LINE, "%s -v ../%s -c ../../%s%s-n -f %s -o %s > %s/%s_out_msg.txt",
                            TDHN_TECDnav_path,  //path to TECdisplay_navigator
                            TECDnav_ipt_data,   //input data
                            cons_meta[cl].fn,   //constraint filename
                            (cons_meta[cl].typ == 'x') ? " -x " : " ", //exclusion option
                            prfx2use,           //prefix for output files
                            TECDnav_out_dir,    //TECdisplay_navigator output directory
                            out_msg_dir_name,   //output message directory name
                            TECDnav_out_dir);   //prefix for output message file
            if (ret >= MAX_LINE || ret < 0) {
                printf("Hfilter: error - error when constructing a TECdisplay_navigator command. aborting...\n");
                abort();
            }
            
            system(TECDnav_command); //run the TECdisplay_Navigator
            
            store_out_filenames(cl, prfx2use, cons_meta); //store output file names
            
            
            //construct the path to the current TECDnav output directory
            //this is needed so that the next layer can find the correct
            //input files
            ret = snprintf(crnt_TECDnav_out, MAX_LINE, "%s/%s", cl_dir, TECDnav_out_dir);
            if (ret >= MAX_LINE || ret < 0) {
                printf("Hfilter: error - error when constructing current TECdisplay_navigator output directory. aborting...\n");
                abort();
            }
            
            chdir(".."); //move up one directory
            Hfilter(cons_meta, prfx2use, cl+1, layr_cnt, crnt_TECDnav_out, layr_list); //process the next layer
        }
    }
    
    return;
}

/* store_out_filenames: store output file name in output_file_names structure */
void store_out_filenames(int cl, char * prefix, constraint_metadata * cons_meta)
{    
    extern output_file_names out_fns[MAX_LAYERS]; //storage for output file names of each layer
    
    int i = 0;            //general purpose index
    int out_file_cnt = 0; //number of output files
    
    char tmp_out_fn[MAX_LINE+1] = {0}; //temp storage for generating output file names
    char * str2appnd = NULL;         //pointer to string that will be appended to the current prefix
    
    int ret = 0; //variable for storing snprintf return value
    
    //set output file count
    //if the current layer constraints file includes matches, set to the number of constraints in the file
    //if the current layer constraints file excludes matches, set to 1
    out_file_cnt = (cons_meta[cl].typ == 'c') ? cons_meta[cl].c_cnt : 1;
    
    //store output file names
    for (i = 0; i < out_file_cnt; i++) {
        
        //check that index of the next stored file name does not exceed the expected number of output file names
        if (out_fns[cl].nxt >= out_fns[cl].f_cnt) {
            printf("store_out_filenames: error - making too many output file names for a constraint\n");
        }
        
        //set the string that will be appended to the current prefix, then append to prefix
        //if the current layer constraint file includes matches, set to the current constraint name
        //if the current layer constraint file excludes matches, set to the constraints file sample name
        str2appnd = (cons_meta[cl].typ == 'c') ? cons_meta[cl].cn[i] : cons_meta[cl].sn;
        ret = snprintf(tmp_out_fn, MAX_LINE, "%s_%s.txt", prefix, str2appnd);
        if (ret >= MAX_LINE || ret < 0) {
            printf("store_out_filenames: error - error when constructing output file name. aborting...\n");
            abort();
        }
        
        //allocate memory for output file name
        if ((out_fns[cl].fn[out_fns[cl].nxt] = malloc((strlen(tmp_out_fn)+1) * sizeof(*(out_fns[cl].fn[out_fns[cl].nxt])))) == NULL) {
            printf("Hfilter: error - memory allocation for output file name failed. aborting...\n");
            abort();
        }
        
        strcpy(out_fns[cl].fn[out_fns[cl].nxt++], tmp_out_fn); //store the output file name
    }
    
    return;
}

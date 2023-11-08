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


void Hfilter(constraint_metadata * cons_meta, char * prev_out_prfx, int cl, int layr_cnt, char * prev_TECDnav_out, char * layr_list[MAX_LAYERS])
{
    extern char TDHN_TECDnav_path[2048];      //TECdisplay_navigator path
    extern char TDHN_merged_out_nm[15];
    
    int i = 0;
    int j = 0;
    
    char pl_dir[256] = {0};
    char cl_dir[256] = {0};
    
    char TECDnav_command[2048] = {0};
    char TECDnav_out_dir[2048] = {0};
    char crnt_TECDnav_out[2048] = {0};
    char crnt_out_prfx[2048] = {0};
    char *prfx2use = NULL;
    char ipt_data_file[2048] = {0};
    int  prev_out_files = 0;
    char prev_layr_type = 0;
    
    char out_msg_dir[2048] = {0};
    char out_msg_dir_name[13] = {"out_messages"};
    struct stat st = {0};
    
    
    
    sprintf(cl_dir, "layer%d", cl+1);
    
    if (stat(cl_dir, &st) == -1  && cl < layr_cnt) {             //if output directory for current layer has not yet been made
        
        //make output directories
        mk_out_dir(cl_dir);                                      //make output directory for current layer
        sprintf(out_msg_dir, "%s/%s", cl_dir, out_msg_dir_name); //construct output message directory name
        mk_out_dir(out_msg_dir);                                 //make output message directory
        
        //print constraint metadata for current layer
        printf("%s\n", cons_meta[cl].fn);                        //print constraint filename
        printf("  sampN: %s\n", cons_meta[cl].sn);               //print constraint sample name
        printf("  type:  %c\n", cons_meta[cl].typ);              //print constraint type
        printf("  count: %d\n", cons_meta[cl].c_cnt);              //print constraint count
        for (i = 0; i < cons_meta[cl].c_cnt; i++) {                //print constraint names
            printf("  %2d: %s\n", i, cons_meta[cl].cn[i]);
        }
        printf("\n");
    }
    
    if (cl == 0 || cons_meta[cl-1].typ == 'x') {
        prev_out_files = 1;
    } else if (cons_meta[cl-1].typ == 'c') {
        prev_out_files = cons_meta[cl-1].c_cnt;
    } else {
        ;
    }
    
    if (cl < layr_cnt) {
        
        for (i = 0; i < prev_out_files; i++) { //for all previous output files
            
            chdir(cl_dir); //change to output directory for current layer
            
            if (cl == 0) {
                sprintf(TECDnav_out_dir, "%s", cons_meta[cl].sn); //construct the TECDnav output directory name
                sprintf(ipt_data_file, "%s",TDHN_merged_out_nm);  //set the input data file to the merged input files
                prfx2use = prev_out_prfx;                         //set the output file prefix to the user-supplied prefix
                
            } else {
                //set input file information for the previous and current layer in layr_list
                //for the previous layer:
                //  if the previous constraint included matches, use the constraint name of the next input data file
                //  if the previous constraint excluded matches, use the sample name of previous constraint
                //for the current layer, use the sample name of the constraint file
                layr_list[cl-1] = (cons_meta[cl-1].typ == 'c') ? &cons_meta[cl-1].cn[i][0] : &cons_meta[cl-1].sn[0];
                layr_list[cl] = &cons_meta[cl].sn[0];
                
                //construct the TECdisplay_navigator output directory name
                strcpy(TECDnav_out_dir, layr_list[0]); //copy the first layer file info to TECDnav_out_dir
                for (j = 1; j <=cl ; j++) {            //append the remaining layer(s) file info to TECDnav_out_dir
                    sprintf(TECDnav_out_dir, "%s_%s",  TECDnav_out_dir, layr_list[j]);
                }
                
                //append the info for the next input data file to the current output file prefix
                sprintf(crnt_out_prfx, "%s_%s", prev_out_prfx, (cons_meta[cl-1].typ == 'c') ? cons_meta[cl-1].cn[i] : cons_meta[cl-1].sn);
                prfx2use = crnt_out_prfx;
                
                //construct input data file name
                sprintf(ipt_data_file, "%s/%s.txt", prev_TECDnav_out, crnt_out_prfx);
            }
            
            //generate a TECdisplay navigator command
            sprintf(TECDnav_command, "%s -v ../%s -c ../../%s -n -f %s -o %s",
                    TDHN_TECDnav_path,  //path to TECdisplay_navigator
                    ipt_data_file,      //input data
                    cons_meta[cl].fn,   //constraint filename
                    prfx2use,           //prefix for output files
                    TECDnav_out_dir);   //TECdisplay_navigator output directory
            
            if (cons_meta[cl].typ == 'x') {     //if running exclude matches mode
                strcat(TECDnav_command, " -x"); //append -x option
            }
            
            //append output message location
            sprintf(TECDnav_command, "%s > %s/%s_out_msg.txt", TECDnav_command, out_msg_dir_name, TECDnav_out_dir);
            
            system(TECDnav_command); //run the TECdisplay_Navigator
            
            store_out_filenames(cl, prfx2use, cons_meta);
            
            
            //construct the path to the current TECDnav output directory
            //this is needed so that the next layer can find the correct
            //input files
            sprintf(crnt_TECDnav_out, "%s/%s", cl_dir, TECDnav_out_dir);
            
            chdir(".."); //move up one directory
            Hfilter(cons_meta, prfx2use, cl+1, layr_cnt, crnt_TECDnav_out, layr_list); //process the next layer
            

        }
    }
    
    return;
    
}


void store_out_filenames(int cl, char * prefix, constraint_metadata * cons_meta)
{
    //TODO add test that not exceeding file name array bounds
    
    extern output_file_names out_fns[MAX_LAYERS];
    
    int i = 0;
    int out_file_cnt = 0;
    
    char tmp_out_fn[2048] = {0};
    char * str2appnd = NULL;
    
    out_file_cnt = (cons_meta[cl].typ == 'c') ? cons_meta[cl].c_cnt : 1;

    
    for (i = 0; i < out_file_cnt; i++) {
        
        str2appnd = (cons_meta[cl].typ == 'c') ? cons_meta[cl].cn[i] : cons_meta[cl].sn;
        
        if (out_fns[cl].nxt > out_fns[cl].f_cnt) {
            printf("error - making too many output file names for a constraint\n");
        }
        
        sprintf(tmp_out_fn, "%s_%s.txt", prefix, str2appnd);
        
        //allocate memory for wt name
        if ((out_fns[cl].fn[out_fns[cl].nxt] = malloc((strlen(tmp_out_fn)+1) * sizeof(*(out_fns[cl].fn[out_fns[cl].nxt])))) == NULL) {
            printf("Hfilter: error - memory allocation for file name failed. aborting...\n");
            abort();
        }
        
        strcpy(out_fns[cl].fn[out_fns[cl].nxt++], tmp_out_fn);
        
    }
}

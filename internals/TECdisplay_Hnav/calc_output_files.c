//
//  calc_output_files.c
//  
//
//  Created by Eric Strobel on 11/7/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "./TECdisplay_Hnav_defs.h"
#include "./TECdisplay_Hnav_structs.h"
#include "./TECdisplay_Hnav_global_vars.h"

#include "calc_output_files.h"

int calc_output_files(int layr_cnt, constraint_metadata * cons_meta)
{
    extern output_file_names out_fns[MAX_LAYERS];
    
    int i = 0;
    int j = 0;
    
    int f_cnt = 1;
    int tot_files = 0;
    
    printf("the provided constraints will generate the following files:\n\n");
    
    for (i = 0; i < layr_cnt; i++) {
        if (cons_meta[i].typ == 'c') {
            f_cnt *= cons_meta[i].c_cnt;
        }
        
        out_fns[i].f_cnt = f_cnt;
        tot_files += f_cnt;
        printf("layer %d: %d files\n", i+1, f_cnt);
        
        //allocate pointers for file names
        if ((out_fns[i].fn = calloc(f_cnt, sizeof(*(out_fns[i].fn)))) == NULL) {
            printf("calc_output_files: error - memory allocation output file names failed. aborting...\n");
            abort();
        }
    }
    
    printf("in total, %d files will be generated.\n\n", tot_files);
    
    
    
    return tot_files;
}

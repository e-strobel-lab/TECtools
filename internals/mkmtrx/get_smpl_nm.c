//
//  get_smpl_nm.c
//  
//
//  Created by Eric Strobel on 10/11/22.
//

#include "get_smpl_nm.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"

/* obtain sample name from prnt_dir input. expected parent directory is:
 <sample_name>/<analysis_name>, which is converted to:
 <sample_name>_analysis_name> */
int get_smpl_nm(char * prnt_dir, char * smpl_nm) {
    
    int i = 0; //general purpose index
    
    char prnt_dir_cpy[MAX_LINE] = {0}; //array to store copy of parent directory
    char * smpl_nm_ptr = NULL;         //pointer to processed sample name string
    
    /* extract sample name from parent directory name */
    strcpy(prnt_dir_cpy, prnt_dir); //copy parent directory, this allows '/' trimming without changing prnt_dir input
    
    int last_indx = strlen(prnt_dir_cpy);
    
    if (prnt_dir_cpy[last_indx-1] == '/') { //if the input parent directory input string ends in a forward slash
        prnt_dir_cpy[last_indx-1] = '\0';   //remove the terminal forward slash from prnt_dir_cpy
        last_indx--;
    }
    
    int found_start = 0;
    for (i = 0; prnt_dir_cpy[i]; i++) {
        //skip past non-directory name characters (e.g. '.' and '/') to find the start of the sample name
        if ((isalpha(prnt_dir_cpy[i]) || isdigit(prnt_dir_cpy[i])) && !found_start) {
            smpl_nm_ptr = &prnt_dir_cpy[i];   //set pointer to the start of the sample name
            found_start = 1;                  //set flag that the start of the sample name was found
        } else if (prnt_dir_cpy[i] == '/') {  //replace '/' characters with '_'
            prnt_dir_cpy[i] = '_';
        }
    }
    
    strcpy(smpl_nm, smpl_nm_ptr); //copy sample name to smpl_nm array
    
    printf("sample name: %s\n", smpl_nm); //print sample name
    
    return 1;
}

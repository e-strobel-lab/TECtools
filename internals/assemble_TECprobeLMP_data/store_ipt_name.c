//
//  store_ipt_name.c
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#include "../global/global_defs.h"
#include "./assemble_TECprobeLMP_data_defs.h"
#include "./assemble_TECprobeLMP_data_structs.h"

#include "store_ipt_name.h"

/* store_ipt_name: allocate memory for input file name storage
 and store input file name */
void store_ipt_name(char ** storage, char * ipt_nm)
{
    //allocate memory for input name
    if (((*storage) = malloc((strlen(ipt_nm)+1) * sizeof(**(storage)))) == NULL) {
        printf("store_ipt_name: error - memory allocation for input name storage failed. aborting...\n");
        abort();
    }
    
    strcpy(*storage, ipt_nm); //copy input name to storage location
}

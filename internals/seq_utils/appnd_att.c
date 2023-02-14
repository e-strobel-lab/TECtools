//
//  appnd_att.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../global/global_defs.h"

#include "appnd_att.h"

/* appnd_att: append attribute code to end of read id in line 1 of fastq file */
int appnd_att(char * read1, char * read2, char * att)
{
    char new_rd_name[READ_MAX][(MAX_LINE*2)] = {{0}}; //array for new read name
    
    char *rd_ptr[READ_MAX] = {NULL};    //pointers to read 1 and 2 to use in for loop below
    char *prefix_ptr = NULL;            //pointer to new read name prefix
    char *suffix_ptr = NULL;            //pointer to new read name suffix
    
    rd_ptr[READ1] = &(read1[0]);    //set read 1 pointer
    rd_ptr[READ2] = &(read2[0]);    //set read 2 pointer
    
    int i = 0;
    int j = 0;
    
    for (i = 0; i < READ_MAX; i++) {    //perform operation for each read
        
        //iterate to first space in read name
        for (j = 0; rd_ptr[i][j] != ' ' && rd_ptr[i][j]; j++) {};
        
        //check that loop ended on a space
        if (rd_ptr[i][j] != ' ') {
            printf("append_att: error - unexpected read id format. aborting...\n");
            abort();
        }
        
        rd_ptr[i][j] = '\0';			//place null to split read string
        prefix_ptr = &rd_ptr[i][0];		//set prefix pointer to start of read
        suffix_ptr = &rd_ptr[i][j+1];	//set suffix pointer to index after the string split point
        
        //check that appending att array to read id will not exceed array bounds
        if (strlen(prefix_ptr) + strlen(att) + strlen(suffix_ptr) + 3 <= MAX_LINE) {
            sprintf(new_rd_name[i], "%s_%s %s",prefix_ptr, att, suffix_ptr); //construct new read name
            strcpy(rd_ptr[i], new_rd_name[i]); //replace old read name with new read name
        } else {
            printf("appnd_att: error - fastq line 1 is too long to add attributes. aborting...\n");
            abort();
        }
    }
    
    return 1;
}

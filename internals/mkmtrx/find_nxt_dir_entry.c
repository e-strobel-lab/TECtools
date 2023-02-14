//
//  find_nxt_dir_entry.c
//  
//
//  Created by Eric Strobel on 10/11/22.
//

#include "find_nxt_dir_entry.h"

#include <stdio.h>
#include <dirent.h>
#include <string.h>

/* find_nxt_dir_entry: search directory stream for entry with user-specified suffix. record
 the name and prefix of the entry in the entry_name and entry_prfx variables, respectively.
 returns 1 if next directory entry is found, 0 if next directory entry is not found */
int find_nxt_dir_entry(DIR * dirp, char * suffix, char * entry_name, char * entry_prfx)
{
    struct dirent * d_entry = NULL;   //pointer to current directory entry
    int entry_len = 0;                //length of current directory entry name
    int suffix_len = strlen(suffix);  //length of suffix string
    char * p_suffix = NULL;           //pointer to location of suffix in directory entry name
    int fnd_entry = 0;                //flag that user-specified directory entry was found
    
    while (!fnd_entry) { //until the next directory entry is found

        if ((d_entry = readdir(dirp)) != NULL) { //read next directory entry
            entry_len = strlen(d_entry->d_name); //record length of directory entry name
            
            if (entry_len > strlen(suffix)) {         //if the directory entry name is longer than the expected suffix
                strcpy(entry_name, d_entry->d_name);  //copy the directory entry name to the entry_name..
                strcpy(entry_prfx, d_entry->d_name);  //and entry_prfx arrays.
                
                if ((p_suffix = strstr(entry_prfx, suffix)) != NULL) { //if the directory entry contains the string suffix
                    if (!p_suffix[suffix_len]) { //match to suffix is at the end of directory entry name
                        *p_suffix = '\0';        //terminate entry_prfx at suffix to recover prefix
                        fnd_entry = 1;           //set fnd_entry to true
                    }
                }
            }
        } else {      //reached the end of the directory without finding match
            return 0; //return false
        }
    }
    return 1; //found next directory entry, return true
}

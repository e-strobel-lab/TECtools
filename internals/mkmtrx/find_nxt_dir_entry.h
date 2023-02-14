//
//  find_nxt_dir_entry.h
//  
//
//  Created by Eric Strobel on 10/11/22.
//

#ifndef find_nxt_dir_entry_h
#define find_nxt_dir_entry_h

#include <stdio.h>
#include <dirent.h>
#include <string.h>

/* find_nxt_dir_entry: search directory stream for entry with user-specified suffix. record
 the name and prefix of the entry in the entry_name and entry_prfx variables, respectively.
 returns 1 if next directory entry is found, 0 if next directory entry is not found */

int find_nxt_dir_entry(DIR * dirp, char * suffix, char * entry_name, char * entry_prfx);

#endif /* find_nxt_dir_entry_h */

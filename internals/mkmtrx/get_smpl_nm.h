//
//  get_smpl_nm.h
//  
//
//  Created by Eric Strobel on 10/11/22.
//

#ifndef get_smpl_nm_h
#define get_smpl_nm_h

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"

/* obtain sample name from prnt_dir input. expected parent directory is:
 <sample_name>/<analysis_name>, which is converted to:
 <sample_name>_analysis_name> */
int get_smpl_nm(char * prnt_dir, char * smpl_nm);

#endif /* get_smpl_nm_h */

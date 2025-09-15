//
//  read_VL_analysis_directories.c
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

#include "../../global/global_defs.h"
#include "../../mkmtrx/mkmtrx_defs.h"
#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"

#include "parse_VL_sample_name.h"

#include "read_VL_analysis_directories.h"

/* read_prnt_directory: read parent shapemapper 2 analysis directory and identify transcript length analysis sub-folders*/
int read_prnt_directory(SM2_analysis_directory * an_dir, int dir_num, sample_names * sn)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    struct dirent *dir = NULL; //directory entry pointer
    
    char an_sffx[10] = {"_analysis"};     //analysis directory suffix
    char tmp_d_name[MAX_NAME+1] = {'\0'}; //temporary directory name
    char tl_dir_nm[MAX_NAME+1] = {'\0'};  //relative path to transcript length directory
    
    int crnt_tl = 0;         //current transcript length
    int profiles_opened = 0; //number of profiles opened
    
    int ret = 0; //variable for storing snprintf return values
    
    DIR * prnt_dir = NULL; //pointer for parent directory
    int subdir_cnt = 0;    //number of subdirectories in parent directory

    int tmp_tl = 0; //temporary transcript length value
    int min_tl = 0; //minimum transcript length value
    
    //open parent directory
    if ((prnt_dir = opendir(an_dir->prnt_dir_nm)) == NULL) {
        printf("read_prnt_directory: error - failed to open parent directory. aborting...\n");
        abort();
    }
    
    //count number of subdirectories that contain substring "_analysis" in parent directory
    //this number is used for memory allocation below
    for (i = 0; (dir = readdir(prnt_dir)) != NULL; i++) {
        if (strstr(dir->d_name, an_sffx)) {  //test for "_analysis" suffix
            subdir_cnt++;                    //if found, increment subdirectory count
            strcpy(tmp_d_name, dir->d_name); //make copy of subdirectory name
            
            //get numerical id from subdirectory name
            if (((tmp_tl = get_nid(tmp_d_name, an_sffx)) < min_tl) || subdir_cnt == 1) {
                min_tl = tmp_tl;
            }
        }
    }
    
    //close parent directory
    if (closedir(prnt_dir) == EOF) {
        printf("read_prnt_directory: error - failed to close parent directory. aborting...\n");
        abort();
    }
    
    //allocate memory for SM2 profile location strings
    if ((an_dir->loc = calloc(subdir_cnt + min_tl, sizeof(*an_dir->loc))) == NULL) {
        printf("read_prnt_directory: error - failed to allocate memory for profile locations. aborting...\n");
        abort();
    }
    
    //allocate memory for SM2 profile data
    if ((an_dir->data = calloc(subdir_cnt + min_tl, sizeof(*an_dir->data))) == NULL) {
        printf("read_prnt_directory: error - failed to allocate memory for profile data. aborting...\n");
        abort();
    }
    
    //open parent shapemapper 2 analysis directory
    if ((prnt_dir = opendir(an_dir->prnt_dir_nm)) == NULL) {
        printf("read_prnt_directory: error - failed to open parent directory. aborting...\n");
        abort();
    }
    
    //read directory contents
    for (i = 0; (dir = readdir(prnt_dir)) != NULL; i++) { //until all directory contents are read
        
        if (isdigit(dir->d_name[0]) &&            //check whether the current directory entry name
            isdigit(dir->d_name[1]) &&            //begins with three digits followed by '_analysis'
            isdigit(dir->d_name[2]) &&
            !strcmp(&dir->d_name[3], an_sffx)) {
            
            //construct relative transcript length directory path
            ret = snprintf(tl_dir_nm, MAX_NAME, "%s/%s", an_dir->prnt_dir_nm, dir->d_name);
            if (ret >= MAX_NAME || ret < 0) {
                printf("read_prnt_directory: error - error when constructing transcript length directory name. aborting...\n");
                abort();
            }
            
            //assess the transcript length of the current analysis directory
            strcpy(tmp_d_name, dir->d_name); //copy the directory entry name to tmp_d_name
            crnt_tl = get_nid(tmp_d_name, an_sffx);
            
            if (crnt_tl < an_dir->min_tl) {  //if the current transcript length is less than min_tl
                an_dir->min_tl = crnt_tl;    //set min_tl to the current transcript length
            }
            
            if (crnt_tl > an_dir->max_tl) {  //if the current transcrip length is greater than max_tl
                an_dir->max_tl = crnt_tl;    //set max_tl to the current transcript length
            }
            
            //read the current transcript length directory
            profiles_opened += read_tl_directory(an_dir, dir_num, crnt_tl, tl_dir_nm, sn);
        }
    }
    
    //check that output directories identified during processing match the subdirectories identified during counting
    if (an_dir->outs_cnt != subdir_cnt) {
        printf("read_prnt_directory: error - output directory count does not match expected count. aborting...\n");
        abort();
    }
    
    //check that the minimum TL identified during processing matches the min TL identified during counting
    if (an_dir->min_tl != min_tl) {
        printf("read_prnt_directory: error - minimum transcript length does not match expected value (%d vs %d). aborting...\n", an_dir->min_tl, min_tl);
        abort();
    }
    
    //close parent directory
    if (closedir(prnt_dir) == EOF) {
        printf("read_prnt_directory: error - failed to close parent directory. aborting...\n");
        abort();
    }

    return profiles_opened; //return the number of shapemapper 2 reactivity profiles opened
}

/* read_tl_directory: read transcript length analysis directory and identify shapemapper 2 output sub-folders */
int read_tl_directory(SM2_analysis_directory * an_dir, int dir_num, int crnt_tl, char * tl_dir_nm, sample_names * sn)
{
    int i = 0; //general purpose index
    
    struct dirent *dir = NULL;   //directory entry pointer
    
    char out_sffx[5] = {"_out"}; //shapemapper 2 output directory suffix
    char * p_out_sffx = NULL;    //pointer to shapemapper 2 output directory suffix in directory entry name
    
    char out_dir_nm[MAX_NAME+1] = {'\0'}; //relative path to SM2 output directory
    int out_dir_cnt = 0;                  //number of output direcotries
    
    int ret = 0; //variable for storing snprintf return value
    
    DIR * crnt_tl_dir = NULL;  //pointer to transcript length directory
        
    //open the current transcript length directory
    if ((crnt_tl_dir = opendir(tl_dir_nm)) == NULL) {
        printf("read_tl_directory: error - failed to open transcript length %d directory. try increasing the file descriptor limit (run: ulimit -n <file descriptor limit> to increase the limit). aborting...\n", crnt_tl);
        abort();
    }
    
    //read directory contents
    for (i = 0; (dir = readdir(crnt_tl_dir)) != NULL; i++) {
        
        p_out_sffx = strstr(dir->d_name, out_sffx); //search for the shapemapper 2 output directory suffix
        
        if (p_out_sffx != NULL && !p_out_sffx[strlen(out_sffx)]) { //found shapemapper 2 output direcotry
            
            //if the sample name for this input has not yet been recorded,
            //record the name of the current shapemapper 2 output directory
            //in the sample names structure
            if (dir_num == sn->cnt) {
                ret = snprintf(sn->ipt[sn->cnt], MAX_NAME, "%s", dir->d_name);
                if (ret >= MAX_NAME || ret < 0) {
                    printf("read_out_directory: error - error when constructing reactivity profile name for transcript length %d. aborting...\n", crnt_tl);
                    abort();
                }
                
                remove_out_suffix(sn->ipt[sn->cnt]); //remove output suffix
                sn->cnt++;                           //increment sample name count
            }
            
            //construct relative shapemapper 2 output directory path
            ret = snprintf(out_dir_nm, MAX_NAME, "%s/%s", tl_dir_nm, dir->d_name);
            if (ret >= MAX_NAME || ret < 0) {
                printf("read_tl_directory: error - error when constructing SM2 output directory name for transcript length %d. aborting...\n", crnt_tl);
                abort();
            }
            
            out_dir_cnt++; //increment output directory count
        }
    }
    
    if (out_dir_cnt == 1) { //if a single SM2 output directory was found
        return (read_SM2out_directory(an_dir, crnt_tl, out_dir_nm)); //read the SM2 output directory
    } else { //if more than one SM2 output directory was found, abort
        printf("read_tl_directory: error - more than one SM2 output directory found in transcript length %d analysis directory. aborting...\n", crnt_tl);
        abort();
    }
}

/* read_SM2out_directory: read shapemapper 2 output directory and open reactivity profile file */
int read_SM2out_directory(SM2_analysis_directory * an_dir, int crnt_tl, char * out_dir_nm)
{
    int i = 0; //general purpose index
    
    struct dirent *dir = NULL; //directory entry pointer
    
    char profile_sffx[13] = {"_profile.txt"}; //reactivity profile file suffix
    char * p_profile_sffx = NULL;             //pointer to reactivity profile file suffix
    
    char profile_nm[MAX_NAME+1] = {'\0'};     //relative path to reactivity profile file
    int profile_cnt = 0;                      //number of reactivity profile files
    
    int ret = 0; //variable for storing snprintf return value
    
    DIR * crnt_out_dir = NULL;  //pointer to output directory
    
    //open the shapemapper 2 output directory
    if ((crnt_out_dir = opendir(out_dir_nm)) == NULL) {
        printf("read_SM2out_directory: error - failed to open SM2 output directory for transcript length %d. aborting...\n", crnt_tl);
        abort();
    } else {
        an_dir->outs_cnt++; //increment out directory opened counter
    }
    
    //read directory contents
    for (i = 0; (dir = readdir(crnt_out_dir)) != NULL; i++) {
        
        p_profile_sffx = strstr(dir->d_name, profile_sffx); //search for the reactivity profile suffix
        
        if (p_profile_sffx != NULL && !p_profile_sffx[strlen(profile_sffx)]) { //found reactivity profile
            
            //generate relative path to reactivity profile file
            ret = snprintf(profile_nm, MAX_NAME, "%s/%s", out_dir_nm, dir->d_name);
            if (ret >= MAX_NAME || ret < 0) {
                printf("read_SM2out_directory: error - error when constructing reactivity profile name for transcript length %d. aborting...\n", crnt_tl);
                abort();
            }
            profile_cnt++; //increment reactivity profile file count
        }
    }
    
    if (profile_cnt == 1) { //if a single reactivity profile was found
        
        //allocate memory for profile relative file path
        if (((an_dir->loc[crnt_tl]) = malloc((strlen(profile_nm)+1) * sizeof(*(an_dir->loc[crnt_tl])))) == NULL) {
            printf("read_SM2out_directory: error - memory allocation for profile relative filepath storage failed. aborting...\n");
            abort();
        }
        strcpy(an_dir->loc[crnt_tl], profile_nm);
        return 1;
        
    } else if (!profile_cnt) { //no reactivity profiles were found
        return 0;
        
    } else if (profile_cnt > 1) { //if more than one reactivity profile was found, abort
        printf("read_SM2out_directory: error - more than one reactivity profile found in transcript length %d SM2 output directory. aborting...\n", crnt_tl);
        abort();
        
    } else { //negative profile count?
        printf("read_SM2out_directory: error - negative reactivity profile in transcript length %d SM2 output directory. how??? aborting...\n", crnt_tl);
        abort();
    }
}

/* set_opnd_profile_bounds: set min and max opened profiles in SM2 analysis directory structure */
void set_opnd_profile_bounds(SM2_analysis_directory * an_dir)
{
    int i = 0;             //general purpose index
    int min_prf_set = 0;   //flag that minimum profile was set
    int max_prf = 0;  //maximum opened profile
    
    for (i = an_dir->min_tl; i <= an_dir->max_tl; i++) {
        if (an_dir->loc[i] != NULL) { //if the current transcript length profile was opened
            if (!min_prf_set) {       //if the minimum profile was not yet set
                an_dir->min_prf = i;  //set the minimum profile to the current transcript length
                min_prf_set = 1;      //turn on flag that minimum profile was set
            }
            an_dir->max_prf = i;      //set the max profile to the current transcript length
        }
    }
    
    return;
}

/* get_nid: get numerical id of variant */
int get_nid(char * str, char * suffix)
{
    int i = 0;               //general purpose index
    char * p_suffix = NULL;  //pointer to suffix that will be truncated
    
    //find suffix in string
    if ((p_suffix = strstr(str, suffix)) == NULL) {
        printf("get_nid: error - suffix '%s' was not found in string '%s'. aborting...\n", suffix, str);
        abort();
    }
    p_suffix[0] = '\0'; //truncate suffix
    
    //check that numerical id is composed of digits
    for (i = 0; str[i]; i++) {
        if (!isdigit(str[i])) {
            printf("get_nid: error - non-digit character in numerical id string '%s'. aborting...\n", str);
            abort();
        }
    }
    
    return atoi(str); //return numerical id
}

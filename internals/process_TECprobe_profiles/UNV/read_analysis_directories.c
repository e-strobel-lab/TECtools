//
//  read_analysis_directories.c
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <dirent.h>

#include "../../global/global_defs.h"
#include "../../mkmtrx/mkmtrx_defs.h"
#include "../process_TECprobe_profiles_defs.h"
#include "../process_TECprobe_profiles_structs.h"

#include "parse_sample_name.h"

#include "read_analysis_directories.h"

const char empty_SM2out[6] = "empty";

/* read_prnt_directory: read parent shapemapper 2 analysis directory and identify target analysis sub-folders*/
int read_prnt_directory(SM2_analysis_directory * an_dir, int dir_num, sample_names * sn)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    struct dirent *dir = NULL; //directory entry pointer
    
    char an_sffx[10] = {"_analysis"};     //analysis directory suffix
    char tmp_d_name[MAX_NAME+1] = {'\0'}; //temporary directory name
    char trg_dir_nm[MAX_NAME+1] = {'\0'};  //relative path to target directory
    
    int crnt_id = 0;         //current target id
    int profiles_opened = 0; //number of profiles opened
    
    int ret = 0; //variable for storing snprintf return values
    
    DIR * prnt_dir = NULL; //pointer for parent directory
    int subdir_cnt = 0;    //number of subdirectories in parent directory

    int tmp_id = 0; //temporary target id value
    int min_id = 0; //minimum target id value
    int max_id = 0; //maximum target id value
    
    //open parent directory
    if ((prnt_dir = opendir(an_dir->prnt_dir_nm)) == NULL) {
        printf("read_prnt_directory: error - failed to open parent directory. aborting...\n");
        abort();
    }
    
    //count number of target analysis  subdirectories in parent directory
    //this number is used for memory allocation below
    for (i = 0; (dir = readdir(prnt_dir)) != NULL; i++) {
        if (test_trg_analysis_dir_format(dir->d_name, an_sffx)) { //test trg analysis dir format
            subdir_cnt++;                          //if found, increment subdirectory count
            strcpy(tmp_d_name, dir->d_name);       //make copy of subdirectory name
            tmp_id = get_nid(tmp_d_name, an_sffx); //get numerical id from subdirectory name
            
            if (tmp_id < min_id || subdir_cnt == 1) { //if assessing 1st id, or tmp_id val is < min_id
                min_id = tmp_id;                      //set min_id to tmp_id
            }
            
            if (tmp_id > max_id || subdir_cnt == 1) { //if assessing 1st id, or tmp_id val is > max_id
                max_id = tmp_id;                      //set max_id to tmp_id
            }
        }
    }
    
    an_dir->sd_cnt = subdir_cnt; //set target analysis subdirectory count

    //close parent directory
    if (closedir(prnt_dir) == EOF) {
        printf("read_prnt_directory: error - failed to close parent directory. aborting...\n");
        abort();
    }
    
    //allocate memory for SM2 profile location strings
    if ((an_dir->loc = calloc(max_id+1, sizeof(*an_dir->loc))) == NULL) {
        printf("read_prnt_directory: error - failed to allocate memory for profile locations. aborting...\n");
        abort();
    }
    
    //allocate memory for SM2 profile data
    if ((an_dir->data = calloc(max_id+1, sizeof(*an_dir->data))) == NULL) {
        printf("read_prnt_directory: error - failed to allocate memory for profile data. aborting...\n");
        abort();
    }
    
    //allocate memory for lookup table
    if ((an_dir->indx = calloc(subdir_cnt+1, sizeof(*an_dir->indx))) == NULL) {
        printf("read_prnt_directory: error - failed to allocate memory for index table. aborting...\n");
        abort();
    }
    an_dir->indx[subdir_cnt] = INT_MAX; //last indx val is a sentinel to signal end of array
    
    //open parent shapemapper 2 analysis directory
    if ((prnt_dir = opendir(an_dir->prnt_dir_nm)) == NULL) {
        printf("read_prnt_directory: error - failed to open parent directory. aborting...\n");
        abort();
    }
    
    //read directory contents
    for (i = 0; (dir = readdir(prnt_dir)) != NULL; i++) { //until all directory contents are read
        
        if (test_trg_analysis_dir_format(dir->d_name, an_sffx)) { //check trg analysis dir format
            
            //construct relative target directory path
            ret = snprintf(trg_dir_nm, MAX_NAME, "%s/%s", an_dir->prnt_dir_nm, dir->d_name);
            if (ret >= MAX_NAME || ret < 0) {
                printf("read_prnt_directory: error - error when constructing target directory name. aborting...\n");
                abort();
            }
            
            //assess the target id of the current analysis directory
            strcpy(tmp_d_name, dir->d_name);        //copy the directory entry name to tmp_d_name
            crnt_id = get_nid(tmp_d_name, an_sffx); //set current target id
            
            //note: min_id and max_id are initialized to INT_MAX and INT_MIN respectively
            if (crnt_id < an_dir->min_id) {  //if the current target id is less than min_id
                an_dir->min_id = crnt_id;    //set min_id to the current target id
            }
            
            if (crnt_id > an_dir->max_id) {  //if the current target id is greater than max_id
                an_dir->max_id = crnt_id;    //set max_id to the current target id
            }
            
            //read the current target directory
            profiles_opened += read_target_directory(an_dir, dir_num, crnt_id, trg_dir_nm, sn);
        }
    }
    
    //check that output directories identified during processing match the subdirectories identified during counting
    if (an_dir->outs_cnt != subdir_cnt) {
        printf("read_prnt_directory: error - output directory count does not match expected count. aborting...\n");
        abort();
    }
    
    //check that the minimum id identified during processing matches the min id identified during counting
    if (an_dir->min_id != min_id) {
        printf("read_prnt_directory: error - minimum id does not match expected value (%d vs %d). aborting...\n", an_dir->min_id, min_id);
        abort();
    }
    
    //check that the maximum id identified during processing matches the max id identified during counting
    if (an_dir->max_id != max_id) {
        printf("read_prnt_directory: error - maximum id does not match expected value (%d vs %d). aborting...\n", an_dir->max_id, max_id);
        abort();
    }
    
    //set target indices
    for (i = 0, j = 0; i <= max_id; i++) { //for every possible target id
        if (an_dir->loc[i] != NULL) {      //if the location string was set for the current target
            an_dir->indx[j++] = i;         //set the target index
        }
    }
    
    //close parent directory
    if (closedir(prnt_dir) == EOF) {
        printf("read_prnt_directory: error - failed to close parent directory. aborting...\n");
        abort();
    }

    return profiles_opened; //return the number of shapemapper 2 reactivity profiles opened
}

/* read_target_directory: read target analysis directory and identify shapemapper 2 output sub-folders */
int read_target_directory(SM2_analysis_directory * an_dir, int dir_num, int crnt_id, char * trg_dir_nm, sample_names * sn)
{
    int i = 0; //general purpose index
    
    struct dirent *dir = NULL;   //directory entry pointer
    
    char out_sffx[5] = {"_out"}; //shapemapper 2 output directory suffix
    char * p_out_sffx = NULL;    //pointer to shapemapper 2 output directory suffix in directory entry name
    
    char out_dir_nm[MAX_NAME+1] = {'\0'}; //relative path to SM2 output directory
    int out_dir_cnt = 0;                  //number of output direcotries
    
    int ret = 0; //variable for storing snprintf return value
    
    DIR * crnt_trg_dir = NULL;  //pointer to current target directory
        
    //open the current target directory
    if ((crnt_trg_dir = opendir(trg_dir_nm)) == NULL) {
        printf("read_target_directory: error - failed to open target %d directory. try increasing the file descriptor limit (run: ulimit -n <file descriptor limit> to increase the limit). aborting...\n", crnt_id);
        abort();
    }
    
    //read directory contents
    for (i = 0; (dir = readdir(crnt_trg_dir)) != NULL; i++) {
        
        p_out_sffx = strstr(dir->d_name, out_sffx); //search for the shapemapper 2 output directory suffix
        
        if (p_out_sffx != NULL && !p_out_sffx[strlen(out_sffx)]) { //found shapemapper 2 output direcotry
            
            //if the sample name for this input has not yet been recorded,
            //record the name of the current shapemapper 2 output directory
            //in the sample names structure
            if (dir_num == sn->cnt) {
                ret = snprintf(sn->ipt[sn->cnt], MAX_NAME, "%s", dir->d_name);
                if (ret >= MAX_NAME || ret < 0) {
                    printf("read_out_directory: error - error when constructing reactivity profile name for target %d. aborting...\n", crnt_id);
                    abort();
                }
                
                remove_id_and_suffix(sn->ipt[sn->cnt], out_sffx); //remove output suffix
                sn->cnt++;                                        //increment sample name count
            }
            
            //construct relative shapemapper 2 output directory path
            ret = snprintf(out_dir_nm, MAX_NAME, "%s/%s", trg_dir_nm, dir->d_name);
            if (ret >= MAX_NAME || ret < 0) {
                printf("read_target_directory: error - error when constructing SM2 output directory name for target %d. aborting...\n", crnt_id);
                abort();
            }
            
            out_dir_cnt++; //increment output directory count
        }
    }
    
    if (out_dir_cnt == 1) { //if a single SM2 output directory was found
        return (read_SM2out_directory(an_dir, crnt_id, out_dir_nm)); //read the SM2 output directory
    } else { //if more than one SM2 output directory was found, abort
        printf("read_target_directory: error - more than one SM2 output directory found in target %d analysis directory. aborting...\n", crnt_id);
        abort();
    }
}

/* read_SM2out_directory: read shapemapper 2 output directory and open reactivity profile file */
int read_SM2out_directory(SM2_analysis_directory * an_dir, int crnt_id, char * out_dir_nm)
{
    extern const char empty_SM2out[6];
    
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
        printf("read_SM2out_directory: error - failed to open SM2 output directory for target %d. aborting...\n", crnt_id);
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
                printf("read_SM2out_directory: error - error when constructing reactivity profile name for target %d. aborting...\n", crnt_id);
                abort();
            }
            profile_cnt++; //increment reactivity profile file count
        }
    }
    
    if (profile_cnt == 1) { //if a single reactivity profile was found
        
        //allocate memory for profile relative file path
        if ((an_dir->loc[crnt_id] = malloc((strlen(profile_nm)+1) * sizeof(*(an_dir->loc[crnt_id])))) == NULL) {
            printf("read_SM2out_directory: error - memory allocation for profile relative filepath storage failed. aborting...\n");
            abort();
        }
        strcpy(an_dir->loc[crnt_id], profile_nm);
        return 1;
        
    } else if (!profile_cnt) { //no reactivity profiles were found
        
        //store location as "empty"
        if ((an_dir->loc[crnt_id] = malloc((strlen(empty_SM2out)+1) * sizeof(*(an_dir->loc[crnt_id])))) == NULL) {
            printf("read_SM2out_directory: error - memory allocation for profile relative filepath storage failed. aborting...\n");
            abort();
        }
        strcpy(an_dir->loc[crnt_id], empty_SM2out);
        return 0;
        
    } else if (profile_cnt > 1) { //if more than one reactivity profile was found, abort
        printf("read_SM2out_directory: error - more than one reactivity profile found in target %d SM2 output directory. aborting...\n", crnt_id);
        abort();
        
    } else { //negative profile count?
        printf("read_SM2out_directory: error - negative reactivity profile in target %d SM2 output directory. how??? aborting...\n", crnt_id);
        abort();
    }
}

/* set_len_range: set minimum and maximum transcript lengths */
void set_len_range(SM2_analysis_directory * an_dir)
{
    int i = 0; //general purpose index
    
    int * ix = &an_dir->indx[0]; //set pointer to target indices
    
    for (i = 0; ix[i] <= an_dir->max_id; i++) {                            //for every target
        if (an_dir->data[ix[i]].trg_nt_cnt < an_dir->len[MIN] || i == 0) { //if 1st target or crnt trg < crnt min
            an_dir->len[MIN] = an_dir->data[ix[i]].trg_nt_cnt;             //set crnt min to crnt trg
        }
        
        if (an_dir->data[ix[i]].trg_nt_cnt > an_dir->len[MAX] || i == 0) { //if 1st target or crnt trg > crnt max
            an_dir->len[MAX] = an_dir->data[ix[i]].trg_nt_cnt;             //set crnt max to crnt trg
        }
    }
}

/* test_trg_analysis_dir_format: test that target analysis directory
 name conforms to expected format. return length of id string. */
int test_trg_analysis_dir_format(char * str, char * suffix)
{
    int i = 0;               //general purpose index
    int id_len = 0;          //length of id string
    char * p_suffix = NULL;  //pointer to suffix
    
    //check that id string, which precedes first underscore, is digits
    for (i = 0; str[i] && str[i] != '_'; i++) {
        
        if (isdigit(str[i])) {       //if char is a digit
            id_len++;                //increment id length
            
        } else if (str[i] != '_') {  //if non-digit char is not an underscore
            return 0;                //return fail
        }
    }
    
    if (str[i] != '_') {  //if loop did not end on an underscore
        return 0;   //return fail
    }
    
    p_suffix = &str[i];               //set pointer to underscore
    if (!strcmp(p_suffix, suffix)) {  //check suffix match
        return 1;                     //if match, return pass
    } else {
        return 0;                     //if not match, return fail
    }
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

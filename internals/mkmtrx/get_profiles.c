//
//  get_profiles.c
//  
//
//  Created by Eric Strobel on 10/11/22.
//

#include "get_profiles.h"

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "./cotrans_mtrx.h"
#include "./mkmtrx_defs.h"
#include "./mkmtrx_structs.h"
#include "./find_nxt_dir_entry.h"
#include "./parse_log.h"

#include "../process_TECprobe_profiles/VL/parse_VL_sample_name.h"

/* get_profiles: read input directory to find shapemapper output files.
 then call parse_profile to store reactivity data in cotrans matrix struct*/
int get_profiles(char * prnt_dir, cotrans_matrix * mtrx, alignment_stats * algn, char * smpl_nm, int rct_typ, int preprocessed, int test_SM2_data)
{
    //printf("in get_profiles\n");
    
    DIR * dirp = NULL;         //pointer to parent directory, which contains all transcript length analysis directories
    DIR * sub_dirp = NULL;     //pointer for reading subdirectories
    DIR * crnt_dp[PRFL_SRCH_DEPTH] = {NULL}; //pointer to directory pointers. used when reading/validating directory hierarchy
    
    /* open parent directory. the parent directory remains
     open for the duration of the function call*/
    if ((dirp = opendir(prnt_dir)) == NULL) {
        printf("get_profiles: error - failed to open parent directory. aborting...\n");
        abort();
    }
    crnt_dp[0] = dirp; //set pointer to parent directory pointer
    
    int i = 0;  //general purpose index
    int j = 0;  //general purpose index

    int transcript_len = 0;          //transcript length for current directory
    //int min_transcript = MIN_INIT;   //minimum observed transcript length
    //int max_transcript = MAX_INIT;   //maximum observed transcript length
    
    char ntry_name[MAX_LINE] = {0};  //full name of next directory entry
    char ntry_prfx[MAX_LINE] = {0};  //prefix of next directory entry, currently used to get transcript length string
    
    char xpctd_suffix[PRFL_SRCH_DEPTH][MAX_LINE] = {"_analysis", "_out", "_profile.txt"}; //identifiers for expected subdirectories
    
    //(MAX_LINE*2)+2 is to silence warning. the string can potentially be longer than this (though it should
    //never be) but a check is done to make sure the string fits in the array
    char file_loc[(MAX_LINE*2)+2] = {0};   //array used to assemble the location of shapemapper output files
    char tmp_file_loc[MAX_LINE] = {0};
    int fnd_d_entry = 0;             //flag that directory entry was found
    
    int proceed = 1;                 //flag to proceed with the next iteration of outermost while loop
    
    int sample_name_set = 0;         //flag that sample name was set
    
    mtrx->sq[0] = '>';
    mtrx->tl[MIN] = MIN_INIT;
    mtrx->tl[MAX] = MAX_INIT;
    
    while (proceed) {
        strcpy(file_loc, prnt_dir); //initialize file_loc with the parent directory string
        
        /* the expected location of shapemapper output files is:
         
         <parent_directory>/<transcript_len>_analysis/<sample_info>_out/<sample_and_analysis_info>profile.txt
         
         the for loop below assembles the shapemapper output file location by sequentially searching
         each directory to confirm the existence of the next expected directory entry and, if found,
         appending the next directory entry to the file_loc string.
         
         */
        
        for (i = 0, fnd_d_entry = 0; i < PRFL_SRCH_DEPTH && proceed; i++) {
            
            /* search current directory for next expected directory
             entry as defined by the xpcted_suffix variable.*/
            fnd_d_entry = find_nxt_dir_entry(crnt_dp[i], xpctd_suffix[i], ntry_name, ntry_prfx);
            
            if (!i && !fnd_d_entry) {               //assessed all parent director entries
                proceed = 0;                        //set proceed to false to exit all loops
                
            } else if (i && !fnd_d_entry) {         //could not find expected subdirectory entry, throw error and abort
                printf("get_profiles: warning - could not find directory entry with suffix %s in directory %s. something may be missing.\n\n", xpctd_suffix[i], file_loc);
                
            } else if (!i && fnd_d_entry) {         //found shapemapper analysis directory
                
                //check that shapemapper analysis directory prefix is composed of digits
                for (j = 0; ntry_prfx[j]; j++) {
                    if (!isdigit(ntry_prfx[j])) {
                        printf("get_profiles: error - unexpected non digit character in prefix of shapemapper analysis directory. aborting...\n");
                        abort();
                    }
                    
                }
                transcript_len = atoi(ntry_prfx);   //set transcript length using ntry prfx
                
                //validate transcript length and compare to min and max transcript values
                if (transcript_len >= MIN_INIT || transcript_len <= MAX_INIT) { //check that transcript length is in bounds
                    printf("get_profiles: error - transcript %d is out of bounds (1-999 nts). aborting...", transcript_len);
                    abort();
                    
                } else { //transcript length is in bounds
                    
                    if (transcript_len < mtrx->tl[MIN]) { //compare transcript_len to current min_transcript
                        mtrx->tl[MIN] = transcript_len;   //if shorter, set new min transcript
                    }
                    
                    if (transcript_len > mtrx->tl[MAX]) { //compare transcript_len to current max_transcript
                        mtrx->tl[MAX] = transcript_len;   //if longer, set new max transcript
                    }
                }
                
            } else if (i && fnd_d_entry) {  //found next subdirectory
                if (closedir(sub_dirp)) {   //close previous subdirectory
                    printf("get_profiles: error - failed to close directory %s. aborting...\n", file_loc);
                    abort();
                }
            }

            
            if (proceed) { //found next directory entry
                
                //set sample name
                if (i == 1 && !sample_name_set) { //if reading SM2 output dir and sample name is not set
                    strcpy(smpl_nm, ntry_name);            //set the sample name
                    remove_id_and_suffix(smpl_nm, "_out"); //remove the id and out suffix from the dir name
                    sample_name_set = 1;                   //set flag that the sample name was set
                }
                
                strcpy(tmp_file_loc, file_loc); //copy current file_loc to tmp_file_loc
                
                //test if new file location fits in array, if so, assemble new file_loc
                if ((strlen(tmp_file_loc) + strlen(ntry_name)) >= MAX_LINE) {
                    printf("get_profiles: error - file location too long. aborting...\n");
                    abort();
                } else {
                    sprintf(file_loc, "%s/%s", tmp_file_loc, ntry_name); //append next directory entry to tmp_file_loc string to assemble new file_loc
                }
                
                //reading transcript analysis directory, parse the log file for alignment details
                //only perform this step if not analyzing test SM2 data or preprocessed data
                if (!i && !test_SM2_data && !preprocessed) {
                    parse_log(file_loc, &algn[transcript_len]);
                }
                
                //when the directory entry search is expected to recover a subdirectory, open the
                //subdirectory. the final directory entry search recovers the shapemapper output file,
                //so it a new directory is not opened on the final loop iteration
                if (i < PRFL_SUBDIRS) {
                    if ((sub_dirp = opendir(file_loc)) == NULL) { //open next directory
                        printf("get_profiles: error - failed to open parent directory. aborting...\n");
                        abort();
                    }
                    //set the crnt_dp pointer for the next loop iteration
                    //to point to the newly opened subdirectory
                    crnt_dp[i+1] = sub_dirp;
                }
            }
        }
        if (proceed && fnd_d_entry) { //found shapemapper output file
            parse_profile(file_loc, transcript_len, mtrx, rct_typ); //store reactivity data in cotrans_matrix struct
        }
    }
    
    if (closedir(dirp)) { //close parent directory
        printf("get_profiles: error - failed to close directory %s. aborting...\n", prnt_dir);
        abort();
    }
    
    //subtract length of sequence string starting at index 1 (index 0 = '>') and add 1
    //this allows accurate nucleotide numbering if a leading segment of the transcript
    //is omitted
    mtrx->nt[MIN] = mtrx->tl[MAX]-strlen(&mtrx->sq[1])+1;
    mtrx->nt[MAX] = mtrx->tl[MAX];
    
    printf("%s\n", mtrx->sq);
    printf("nucleotides:\tmin=%3d, max=%3d\n", mtrx->nt[MIN], mtrx->nt[MAX]);
    printf("transcripts:\tmin=%3d, max=%3d\n", mtrx->tl[MIN], mtrx->tl[MAX]);
    return 1;
}

/* parse_profile: shapemapper output file to validate file integrity
 and copy reactivity values into cotrans matrix struct */
int parse_profile(char * prfl_loc, int len, cotrans_matrix * mtrx, int rct_typ)
{
    
    FILE * prfl_p = NULL;         //pointer to profile file
    
    int i = 0;                    //general purpose index
    int j = 0;                    //general purpose index
    char line[MAX_LINE] = {0};    //array to store lines from shapemapper output file
    int include_row = 0;          //flag to include current row in cotrans matrix
    int crnt_col = 0;             //number of the current column being evaluated
    
    char tmp_val[MAX_LINE] = {0}; //array to store values from shapemapper output file
    char * val2cpy = NULL;        //pointer to the value that will be copied into cotrans_matrix, used to mask nan as 0.0
    char zero[4] = {"0.0"};       //mask for nan values
    int val_cnt[3] = {0};         //number of values copied to the cotrans matrix from the current shapemapper output file
    val_cnt[0] = val_cnt[1] = val_cnt[2] = 1;
    int * p_val_cnt = NULL;       //pointer to current matrix value counter
    
    int store_column_val = 0;     //flag that current shapemapper output column has a value that will be stored
    char ** array2use = NULL;     //pointer to matrix location that will be used to store value
    
    int sq_index = 1;             //index used when storing sequence
    
    //open shapemapper output file
    if ((prfl_p = fopen(prfl_loc, "r")) == NULL) {
        printf("parse_profile: warning - failed to open profile file in location %s. aborting\n", prfl_loc);
        abort();
    }
    
    
    for (i = 0; fgets(line, MAX_LINE, prfl_p) != NULL; i++) { //until all lines of the shapemapper output file have been read
        if (!i) { //if reading the first line of the shapemapper output file (which contains column headers)
            if (vldt_SMO_hdr(line) == -1) { //validate that column headers that will be used are correct
                printf("parse_profile: error - unexpected format for shapemapper output file %s. aborting...\n", prfl_loc);
                abort();
            }
            
        } else { //if reading data lines from the shapemapper ouput file
            
            //sequentially read the values from each column of the shapemapper output file
            
            //determine if the current row should be included in the reactivity matrix by
            //checking whether the the sequence is uppercase or lowercase. lowercase rows
            //correspond to the SC1 structure cassette and are excluded from the matrix.
            //uppercase rows correspond to the target RNA and are included in the matrix.
            
            //if the current row is to be included, copy the reactivity_profile value to
            //the reactivity matrix
            
            for (include_row = 1, crnt_col = 0, i = 0; line[i] && include_row; crnt_col++, i++) {
                
                store_column_val = 0;
                array2use = NULL;
                p_val_cnt = NULL;
                
                //copy value from the current column to tmp_val array
                for (j = 0; line[i] && !isspace(line[i]); i++, j++) {
                    tmp_val[j] = line[i];
                }
                tmp_val[j] = '\0';
                
                
                if (crnt_col == SEQUENCE && islower(tmp_val[0])) { //if reading the sequence column, and the char is lower case...
                    include_row = 0;                               //do not include the current row (exit loop, proceed to next row)
                    
                } else if (crnt_col == SEQUENCE && len == mtrx->tl[MAX]) {
                    mtrx->sq[sq_index++] = tmp_val[0];
                    
                } else if (crnt_col == MOD_EFFECTIVE) {            //if reading the modified_effective_depth column...
                    store_column_val = 1;                          //set flag to store column value
                    array2use = &mtrx->mod[len][val_cnt[MOD]];     //set pointer to corresponding mod array location
                    p_val_cnt = &val_cnt[MOD];                     //set pointer to mod value counter
                    
                } else if (crnt_col == UNT_EFFECTIVE) {            //if reading the untreated_effective_depth column...
                    store_column_val = 1;                          //set flag to store column value
                    p_val_cnt = &val_cnt[UNT];                     //set pointer to corresponding unt array location
                    array2use = &mtrx->unt[len][val_cnt[UNT]];     //set pointer to unt value counter
                    
                } else if (crnt_col == rct_typ) {                  //if reading the output reactivity column...
                    store_column_val = 1;                          //set flag to store column value
                    array2use = &mtrx->vals[len][val_cnt[RCT]];    //set pointer to corresponding reactivity arrry loc
                    p_val_cnt = &val_cnt[RCT];                     //set pointer to rct vlaue counter
                }
                
                if (store_column_val) {
                    //allocate memory for storing the value (as a string of chars)
                    if (((*array2use) = malloc((strlen(tmp_val)+1) * sizeof(**array2use))) == NULL) {
                        printf("store_mtrx: error - memory allocation for matrix field value failed. aborting...\n");
                        abort();
                    }
                    
                    //set val2cpy to tmp_val if value is numerical, or 0.0 if value is nan
                    val2cpy = (strcmp(tmp_val, "nan")) ? &tmp_val[0] : &zero[0];
                    strcpy(*array2use, val2cpy); //copy value to matrix
                    (*p_val_cnt)++;              //track number of values recorded in matrix
                }
            }
        }
    }
    
    if (len == mtrx->tl[MAX]) {
        mtrx->sq[sq_index++] = '\0'; //terminate sequence string
    }
    
    
    
    //printf("%s\n", *array2use);
    
    //confirm number of values obtained from shapemapper output matches transcript lengths
    for (i = 0; i < 3; i++) {
        if (val_cnt[i]-1 != len) { //if number of values recorded in matrix does not match transcript length, throw error and abort
            printf("parse_profile: error - number of included ");
            switch (i) {
                case MOD: printf("modified effective depth"); break;
                case UNT: printf("untreated effective depth"); break;
                case RCT: printf("reactivity"); break;
                default:
                    break;
            }
            
            printf(" values in shapemapper output file does not match the transcript length. aborting...");
            abort();
        }
    }
    
    //allocate memory for storing the number of values in the current row
    char len_str[MAX_LINE] = {0};            //length string
    sprintf(len_str, "%d", val_cnt[RCT]-1);  //convert length from integer to string
    
    array2use = &mtrx->vals[len][0]; //set storage location to vals index zero of current row and allocate memory
    if (((*array2use) = malloc((strlen(len_str)+1) * sizeof(**array2use))) == NULL) {
        printf("store_mtrx: error - memory allocation for matrix field value failed. aborting...\n");
        abort();
    }
    strcpy(*array2use, len_str);    //store length string
    
    if (fclose(prfl_p) == EOF) { //close shapemapper output file
        fprintf(prfl_p, "parse_profile: error - failed to close shapemapper output file. Aborting program...\n");
        abort();
    }
    
    return 1;
}


/* vldt_SMO_hdr: validate shapemapper output file header by checking first column,
 columns that will be used in parse_profiles, and the number of columns */
int vldt_SMO_hdr(char * line)
{
    int i = 0;       //general purpose index
    int j = 0;       //general purpose index
    int col_num = 0; //number of columns

    char tmp_val[MAX_LINE] = {0}; //array to store column headers
    
    for (col_num = 0, i = 0; line[i]; col_num++, i++) { //iterate through column header line, counting columns
        
        for (j = 0; line[i] && !isspace(line[i]); i++, j++) { //get column header string
            tmp_val[j] = line[i];
        }
        tmp_val[j] = '\0';
                
        if (col_num == NUCLEOTIDE && strcmp(tmp_val, "Nucleotide")) {    //check Nucleotide column header
            return -1; //error
        } else if (col_num == SEQUENCE && strcmp(tmp_val, "Sequence")) { //check Sequence column header
            return -1; //error
        } else if (col_num == MOD_EFFECTIVE && strcmp(tmp_val, "Modified_effective_depth")) { //check Modified_effective_depth column header
            return -1; //error
        } else if (col_num == UNT_EFFECTIVE && strcmp(tmp_val, "Untreated_effective_depth")) { //check Untreated_effective_depth column header
            return -1; //error
        } else if (col_num == REACTIVITY_PROFILE && strcmp(tmp_val, "Reactivity_profile")) { //check Reactivity_profile column header
            return -1; //error
        } else if (col_num == HQ_PROFILE && strcmp(tmp_val, "HQ_profile")) { //check HQ_profile column header
            return -1; //error
        } else if (col_num == NORM_PROFILE && strcmp(tmp_val, "Norm_profile")) { //check Norm_profile column header
            return -1; //error
        }
    }
    
    if (col_num >= 27 && col_num <= 29) { //check that file contains minimum expected number of columns
        return col_num; //passed all tests, return column number
    } else {
        return -1; //error
    }
}


/* check_contiguity: check that reactivity matrix is
 contiguous from minimum to maximum transcript lengths */
int check_contiguity(cotrans_matrix * mtrx)
{
    int i = 0;                //general purpose index
    int missing_len[MAX_ROW]; //array to map what transcript lengths are missing
    int contiguous = 1;       //flag that matrix is contiguous
    
    for (i = mtrx->tl[MIN]; i <= mtrx->tl[MAX] && i <= MAX_ROW; i++) { //for each row of the reactivity matrix
        if (mtrx->vals[i][1] == NULL) { //if the vals array is NULL
            missing_len[i] = 1;         //record this missing length
            contiguous = 0;             //set the contiguous flag to false
        }
    }
    
    if (!contiguous) { //if the matrix is not contiguous, throw error and abort
        printf("get_profiles: error - matrix is not contiguous. missing lengths: ");
        for (i = mtrx->tl[MIN]; i <= mtrx->tl[MAX] && i <= MAX_ROW; i++) {
            if (missing_len[i]) {
                printf("%d ", i);
            }
        }
        printf("\n");
    }
    
    return 1;
}

//TODO: Add option for providing custom downstream constant sequence instead of HP4
int set_enriched_lengths(cotrans_matrix *mtrx, int incld_up2, int excld_trmnl)
{
    char hp4[34] = {"UGAUUCGUCAGGCGAUGUGUGCUGGAAGACAUU"};
    char *p_hp4 = NULL;
    
    int len = strlen(mtrx->sq);            //length of longest transcript in matrix
    int strt_intrnl_Nrchd = mtrx->tl[MIN]; //start of enriched internal transcripts
    int end_intrnl_Nrchd = 0;              //end of enriched internal transcripts
    int strt_trmnl_Nrchd = 0;              //start of terminally enriched transcripts
    int end_trmnl_Nrchd = 0;               //end of terminally enriched transcripts
    
    int i = 0; //general purpose index
    
    p_hp4 = strstr(mtrx->sq, hp4); //pointer to HP4 sequence
    if (p_hp4 != NULL && p_hp4[strlen(hp4)] == '\0') {
        
        if (incld_up2) {                       //if include-up-to option was provided
            if (incld_up2 > len-44) {          //test if incld_up2 is greater than the default cutoff
                end_intrnl_Nrchd = incld_up2;  //extend the enriched lengths to incld_up2
            } else {
                printf("set_enriched_lengths: error - include-up-to value (%d) is less than or equal to the default internal transcript length cutoff (%d). if you are trying to exclude terminally roadblocked transcripts, use the exclude-term option instead\n. aborting...", incld_up2, len-44);
                abort();
            }
        } else {
            end_intrnl_Nrchd = len-44; //use default internal transcript cutoff
        }
        
        strt_trmnl_Nrchd = len - 15; //set start of terminal roadblock enrichment
        end_trmnl_Nrchd = len - 12;  //set end of terminal roadblock enrichment
        
        //set enriched stall site bounds, for use when making rdat files
        mtrx->last_iStl = len-44;
        mtrx->frst_tStl = len-15;
        mtrx->last_tStl = len-12;

        //throw error if whitelisting terminally roadblocked transcripts
        //while also excluding terminally roadblocked transcripts
        if (excld_trmnl && (end_intrnl_Nrchd >= strt_trmnl_Nrchd)) {
            printf("set_enriched_lengths: error - conflicting settings for include-up-to and exclude-term options. include-up-to length %d requires that terminally roadblocked transcripts (lengths %d to %d) are included. aborting...\n", incld_up2, strt_trmnl_Nrchd, end_trmnl_Nrchd);
        }
        
        
        
        for (i = 0; i <= mtrx->tl[MAX]; i++) {
            if (i >= strt_intrnl_Nrchd && i <= end_intrnl_Nrchd) {
                mtrx->Nrchd[i] = 1;
            } else if (i >= strt_trmnl_Nrchd && i <= end_trmnl_Nrchd && !excld_trmnl) {
                mtrx->Nrchd[i] = 1;
            } else {
                mtrx->Nrchd[i] = 0;
            }
        }
        /*
        printf("len = %d\ninclude %d to %d and %d to %d\n", len, strt_intrnl_Nrchd, end_intrnl_Nrchd, strt_trmnl_Nrchd, end_trmnl_Nrchd); */
        
    } else {
        printf("set_enriched_lengths: error - could not find HP4 sequence");
    }
    
    return 1;
}

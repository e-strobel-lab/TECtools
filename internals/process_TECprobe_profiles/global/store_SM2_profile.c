//
//  store_SM2_profile.c
//  
//
//  Created by Eric Strobel on 2/10/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "../../global/global_defs.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "../../utils/gen_utils.h"
#include "../../seq_utils/isRNAbase.h"

#include "store_SM2_profile.h"

/* store_SM2_profile: store shapemapper 2 profile in SM2_profile struct */
void store_SM2_profile(struct SM2_profile * prf, char * filepath)
{
    FILE * fp = NULL;          //file pointer for opening profile file
    char line[MAX_LINE] = {0}; //storage for file lines
    int file_lines = 0;        //file line count
    
    int EOF_reached = 0;     //flag that EOF was reached
    int found_trgt_start = 0; //flag that the start of the target RNA (unmasked sequence) was found
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
    
    char tmp_str[MAX_FIELD] = {0}; //storage for data field strings
    int field = 0;                 //field index
    
    //open the input file to figure out line count
    if ((fp = fopen(filepath, "r")) == NULL) {
        printf("store_SM2_profile: error - could not open %s. Aborting program...\n", filepath);
        abort();
    }
    
    //determine the number of lines in the file
    for (i = 0, EOF_reached = 0; !EOF_reached; i++) {
        if (!get_line_local(line, fp)) {
            EOF_reached = 1;
        }
        
        if (line[0]) {
            file_lines++;
        }
    }
    
    //close the input file
    if (fclose(fp) == EOF) {
        printf("print_input_filenames: error - failed to close output file. Aborting program...\n");
        abort();
    }
        
    line[0] = '\0';                               //zero index [0] of line array
    allocate_SM2_profile_memory(prf, file_lines); //allocate SM2 profile memory using file line count
    
    //open the input file again, this time to parse/store the data
    if ((fp = fopen(filepath, "r")) == NULL) {
        printf("store_SM2_profile: error - could not open %s. Aborting program...\n", filepath);
        abort();
    }
        
    //until the EOF is reached
    for (i = 0, EOF_reached = 0; !EOF_reached; i++) {
        line[0] = '\0'; //zero index [0] of line array
        
        if (!get_line_local(line, fp)) { //get line and set EOF if found
            EOF_reached = 1;
        }
                
        if (line[0]) {                  //if line contains data, process line
            if (!i) {                   //if reading first line of file...
                validate_header(line);  //...validate the header line
            } else {                    //otherwise, process data line
                
                //parse each field from the data line and store
                //the information the SM2 profile structure
                for (j = 0, field = 0; field < PRFL_CLMNS; field++, j++) {
                    
                    //parse field string
                    for (k = 0, tmp_str[0] = '\0'; line[j] != '\t' && line[j]; k++) {
                        tmp_str[k] = line[j++];
                    }
                    tmp_str[k] = '\0';
                    
                      
                    //store field data in SM2_profile struct
                    switch (field) {
                        case NUCLEOTIDE:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->nucleotide[i-1] = atoi(tmp_str);
                            break;
                            
                        case SEQUENCE:
                            check_seq_str(tmp_str);
                            prf->sequence[i-1] = tmp_str[0];
                            prf->tot_nt_cnt++;
                            if (isupper(tmp_str[0])) {
                                prf->trg_nt_cnt++;
                                if (!found_trgt_start) {
                                    prf->trgt_start = i-1;
                                    found_trgt_start = 1;
                                }
                            }
                            break;
                            
                        case MOD_MUTATIONS:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->mod_mutations[i-1] = atoi(tmp_str);
                            break;
                            
                        case MOD_READ_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->mod_read_depth[i-1] = atoi(tmp_str);
                            if (prf->mod_read_depth[i-1] > 0) {
                                prf->chnls.mod = 1;
                            }
                            break;
                            
                        case MOD_EFF_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->mod_eff_depth[i-1] = atoi(tmp_str);
                            break;
                            
                        case MOD_RATE:
                            check_float_str(tmp_str, ABORT_FAILURE);
                            prf->mod_rate[i-1] = strtod(tmp_str, NULL);
                            break;
                            
                        case MOD_OFF_TARGET_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->mod_off_target_depth[i-1] = atoi(tmp_str);
                            break;
                            
                        case MOD_LOW_MAPQ_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->mod_low_mapq_depth[i-1] = atoi(tmp_str);
                            break;
                            
                        case MOD_MAPPED_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->mod_mapped_depth[i-1] = atoi(tmp_str);
                            break;
                            
                        case UNT_MUTATIONS:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->unt_mutations[i-1] = atoi(tmp_str);
                            break;
                            
                        case UNT_READ_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->unt_read_depth[i-1] = atoi(tmp_str);
                            if (prf->unt_read_depth[i-1] > 0) {
                                prf->chnls.unt = 1;
                            }
                            break;
                            
                        case UNT_EFF_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->unt_eff_depth[i-1] = atoi(tmp_str);
                            break;
                            
                        case UNT_RATE:
                            check_float_str(tmp_str, ABORT_FAILURE);
                            prf->unt_rate[i-1] = strtod(tmp_str, NULL);
                            break;
                            
                        case UNT_OFF_TARGET_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->unt_off_target_depth[i-1] = atoi(tmp_str);
                            break;
                            
                        case UNT_LOW_MAPQ_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->unt_low_mapq_depth[i-1] = atoi(tmp_str);
                            break;
                            
                        case UNT_MAPPED_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->unt_mapped_depth[i-1] = atoi(tmp_str);
                            break;
                            
                        case DEN_MUTATIONS:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->den_mutations[i-1] = atoi(tmp_str);
                            break;
                            
                        case DEN_READ_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->den_read_depth[i-1] = atoi(tmp_str);
                            if (prf->den_read_depth[i-1] > 0) {
                                prf->chnls.den = 1;
                            }
                            break;
                            
                        case DEN_EFF_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->den_eff_depth[i-1] = atoi(tmp_str);
                            break;
                            
                        case DEN_RATE:
                            check_float_str(tmp_str, ABORT_FAILURE);
                            prf->den_rate[i-1] = strtod(tmp_str, NULL);
                            break;
                            
                        case DEN_OFF_TARGET_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->den_off_target_depth[i-1] = atoi(tmp_str);
                            break;
                            
                        case DEN_LOW_MAPQ_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->den_low_mapq_depth[i-1] = atoi(tmp_str);
                            break;
                            
                        case DEN_MAPPED_DEPTH:
                            check_int_str(tmp_str, ABORT_FAILURE);
                            prf->den_mapped_depth[i-1] = atoi(tmp_str);
                            break;
                            
                        case REACTIVITY_PROFILE:
                            check_float_str(tmp_str, ABORT_FAILURE);
                            prf->reactivity_profile[i-1] = strtod(tmp_str, NULL);
                            break;
                            
                        case STD_ERR:
                            check_float_str(tmp_str, ABORT_FAILURE);
                            prf->std_err[i-1] = strtod(tmp_str, NULL);
                            break;
                            
                        case HQ_PROFILE:
                            check_float_str(tmp_str, ABORT_FAILURE);
                            prf->hq_profile[i-1] = strtod(tmp_str, NULL);
                            break;
                            
                        case HQ_STDERR:
                            check_float_str(tmp_str, ABORT_FAILURE);
                            prf->hq_stderr[i-1] = strtod(tmp_str, NULL);
                            break;
                            
                        case NORM_PROFILE:
                            check_float_str(tmp_str, ABORT_FAILURE);
                            prf->norm_profile[i-1] = strtod(tmp_str, NULL);
                            break;
                            
                        case NORM_STDERR:
                            check_float_str(tmp_str, ABORT_FAILURE);
                            prf->norm_stderr[i-1] = strtod(tmp_str, NULL);
                            break;
                            
                        default:
                            printf("store_SM2_profile: error - exceeded expected number of profile columns. aborting...\n");
                            abort();
                            break;
                    }
                }
            }
        }
    }
        
    prf->sequence[i-1] = '\0'; //terminate sequence string
        
    //close the input file
    if (fclose(fp) == EOF) {
        printf("print_input_filenames: error - failed to close output file. Aborting program...\n");
        abort();
    }
}


/* get_line_local: get line from file, place into array, remove trailing newline, and return
 line length if successful. local version that allows files to end on non-newline characters */
int get_line_local(char *line, FILE *ifp)
{
    /* function gets line and replaces terminal newline with null character.
     there is no need for buffering in this case because lines that exceed
     MAX_LINE should not exist and if they do, they are an error.
     the only acceptable mode of failure is to reach the end of the file
     without getting any preceeding characters */
    
    int i = 0;
    char c = 0;
    
    for (i = 0; (c = fgetc(ifp)) != '\n' && c != EOF && c &&  i < MAX_LINE; i++) {
        line[i] = c;
    }
    
    if (c == '\n' && i != MAX_LINE) {
        line[i] = '\0';       //remove trailing newline
        return i;             //success
    } else if (c == EOF) {    //reached end of file
        if (i == 0) {         //EOF is expected at the start of a line
            return 0;
        } else {              //last line did not contain newline
            line[i] = '\0';   //append terminating null character
            return 0;
        }
    }else if (!c) {                //unexpected null character
        printf("get_line_local: error - unanticipated null character\n");
        abort();
    } else if (i == MAX_LINE) {    //unexpected long line
        printf("get_line_local: error - unanticipated long line\n");
        abort();
    } else {
        printf("get_line_local: error - reached unreachable code. figure out why. aborting...\n");
        abort();
    }
}


/* validate_header: confirm that column headers match expectations */
void validate_header(char * hdr)
{
    //list of expected headers
    char xpctd_hdrs[PRFL_CLMNS][MAX_NAME] = {
        "Nucleotide",
        "Sequence",
        "Modified_mutations",
        "Modified_read_depth",
        "Modified_effective_depth",
        "Modified_rate",
        "Modified_off_target_mapped_depth",
        "Modified_low_mapq_mapped_depth",
        "Modified_mapped_depth",
        "Untreated_mutations",
        "Untreated_read_depth",
        "Untreated_effective_depth",
        "Untreated_rate",
        "Untreated_off_target_mapped_depth",
        "Untreated_low_mapq_mapped_depth",
        "Untreated_mapped_depth",
        "Denatured_mutations",
        "Denatured_read_depth",
        "Denatured_effective_depth",
        "Denatured_rate",
        "Denatured_off_target_mapped_depth",
        "Denatured_low_mapq_mapped_depth",
        "Denatured_mapped_depth",
        "Reactivity_profile",
        "Std_err",
        "HQ_profile",
        "HQ_stderr",
        "Norm_profile",
        "Norm_stderr"
    };
    
    char crrnt_hdr[MAX_NAME] = {0}; //storage for current header string
    
    int i = 0;   //general purpose index
    int j = 0;   //general purpose index
    int col = 0; //column index
    
    int mtchs = 0;       //number of column header matches
    int hanging_tab = 0; //flag that line contains a hangingtab
    
    for (i = 0, j = 0, col = 0; hdr[i] && hdr[i] != '\r'; i++) {
        if (hdr[i] == '\t') {
            crrnt_hdr[j] = '\0';
            j = 0;
            
            if (hdr[i+1] == '\n' ||  //adjust column count when
                hdr[i+1] == '\r' ||  //there is a hanging tab
                hdr[i+1] == '\0' ) {
                hanging_tab = 1;
            }
            
            
            if (strcmp(crrnt_hdr, xpctd_hdrs[col])) { //if current header does not match expected header, abort
                printf("validate header: error - header does not match expected value.\ncurrent  header: %s\nexpected header: %s\naborting...\n", crrnt_hdr, xpctd_hdrs[col]);
                abort();
            } else {
                mtchs++;
            }
            
            col++; //increment column index
            
        } else {
            crrnt_hdr[j++] = hdr[i]; //copy char to crrnt_hdr
        }
    }
    
    if (!hanging_tab) {       //if the line did not contain a hanging tab
        crrnt_hdr[j] = '\0';  //terminate the input values string
        
        if (strcmp(crrnt_hdr, xpctd_hdrs[col])) { //if last column header does not match expected header, abort
            printf("validate header: error - header does not match expected value.\ncurrent  header: %s\nexpected header: %s\naborting...\n", crrnt_hdr, xpctd_hdrs[col]);
            abort();
        } else {
            mtchs++; //increment match count
        }
        
        col++; //increment column index
    }
        
    if (mtchs != PRFL_CLMNS) { //test whether column count matches expected number
        printf("merge_profiles: error - unexpected column number (%d). aborting...\n", col);
        abort();
    }
}


void allocate_SM2_profile_memory(struct SM2_profile * prf, int data_lines)
{
    int fail = 0;
    
    if ((prf->nucleotide = calloc(data_lines, sizeof(*prf->nucleotide))) == NULL) {
        fail = 1;
    }
    
    if ((prf->sequence = calloc(data_lines+1, sizeof(*prf->sequence))) == NULL) {
        fail = 1;
    }
    
    if ((prf->mod_mutations = calloc(data_lines, sizeof(*prf->mod_mutations))) == NULL) {
        fail = 1;
    }
        
    if ((prf->mod_read_depth = calloc(data_lines, sizeof(*prf->mod_read_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->mod_eff_depth = calloc(data_lines, sizeof(*prf->mod_eff_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->mod_rate = calloc(data_lines, sizeof(*prf->mod_rate))) == NULL) {
        fail = 1;
    }
        
    if ((prf->mod_off_target_depth = calloc(data_lines, sizeof(*prf->mod_off_target_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->mod_low_mapq_depth = calloc(data_lines, sizeof(*prf->mod_low_mapq_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->mod_mapped_depth = calloc(data_lines, sizeof(*prf->mod_mapped_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->unt_mutations = calloc(data_lines, sizeof(*prf->unt_mutations))) == NULL) {
        fail = 1;
    }
        
    if ((prf->unt_read_depth = calloc(data_lines, sizeof(*prf->unt_read_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->unt_eff_depth = calloc(data_lines, sizeof(*prf->unt_eff_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->unt_rate = calloc(data_lines, sizeof(*prf->unt_rate))) == NULL) {
        fail = 1;
    }
        
    if ((prf->unt_off_target_depth = calloc(data_lines, sizeof(*prf->unt_off_target_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->unt_low_mapq_depth = calloc(data_lines, sizeof(*prf->unt_low_mapq_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->unt_mapped_depth = calloc(data_lines, sizeof(*prf->unt_mapped_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->den_mutations = calloc(data_lines, sizeof(*prf->den_mutations))) == NULL) {
        fail = 1;
    }
        
    if ((prf->den_read_depth = calloc(data_lines, sizeof(*prf->den_read_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->den_eff_depth = calloc(data_lines, sizeof(*prf->den_eff_depth))) == NULL) {
        fail = 1;
    }
    
    if ((prf->den_rate = calloc(data_lines, sizeof(*prf->den_rate))) == NULL) {
        fail = 1;
    }
        
    if ((prf->den_off_target_depth = calloc(data_lines, sizeof(*prf->den_off_target_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->den_low_mapq_depth = calloc(data_lines, sizeof(*prf->den_low_mapq_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->den_mapped_depth = calloc(data_lines, sizeof(*prf->den_mapped_depth))) == NULL) {
        fail = 1;
    }
        
    if ((prf->reactivity_profile = calloc(data_lines, sizeof(*prf->reactivity_profile))) == NULL) {
        fail = 1;
    }
        
    if ((prf->std_err = calloc(data_lines, sizeof(*prf->std_err))) == NULL) {
        fail = 1;
    }
        
    if ((prf->hq_profile = calloc(data_lines, sizeof(*prf->hq_profile))) == NULL) {
        fail = 1;
    }
        
    if ((prf->hq_stderr = calloc(data_lines, sizeof(*prf->hq_stderr))) == NULL) {
        fail = 1;
    }
        
    if ((prf->norm_profile = calloc(data_lines, sizeof(*prf->norm_profile))) == NULL) {
        fail = 1;
    }
        
    if ((prf->norm_stderr = calloc(data_lines, sizeof(*prf->norm_stderr))) == NULL) {
        fail = 1;
    }
    
    if ((prf->recalc_norm_profile = calloc(data_lines, sizeof(*prf->recalc_norm_profile))) == NULL) {
        fail = 1;
    }
        
    if ((prf->dataset_norm_profile = calloc(data_lines, sizeof(*prf->dataset_norm_profile))) == NULL) {
        fail = 1;
    }
        
    if (fail) {
        printf("allocate_SM2_profile_memory: error - memory allocation failed. aborting...\n");
        abort();
    }
}

/* check_seq_str: check that sequence string is an RNA base and only a single character */
void check_seq_str(char * str)
{
    if (isRNAbase(str[0]) && !str[1]) {
        return;
    } else {
        printf("check_seq_str: unexpected format for profile sequence string. aborting...\n");
        abort();
    }
}


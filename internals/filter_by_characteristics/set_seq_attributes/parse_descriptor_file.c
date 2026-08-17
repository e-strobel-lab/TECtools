//
//  parse_descriptor_file.c
//  
//
//  Created by Eric Strobel on 11/12/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"
#include "../../utils/gen_utils.h"
#include "../../seq_utils/isIUPACbase.h"

#include "parse_descriptor_file.h"

/* parse_descriptor_file: parse descriptor file and store contents in descriptor structure */
int parse_descriptor_file(char * des_path, descriptor ** des, int * typ_cnt)
{
    int i = 0; //general purpose index
    
    FILE * dfp = NULL; //descriptor file pointer
    
    int xpctd = 0; //expected number of descriptor file lines
        
    //count number of descriptor file lines
    count_descriptors(des_path, &xpctd, typ_cnt);
    
    //allocate memory for descriptors
    if ((*des = calloc(xpctd, sizeof(**des))) == NULL) {
        printf("parse_descriptor_file: error - failed to allocate memory for descriptors. aborting...");
        abort();
    }
    
    //open descriptor file for parsing
    if ((dfp = fopen(des_path, "r")) == NULL) {
        printf("parse_descriptor_file: error - failed to open descriptor file. aborting...\n");
        abort();
    }
    
    char line[MAX_LINE+1] = {0}; //line storage
    
    char * substring = {0}; //pointer to descriptor line substring
    int typ = TYPE_INIT;    //descriptor type, used to simplify code below
        
    //process descriptor file line
    for (i = 0; get_line(line, dfp) && i < xpctd; i++) {
        
        //split line into substrings and store in action structure
        substring = store_field(&((*des)[i].typ_s), line);      //split descriptor type
        substring = store_field(&((*des)[i].nm), substring);    //split descriptor name
        substring = store_field(&((*des)[i].act_s), substring); //split descriptor action
        substring = store_field(&((*des)[i].val_s), substring); //split descriptor value
        
        //check that the descriptor line did not contain any additional fields
        if (substring[0]) {
            printf("parse_descriptor_file: error - descriptor line contains more than the expected number of fields. aborting...\n");
            abort();
        }
        
        printf("\n--------------\n\n%s\n%s\n%s\n%s\n", (*des)[i].nm, (*des)[i].typ_s, (*des)[i].act_s, (*des)[i].val_s);
    
        //parse action field and set action variables
        typ = (*des)[i].typ = ret_type((*des)[i].typ_s);                           //set descriptor type
        parse_descriptor_action((*des)[i].typ, (*des)[i].act_s, &((*des)[i].act)); //parse descriptor action
        
        //parse value field and set value variables
        if (typ == NUC_ID) {
            parse_nuc_iden((*des)[i].val_s, &((*des)[i].sq), &((*des)[i].wndw[SEG1]), &((*des)[i].wndw_cnt));
            //printf("\n%d %s\n", ((*des)[i].wndw[SEG1].b1), ((*des)[i].sq));
            
        } else if (typ == PRX_DG || typ == DST_DG || typ == SS_LEN) {
            (*des)[i].wndw_cnt = parse_window((*des)[i].val_s, &((*des)[i].wndw[SEG1]));
            //printf("wndw1: %d\t%d\n", (*des)[i].wndw[SEG1].b1, (*des)[i].wndw[SEG1].b2);
            //printf("wndw2: %d\t%d\n", (*des)[i].wndw[SEG2].b1, (*des)[i].wndw[SEG2].b2);
            
            //check that the number of values matches the expected number for the descriptor
            //proximal dG and subseq len descriptors must have 2 values that define a single window
            //distal dG must have 4 values that define two windows
            if ((typ == PRX_DG || typ == SS_LEN) && (*des)[i].wndw_cnt != 1) {
                printf("parse_descriptor: unexpected format for descriptor line. aborting...\n");
                abort();
            } else if (typ == DST_DG && (*des)[i].wndw_cnt != 2) {
                printf("parse_descriptor: unexpected format for descriptor line. aborting...\n");
                abort();
            }
        }
        print_descriptor(&(*des)[i]);
    }
        
    //close descriptor file
    if (fclose(dfp) == EOF) {
        printf("parse_descriptor_file: error - failed to close descriptor file. aborting...\n");
        abort();
    }
        
    return i;
}

/* count_descriptors: count number of lines in descriptor file */
void count_descriptors(char * des_path, int * xpctd, int * typ_cnt)
{
    FILE * dfp = NULL; //descriptor file pointer
    
    char line[MAX_LINE+1] = {0}; //line storage
    int i = 0;                   //general purpose index
    
    //open descriptor file
    if ((dfp = fopen(des_path, "r")) == NULL) {
        printf("parse_descriptor_file: error - failed to open descriptor file. aborting...\n");
        abort();
    }
    
    //count number of each descriptor type
    while (get_line(line, dfp)) {
        (*xpctd)++; //increment number of expected descriptors
        
        for (i = 0; !isspace(line[i]) && line[i] && i < MAX_LINE; i++) {;} //iterate to the end of the descriptor type
        line[i] = '\0';                                                    //terminate descriptor type string
        typ_cnt[ret_type(line)]++;                                         //increment count of descriptor type
    }
    
    //close descriptor file
    if (fclose(dfp) == EOF) {
        printf("parse_descriptor_file: error - failed to close descriptor file. aborting...\n");
        abort();
    }
    
    return;
}

/* store_field: store field from descriptor line and return pointer to next field */
char * store_field(char ** field, char * str)
{
    int i = 0;                        //general purpose index
    char tmp_field[MAX_LINE+1] = {0}; //temp storage for field string
        
    //store field string in tmp_field
    for (i = 0; !isspace(str[i]) && str[i] && i < MAX_LINE; i++) {
        tmp_field[i] = str[i];
    }
    tmp_field[i] = '\0';
        
    //allocate memory for field string storage
    if ((*field = malloc((strlen(tmp_field)+1) * sizeof(**field))) == NULL) {
        printf("store_field: error - failed to allocate memory for storing descriptor field. aborting...\n");
        abort();
    }
    
    strcpy(*field, tmp_field); //store field string
        
    if (i == MAX_LINE) { //check that MAX_LINE was not reached
        printf("store_field: error - field value is too long. aborting...\n");
        abort();
    } else if ((*field)[0] == '\0') { //check that field contained information
        printf("store_field: error - field contained no value\n. aborting...\n");
        abort();
    } else {
        while (isspace(str[i])) {i++;} //bypass trailing spaces
        return &str[i];                //return pointer to next non-space character
    }
}

/* ret_type: return descriptor type value */
int ret_type(char * str)
{
    char nt_identity_typ[12] = "nt_identity"; //nucleotide identity type string
    char prox_hlx_dG_typ[12] = "prox_hlx_dG"; //proximal helix deltaG type string
    char dist_hlx_dG_typ[12] = "dist_hlx_dG"; //distal helix deltaG type string
    char sub_seq_len_typ[12] = "sub_seq_len"; //subsequence length type string
    
    //check if input string matches any of the type strings defined above
    //if match is found, return corresponding integer
    //if no match is found, throw error and abort
    if (!strcmp(str, nt_identity_typ)) {
        return NUC_ID;
        
    } else if (!strcmp(str, prox_hlx_dG_typ)) {
        return PRX_DG;
        
    } else if (!strcmp(str, dist_hlx_dG_typ)) {
        return DST_DG;
        
    } else if (!strcmp(str, sub_seq_len_typ)) {
        return SS_LEN;
        
    } else {
        printf("ret_type: error - unrecognized descriptor type in line:\n%s\naborting...\n", str);
        abort();
    }
}

/* parse_descriptor action: parse descriptor file action string and store information in action structure */
void parse_descriptor_action(int typ, char * str, des_action * act)
{
    static int ptyp_all_set = 0; //flag that prediction type all was set for a descriptor
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
        
    char * substr = NULL; //pointer to substring
     
    char include_act[9] = "include:"; //include action string
    char exclude_act[9] = "exclude:"; //exclude action string
    char split_act[7] = "split:";     //split action string
    char bin_act[5] = "bin:";         //bin action string
    
    //determine action by comparing the head of each string to the action strings defined above
    //if match is found, set action and set pointer to character after the action string.
    //if no match is found, throw error and abort
    if (!memcmp(str, include_act, strlen(include_act))) {
        act->act = INCLUDE;
        substr = &str[strlen(include_act)];
        
    } else if (!memcmp(str, exclude_act, strlen(exclude_act))) {
        act->act = EXCLUDE;
        substr = &str[strlen(exclude_act)];
    
    } else if (!memcmp(str, split_act, strlen(split_act))) {
        act->act = SPLIT;
        substr = &str[strlen(split_act)];
        
    }  else if (!memcmp(str, bin_act, strlen(bin_act))) {
        act->act = BIN;
        substr = &str[strlen(bin_act)];
        
    } else {
        printf("parse_descriptor_action: error - unrecognized action. aborting...\n");
        abort();
    }
        
    //prediction type strings
    char mfe_ptyp[5] = "mfe:";  //minimum free energy prediction type
    char all_ptyp[5] = "all:";  //all structures prediction type
    
    //action operator strings
    char is_match_op[6] = "match";                   //match operator
    char not_match_op[10] = "not_match";             //not match operator
    
    char is_equal2_op[2] = "=";                      //equal to operator
    char not_eql_2_op[3] = "!=";                     //not equal to operator
    char less_than_op[2] = "<";                      //less than operator
    char lessOReql_op[3] = "<=";                     //less than or equal to operator
    char grtr_than_op[2] = ">";                      //greater than operator
    char grtrOReql_op[3] = ">=";                     //greater than or equal to operator
    
    char eql_width_op[13] = "equal_width:";          //equal width binning operator
    char eql_frqnc_op[17] = "equal_frequency:";      //equal frequency binning operator
    
    char hlx_pttrn_op[15] = "helix_pattern:";        //helix pattern operator
    
    char stochastic_val[11] = "stochastic";           //stochastic value string
    char min_free_energy_val[16] = "min_free_energy"; //minimum free energy value string
    
    char * tmp_val[2] = {NULL}; //pointer to value strings
    
    int fnd_end = 0; //flag that end of string was found
    
    if (typ == NUC_ID) { //descriptor type is nucleotide identity
        
        //nucleotide identity actions can only be INCLUDE, EXCLUDE, or SPLIT
        //and must use either the IS_MATCH or NOT_MATCH operators
        
        //TODO: not currently allowing split, revisit
        
        if (act->act == INCLUDE || act->act == EXCLUDE /*|| act->act == SPLIT*/) { //if action is permitted
            
            //check that operator is IS_MATCH or NOT_MATCH,
            //otherwise, throw error and abort
            if (!memcmp(substr, is_match_op, strlen(is_match_op))) {
                act->op[0] = IS_MATCH;
                substr = &substr[strlen(is_match_op)];
                act->op_cnt++;
                
            } else if (!memcmp(substr, not_match_op, strlen(not_match_op))) {
                act->op[0] = NOT_MATCH;
                substr = &substr[strlen(not_match_op)];
                act->op_cnt++;
        
            } else {
                printf("parse_descriptor_action: error - operator is incompatible with nucleotide identity descriptor type. aborting...\n");
                abort();
            }
            
            //check that no additional characters are present after the operator
            if (substr[i]) {
                printf("parse_descriptor_action: error - nucleotide identity descriptor strings should only include an action and an operator, not a value. aborting...\n");
                abort();
                
            } else {
                act->val_cnt = 0; //set value count to 0
                return;           //action processed, return
            }
            
        } else { //if action is not permitted, throw error and abort
            printf("parse_descriptor_action: error - action is incompatible with nucleotide identity descriptor type. aborting...\n");
            abort();
        }
        
    } else if (act->act == BIN) { //action type is bin
                
        if (!memcmp(substr, mfe_ptyp, strlen(mfe_ptyp))) { //if prediction type is minimum free energy
            act->ptyp = MFE_STRUCT;                        //set action ptyp to MFE_STRUCT
            substr = &substr[strlen(mfe_ptyp)];            //bypass mfe_ptyp substring
            
        } else if (!memcmp(substr, all_ptyp, strlen(all_ptyp))) { //if prediction type is all
            if (!ptyp_all_set) {                                  //if no other prediction type was already set to all
                act->ptyp = ALL_STRUCTS;                          //set action ptyp to ALL_STRUCTS
                substr = &substr[strlen(all_ptyp)];               //bypass all_ptyp substring
                ptyp_all_set = 1;                                 //set flag that a ptyp was set to all
                
            } else { //if more than one prediction typ was set to all, throw error and abort
                printf("parse_descriptor_action: error - more than one descriptor contains prediction type 'all'. aborting...\n");
                abort();
            }
            
        } else { //unrecognized prediction type, throw error and abort
            printf("parse_descriptor_action: error - prediction type is incompatible with bin action. aborting...\n");
            abort();
        }
        
        //check that operator is EQL_WIDTH, EQL_FRQNC, or HLX_PTTRN
        //otherwise, throw error and abort
        if (!memcmp(substr, eql_width_op, strlen(eql_width_op))) {
            act->op[0] = EQL_WIDTH;
            substr = &substr[strlen(eql_width_op)];
            act->op_cnt++;
            
        } else if (!memcmp(substr, eql_frqnc_op, strlen(eql_frqnc_op))) {
            act->op[0] = EQL_FRQNC;
            substr = &substr[strlen(eql_frqnc_op)];
            act->op_cnt++;
    
        } else if (!memcmp(substr, hlx_pttrn_op, strlen(hlx_pttrn_op))) {
            act->op[0] = HLX_PTTRN;
            substr = &substr[strlen(hlx_pttrn_op)];
            act->op_cnt++;
            
        } else {
            printf("parse_descriptor_action: error - operator is incompatible with bin action. aborting...\n");
            abort();
        }
        
        //check additional characters are present after the operator
        if (substr[0]) {
            //note: value string is parsed later
            tmp_val[0] = substr; //set tmp_val[0] pointer to the start of substring
            act->val_cnt++;      //increment value count
        } else {
            printf("parse_descriptor_action: error - descriptor action 'bin' requires a value following the operator. aborting...\n");
            abort();
        }
                
    } else if ((typ == PRX_DG       || typ == DST_DG       || typ == SS_LEN) &&
               (act->act == INCLUDE || act->act == EXCLUDE || act->act == SPLIT)) {
        
        //proximal dG, distal dG, and subseq length descriptors with actions
        //INCLUDE, EXCLUDE, and SPLIT are all processed in the same way
        
        //these actions may include up to two operators and associated values,
        //which are parsed by the loop below
        
        for (i = 0, fnd_end = 0; i < 2 && !fnd_end; i++) {
            
            //check that the head of the substring matches one of the operator
            //strings defined above, increment action operator count
            
            if (!memcmp(substr, is_equal2_op, strlen(is_equal2_op))) {
                act->op[i] = IS_EQUAL2;
                substr = &substr[strlen(is_equal2_op)];
                act->op_cnt++;
                
            } else if (!memcmp(substr, not_eql_2_op, strlen(not_eql_2_op))) {
                act->op[i] = NOT_EQL_2;
                substr = &substr[strlen(not_eql_2_op)];
                act->op_cnt++;
                
            } else if (!memcmp(substr, lessOReql_op, strlen(lessOReql_op))) {
                act->op[i] = LESSorEQL;
                substr = &substr[strlen(lessOReql_op)];
                act->op_cnt++;
                
            } else if (!memcmp(substr, grtrOReql_op, strlen(grtrOReql_op))) {
                act->op[i] = GRTRorEQL;
                substr = &substr[strlen(grtrOReql_op)];
                act->op_cnt++;
                
            } else if (!memcmp(substr, less_than_op, strlen(less_than_op))) {
                act->op[i] = LESS_THAN;
                substr = &substr[strlen(less_than_op)];
                act->op_cnt++;
                
            } else if (!memcmp(substr, grtr_than_op, strlen(grtr_than_op))) {
                act->op[i] = GRTR_THAN;
                substr = &substr[strlen(grtr_than_op)];
                act->op_cnt++;
                
            } else {
                printf("parse_descriptor_action: unrecognized operator. aborting...\n");
                abort();
            }
            
            //parse value for current operator
            tmp_val[i] = substr;                                //point tmp_val[i] to substring
            for (j = 0; substr[j] && substr[j] != ','; j++) {;} //iterate to null char or comma
            
            if (j && !substr[j]) { //if j was incremented and null char was reached
                act->val_cnt++;    //increment action value count
                fnd_end = 1;       //set fnd_end flag to indicate that the end of the string was found
                
            } else if (j && substr[j] == ',') { //if j was incremented and a comma was reached
                if (!i) {                  //if in the first loop iteration
                    substr[j] = '\0';      //substitute the comma with a null char to terminate tmp_val[0]
                    substr = &substr[j+1]; //set pointer to the char after the comma
                    act->val_cnt++;        //set number of values detected in the current action
                    
                } else { //if comma was found on >0 iteration, throw error and abort
                    printf("parse_descriptor_action: error - unexpected format for action substring. aborting...\n");
                    abort();
                }
                
            } else { //it should not be possible to reach this error
                printf("parse_descriptor_action: descriptor operator that requires an associated value did not contain a value. aborting...");
                abort();
            }
        }
    } else {
        printf("parse_descriptor_action: error - unexpected input in descriptor file. aborting...\n");
        abort();
    }
    
    //check that value count matches operator count. when parsing descriptors where
    //value and operator counts are not supposed to match, return is called before
    //reaching this check
    if (act->val_cnt != act->op_cnt) {
        printf("parse_descriptor_action: action value count must match action operator count. aborting...\n");
        abort();
    }
    
    //convert value strings to integer or double, depending on descriptor type and action
    //descriptor type 'subseq length' and action 'bin' (excluding bin actions for helix
    //pattern type) are always associated with integer values. all other types/actions
    //that are associated with values use double
    
    char * endptr = NULL; //pointer for error checking strtol and strtod
    errno = 0;            //initialize errno to 0
    
    //if subseq_len type or bin action with equal width or freq operator, store value as integer
    if (typ == SS_LEN || (act->act == BIN && (act->op[0] == EQL_WIDTH || act->op[0] == EQL_FRQNC))) {
        
        act->val_typ = LONG_INTEGER; //set value type to LONG_INTEGER
        
        //allocate memory for storing value(s)
        if ((act->val = calloc(act->val_cnt, sizeof(long))) == NULL) {
            printf("parse_descriptor_action: error - failed to allocate memory for action value string. aborting...\n");
            abort();
        }
        
        for (i = 0; i < act->val_cnt; i++) {   //for each value string
            for (j = 0; tmp_val[i][j]; j++) {  //check that all chars are digits,
                if (!isdigit(tmp_val[i][j])) { //if not, throw error and abort
                    printf("parse_descriptor_action: error - integer value string must be composed of digits. aborting...\n");
                    abort();
                }
            }
            
            ((long *)act->val)[i] = strtol(tmp_val[i], &endptr, 10); //store value as long int
            
            if (endptr == &tmp_val[i][0]) { //if conversion failed, throw error and abort
                printf("parse_descriptor_action: error - conversion of value string to long integer failed. aborting...\n");
                abort();
            }
        }
        
    } else if (act->act == BIN && act->op[0] == HLX_PTTRN) { //if action is BIN and operation is HLX_PATTERN
        if (!strcmp(tmp_val[0], min_free_energy_val)) {      //if value type is minimum free energy
            act->val_typ = MIN_FREE_NRG;                     //set value type to MIN_FREE_NRG
            
        } else if (!strcmp(tmp_val[0], stochastic_val)) {    //if falue type is stochastic
            act->val_typ = STOCHASTIC;                       //set value type to STOCHASTIC
        
        } else { //unrecognized value type, throw error and abort
            printf("parse_descriptor_file: error - unrecognized value for helix_pattern operator. expected either 'min_free_enery' or 'stochastic'. aborting...\n");
            abort();
        }
        
    } else { //for all other types/actions that are associated with a value, store value as double
        
        act->val_typ = FLOATING_PT; //set value type to floating point (double)
        
        //allocate memory for values
        if ((act->val = calloc(act->val_cnt, sizeof(double))) == NULL) {
            printf("parse_descriptor_action: error - failed to allocate memory for action value string. aborting...\n");
            abort();
        }
        
        for (i = 0; i < act->val_cnt; i++) {    //for each value string
            for (j = 0; tmp_val[i][j]; j++) {   //check that chars are digits, '-' or '.'
                
                //TODO: could allow only one '.' and require that '-' is the first char
                
                if (!isdigit(tmp_val[i][j])  && //if not, throw error and abort
                    tmp_val[i][j] != '-'     &&
                    tmp_val[i][j] != '.') {
                    printf("parse_descriptor_action: error - floating point value string must be composed of digits, '.' or '-' chars. aborting...\n");
                    abort();
                }
            }
            
            ((double *)act->val)[i] = strtod(tmp_val[i], &endptr); //store value as double
            
            if (endptr == &tmp_val[i][0]) { //if conversion failed throw error and abort
                printf("parse_descriptor_action: error - conversion of value string to long integer failed. aborting...\n");
                abort();
            }
        }
    }
    
    printf("\n");
    return;
}

/* parse_nuc_iden: parse nucleotide identity string */
void parse_nuc_iden(char * str, char ** sq, nt_window * wndw, int * wndw_cnt)
{
    //NOTE: for nucleotide identity descriptors, bound 1 is the index at
    //which the string starts and bound 2 is the length of the string
    
    int v = 0; //value index
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int len = strlen(str);          //set input string length
    char tmp_str[MAX_LINE+1] = {0}; //array for storing substring
    
    if (len > MAX_LINE) { //throw error if input string is too long to store in tmp_str
        printf("parse_nuc_iden: error - input string is too long. aborting...\n");
        abort();
    }
    
    //parse nucleotide identity string
    for (v = 0, i = 0; v < 2; v++) { //for each value up to a maximum of 2
        
        for (j = 0; str[i] != ':' && str[i] && i < len; i++, j++) { //read string until terminating char
            
            //in iteration 0, throw error if char is not a digit
            if (v == 0 && !isdigit(str[i])) {
                printf("parse_nuc_iden: error - nucleotide position must be composed of digits. aborting...\n");
                abort();
        
            //in iteration 1, throw error if char is not an IUPAC base
            } else if (v == 1 && !isIUPACbase(str[i])) {
                printf("parse_nuc_iden: error - sequence must be composed of IUPAC bases. aborting...\n");
                abort();
            
            //if char is expected type, copy to tmp_str
            } else {
                tmp_str[j] = str[i];
            }
        }
        tmp_str[j] = '\0'; //terminate tmp_str
        
        if (v == 0) { //in iteration 0, check that loop ended on a ':' and that there are more chars
            if (str[i] != ':' || i == len) { //if not, throw error and abort
                printf("parse_nuc_iden: error- unexpected format for nucleotide identity descriptor. aborting...\n");
                abort();
            } else {                        //if so
                wndw->b1 = atoi(tmp_str);   //set bound 1 to int encoded by tmp_str
                (*wndw_cnt)++;              //increment window count
                i++;                        //increment i
            }
            
        } else if (v == 1) { //in iteration 1, check that loop ended on'\0' and that there are no more chars
            if (str[i] || i != len) { //if not, throw error and abort
                printf("parse_nuc_iden: error- unexpected format for nucleotide identity descriptor. aborting...\n");
                abort();
                
            } else { //if so, allocate mem for sq string and copy tmp_str to sq
                if ((*sq = malloc((strlen(tmp_str)+1) * sizeof(**sq))) == NULL) {
                    printf("parse_nuc_iden: error - failed to allocate memory for nucleotide identity descriptor. aborting...\n");
                    abort();
                }
                strcpy(*sq, tmp_str);
                wndw->b2 = strlen(tmp_str); //set bound 2 to string length
                return;
            }
        }
    }
    
    //if there are additional characters after the expected values, throw error and abort
    //(this may not be reachable)
    printf("parse_nuc_iden: error - more fields than expected in nucleotide identity string. aborting...\n");
    abort();
}

/* parse_window: parse window field in descriptor file line */
int parse_window(char * str, nt_window * wndw)
{
    int w = 0; //window index
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int len = strlen(str); //set input string length
    
    char tmp_str[MAX_LINE+1] = {0}; //array for storing substring
    int  tmp_val[4] = {0};          //array for storing substring values
    int  wndw_bnd_cnt = 0;          //window bound count
    
    if (len > MAX_LINE) {
        printf("parse_nuc_iden: error - input string is too long. aborting...\n");
        abort();
    }
    
    //parse window bound strings. the loop always runs for four iterations,
    //which is the maximum number of window bounds. when there are fewer than
    //four window bounds, the latter two remain zero
    
    for (w = 0, i = 0; w < 4; w++) {

        //read and store window bound in tmp_str
        for (j = 0; isdigit(str[i]) && str[i] && i < len; i++, j++) {
            tmp_str[j] = str[i];
        }
        tmp_str[j] = '\0';
        
        if (str[i] == ',' || str[i] == ';' || !str[i]) { //if a window bound terminator was reached,
            tmp_val[wndw_bnd_cnt++] = atoi(tmp_str);     //convert tmp_str to integer and store in tmp_val
            
            if (!str[i]) {                  //if terminating null char was reached
                wndw[SEG1].b1 = tmp_val[0]; //store all window bound values. if
                wndw[SEG1].b2 = tmp_val[1]; //only one window was encoded, the latter
                wndw[SEG2].b1 = tmp_val[2]; //two values will be zero
                wndw[SEG2].b2 = tmp_val[3];
                
                //TODO: decide how string is stored and make sure indexing is done correctly.
                //TODO: with the current code, the string should start at index 1 with some placeholder at position 0 of the array
                
                if (wndw_bnd_cnt % 2) {
                    printf("parse_window: error - expected window bound count to be even. aborting...\n");
                    abort();
                } else {
                    return wndw_bnd_cnt/2; //return the number of windows
                }

            } else { //otherwise
                i++; //increment past the terminating ',' or ';' char
            }
            
        } else {
            printf("parse_window: error - unexpected window descriptor value. aborting...\n");
            abort();
        }
    }
    
    printf("parse_window: error - more fields than expected in window string. aborting...\n");
    abort();
}

/* print_descriptor: print descriptor parameters */
void print_descriptor(descriptor * des)
{
    int i = 0; //general purpose index
    
    //print descriptor type
    switch (des->typ) {
        case NUC_ID: printf("nt_identity\t"); break;
        case PRX_DG: printf("prox_hlx_dG\t"); break;
        case DST_DG: printf("dist_hlx_dG\t"); break;
        case SS_LEN: printf("sub_seq_len\t"); break;
        default:
            printf("print_descriptor: unrecognized descriptor type. aborting...\n");
            abort();
            break;
    }
    
    //print descriptor action
    switch (des->act.act) {
        case INCLUDE: printf("include:"); break;
        case EXCLUDE: printf("exclude:"); break;
        case SPLIT: printf("split:"); break;
        case BIN: printf("bin:"); break;
        default:
            printf("print_descriptor: unrecognized descriptor action. aborting...\n");
            break;
    }
    
    //print descriptor operators
    for (i = 0; i < des->act.op_cnt; i++) {
        
        //print current operator
        switch (des->act.op[i]) {
            case IS_MATCH: printf("match\t"); break;
            case NOT_MATCH: printf("not_match\t"); break;
            case STOCH_FLD: printf("stochastic\t"); break;
            case MFNRG_FLD: printf("min_free_energy\t"); break;
            case IS_EQUAL2: printf("="); break;
            case NOT_EQL_2: printf("!="); break;
            case GRTRorEQL: printf(">="); break;
            case LESSorEQL: printf("<="); break;
            case GRTR_THAN: printf(">"); break;
            case LESS_THAN: printf("<"); break;
            case EQL_WIDTH: printf("equal_width:"); break;
            case EQL_FRQNC: printf("equal_frequency:"); break;
                
            default:
                printf("print_descriptor: error - unrecognized operator. aborting...\n");
                abort();
                break;
        }
        
        //print operator value
        if (i < des->act.val_cnt) {                           //if there is a value to print
            if (des->act.val_typ == LONG_INTEGER) {           //if the value type is long int
                printf("%ld", ((long *)des->act.val)[i]);     //print value as a long integer
                
            } else if (des->act.val_typ == FLOATING_PT) {     //if the value type is floating point
                printf("%5.2f", ((double *)des->act.val)[i]); //print the value as a double
            
            } else { //unrecognized value type, throw error and abort
                printf("print_descriptor: error - unrecognized value type. aborting...\n");
                abort();
                break;
            }
            
            if (i+1 < des->act.op_cnt) { //if there is another operator to print
                printf(",");             //print a comma
            } else {
                printf("\t");            //otherwise print a tab
            }
        }
    }
    
    //print window(s)
    for (i = 0; i < des->wndw_cnt; i++) {                   //for each window
        printf("%d,%d", des->wndw[i].b1, des->wndw[i].b2);  //print the window bounds
        if (i+1 < des->wndw_cnt) {                          //if there is another window
            printf(";");                                    //print a semi-colon
        } else if (des->typ == NUC_ID) {                    //if descriptor type is nucleotide identity
            printf("\t%s\n", des->sq);                      //print the sequence
        } else {
            printf("\n");                                   //if nothing else to print, print a newline
        }
    }
    
    return;
}



//
//  filter_values.c
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#include <stdio.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"
#include "../seq_utils/ispair.h"
#include "../seq_utils/basemap.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/isDNAbase.h"

#include "../TECdisplay_mapper/TECdisplay_output_column_headers.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"
#include "./read_vbase.h"
#include "./search_4_vbase_match.h"

#include "filter_values.h"

/* filter_values: read input values file and send variants that match each set of
 constraints to the respective output file */
int filter_values(FILE * ipt, constraints * cons, int cons_cnt, basemap * bmap, char * out_dir_nm, int nonstandard, char * out_prefix) {
        
    int i = 0; //general purpose index
    int vi = 0; //variant index, counts the number of variants in values file
    
    //output variables
    FILE * ofp[MAX_CONSTRAINTS];          //output file pointers
    char tmp_nm[2048] = {0};              //array to store name of output files
    
    //general input variables
    int smpl_cnt = 0;                     //number of samples in merged values file
    char line[MAX_LINE+1] = {0};          //array to store value file line
    char tmp_line[MAX_LINE+1] = {0};      //array for storing temporary copy of line
    
    //vbase/values input variables
    char *p_id = NULL;                    //pointer to start of variant id in line array
    char *p_vals = NULL;                  //pointer to start of variant values in line array
    int tot_vals = 0;                     //number of tab-separated value fields observed in the line
    sample_values vals[MAX_SAMPLE] = {0}; //array of variant names. NOTE: this structure was previously
                                          //used to store data too, but this functionality is not in use
    
    //match tracking variables
    int match[MAX_CONSTRAINTS] = {0};     //array to track which constraints a given variant matches
    int fnd_match = 0;                    //flag that variant is a match to current constraint being assessed
    int match_cnt = 0;                    //the number of variants with at least 1 constraint match
    
    /**** parse first line of values file to get column headers ****/
    get_line(line, ipt); //get first line of values file
    
    //non-standard format data that is organized in excel likely has a
    //carriage return at the end of the string, which should be removed
    if (line[strlen(line)-1] == '\r') {
        line[strlen(line)-1] = '\0';
    }
    
    strcpy(tmp_line, line);
    
    if (!nonstandard) {  //input data is standard TECdisplay format
        smpl_cnt = get_sample_info(tmp_line, &vals[0], &tot_vals); //parse line for column headers
    }
    
    /**** open output file for each constraint ****/
    for (i = 0; i < cons_cnt; i++) {
        
        //generate output file
        if (out_prefix[0]) {
            sprintf(tmp_nm, "%s/%s_%s.txt", out_dir_nm, out_prefix, cons[i].nm);
        } else {
            sprintf(tmp_nm, "%s/%s.txt", out_dir_nm, cons[i].nm);
        }

        if ((ofp[i] = fopen(tmp_nm, "w")) == NULL) {
            printf("filter_values: error - could not open output file. Aborting program...\n");
            abort();
        }
        
        print_header_line(ofp[i], line, out_dir_nm, cons[i].nm, smpl_cnt, nonstandard);
    }
        
    /**** determine if each variant matches any of the supplied constraints ****/
    for (vi = 0; get_line(line, ipt); vi++) { //for every variant in the merged input file
        
        fnd_match = test_variant_match(bmap, cons, line, &p_id, &p_vals, &match[0], cons_cnt, tot_vals, nonstandard);
        
        if (fnd_match) {  //found match
            match_cnt++;  //track number of variants with at least 1 match
            
            //iterate through match array. if match to a constraint was found,
            //print the input line to the corresponding constraints file
            for (i = 0; i < cons_cnt; i++) {
                if (match[i]) {
                    fprintf(ofp[i], "%s\t%s\n", p_id, p_vals); //print match to corresponding constraints file
                }
            }
        }
        
        line[0] = '\0'; //zero the first character of the line array
    }
    
    printf("\n%d variants with >=1 match\n", match_cnt); //print number of variants with at least one match to screen
    
    /* close output files */
    for (i = 0; i < cons_cnt; i++) {
        if (fclose(ofp[i]) == EOF) {
            printf("filter_values: error - error occurred when closing output file %s.txt. Aborting program...\n", cons[i].nm);
            abort();
        }
    }
    
    return 1;
}

/* exclude_matches: read input values file and send variants that do not match any
 constraints to the output file */
int exclude_matches(FILE * ipt, constraints * cons, int cons_cnt, char * cons_nm, basemap * bmap, char * out_dir_nm, int nonstandard, char * out_prefix) {
        
    int i = 0; //general purpose index
    int vi = 0; //variant index, counts the number of variants in values file
    
    //output variables
    FILE * ofp;                           //output file pointers
    char tmp_nm[2048] = {0};              //array to store name of output files
    
    //general input variables
    int smpl_cnt = 0;                     //number of samples in merged values file
    char line[MAX_LINE+1] = {0};          //array to store value file line

    //vbase/values input variables
    char *p_id = NULL;                    //pointer to start of variant id in line array
    char *p_vals = NULL;                  //pointer to start of variant values in line array
    int tot_vals = 0;                     //number of tab-separated value fields observed in the line
    sample_values vals[MAX_SAMPLE] = {0}; //array of variant names. NOTE: this structure was previously
                                          //used to store data too, but this functionality is not in use
    //match tracking variables
    int match[MAX_CONSTRAINTS] = {0};     //array to track which constraints a given variant matches
    int fnd_match = 0;                    //flag that variant is a match to current constraint being assessed
    int nomatch_cnt = 0;                  //the number of variants without a constraint match
    
    /**** parse first line of values file to get column headers ****/
    get_line(line, ipt); //get first line of values file
    
    if (!nonstandard) {  //input data is standard TECdisplay format
        smpl_cnt = get_sample_info(line, &vals[0], &tot_vals); //parse line for column headers
    }
    
    /**** open output file ****/
    
    //generate output file
    if (out_prefix[0]) {
        sprintf(tmp_nm, "%s/%s_%s.txt", out_dir_nm, out_prefix, cons_nm);
    } else {
        sprintf(tmp_nm, "%s/%s.txt", out_dir_nm, cons_nm);
    }
    
    
    if ((ofp = fopen(tmp_nm, "w")) == NULL) {
        printf("filter_values: error - could not open output file. Aborting program...\n");
        abort();
    }
    
    print_header_line(ofp, line, out_dir_nm, cons_nm, smpl_cnt, nonstandard);
        
    /**** determine if each variant matches any of the supplied constraints ****/
    for (vi = 0; get_line(line, ipt); vi++) { //for every variant in the merged input file
        
        fnd_match = test_variant_match(bmap, cons, line, &p_id, &p_vals, &match[0], cons_cnt, tot_vals, nonstandard);
        
        if (!fnd_match) {   //no match was found
            nomatch_cnt++;  //track number of variants with no matches
            fprintf(ofp, "%s\t%s\n", p_id, p_vals); //print variant to output file
        }
        
        line[0] = '\0'; //zero the first character of the line array
    }
    
    printf("\n%d variants with no matches were included in the output file\n", nomatch_cnt); //print number of variants with at least one match to screen
    
    /* close output file */
    if (fclose(ofp) == EOF) {
        printf("filter_values: error - error occurred when closing output file %s.txt. Aborting program...\n", cons[i].nm);
        abort();
    }
    
    return 1;
}

/* print_header_line: print header line to file */
void print_header_line(FILE * ofp, char * line, char * out_dir_nm, char * cons_hdr, int smpl_cnt, int nonstandard)
{
    extern const char TECdsply_clmn_hdrs[4][32]; //column headers from TECdisplay_output_column_headers.c
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
    int l = 0; //general purpose index
    
    char tmp_hdr[MAX_COL_NM+1] = {0};     //array to temporarily store headers
    char * tmp_std_hdr = NULL;            //temp pointer to standard header string
    char * std_hdr_ptr = NULL;            //pointer to standard header string
    int std_hdr_mtchs = 0;                //number of matches to standard header string
    int line_end = 0;                     //flag that end of line was reached
    
    for (j = 0, k = 0, line_end = 0; !line_end; j++) {
        
        if (line[j] && line[j] != '\t') { //store tab-delimited column header
            tmp_hdr[k++] = line[j];
            
        } else {
            tmp_hdr[k] = '\0'; //terminate tmp_hdr string
            
            //check for each standard header
            for (l = 1, std_hdr_mtchs = 0; l < TDSPLY_HDR_CNT && std_hdr_mtchs <= 1; l++) {
                
                //search tmp_hdr for the current standard header string
                if ((tmp_std_hdr = strstr(tmp_hdr, TECdsply_clmn_hdrs[l])) != NULL) {
                    
                    //check that the standard header is preceded by a '_' and followed by a '\0'
                    if (tmp_std_hdr[-1] == '_' && !tmp_std_hdr[strlen(TECdsply_clmn_hdrs[l])]) {
                        
                        //check that the standard header is part of a longer string
                        if (strlen(tmp_hdr) > (strlen(tmp_std_hdr)+1)) {
                            tmp_std_hdr[-1] = '\0';     //change preceding underscore to terminating null
                            std_hdr_ptr = tmp_std_hdr;  //set std_hdr_ptr to tmp_std_hdr
                            std_hdr_mtchs++;            //increment standard header matches
                            
                        } else {
                            printf("exclude_matches: error - expected standard TECdisplay header to be part of a longer string. aborting...\n");
                            abort();
                        }
                    }
                }
            }
            
            //print header
            if (std_hdr_mtchs == 0) {                                     //no standard header was found
                fprintf(ofp, "%s_%s", tmp_hdr, cons_hdr);                 //print column header
                
            } else if (std_hdr_mtchs == 1) {                              //one standard header was found
                fprintf(ofp, "%s_%s_%s", tmp_hdr, cons_hdr, std_hdr_ptr); //print column header
                
            } else { //error, more than one standard header was found
                printf("exclude_matches: error - nonstandard header contains multiple different standard header strings. aborting...");
                abort();
            }
            
            //print tab separator or turn on line end flag
            if (line[j] == '\t') {     //found tab separator
                fprintf(ofp, "\t");    //print tab to file
                
            } else if (!line[j]) {     //found terminal null
                line_end = 1;          //turn on line end flag
                
            } else {
                printf("filter_values: error - this should be unreachable.\n");
            }
            
            k = 0; //reset k to 0 for next column header
        }
    }
    
    fprintf(ofp, "\n"); //print newline to end line
    
    return;
}

int test_variant_match(basemap * bmap, constraints * cons, char * line, char ** p_id, char ** p_vals, int * match, int cons_cnt, int tot_vals, int nonstandard)
{
    int i = 0; //general purpose index
    
    //vbase input variables
    char *p_vbase = NULL;                 //pointer to current vbase string in variant id
    base_params vbases[MAX_VBASES] = {0}; //storage for variable bases in variant ids
    int crnt_vbase_cnt = 0;               //number of variable bases in current variant id
    
    //match tracking variables
    int fnd_match = 0;                    //flag that variant is a match to current constraint being assessed

    //initialize p_id and p_vals to NULL
    *p_id = NULL;
    *p_vals = NULL;
    
    //initialize match array to 1. failed matches are set
    //to zero below in the test_bases function
    for (i = 0; i < cons_cnt; i++) {match[i] = 1;}

    //split input line into id and values strings
    for (i = 0; line[i] != '\t' && line[i] && i < MAX_LINE; i++) {;} //iterate to first tab
    if (line[i] == '\t') {
        line[i] = '\0';              //set first tab to null char to terminate variant id
        *p_id = p_vbase = &line[0];  //set id and vbase pointers to start of line
        *p_vals = &line[i+1];        //set values pointer to the index after the first tab
    } else if (line[i]) {
        printf("test_variant_match: unexpected format for data line in merged values file. aborting...\n");
        abort();
    }
    
    //parse variant id for vbases
    //NOTE: MAX_VBASES is used as the limit instead of SEQ2BIN_MAX_KEY because
    //CONSTANT constraints can be present in variant ids alongside VARIABLE constraints
    for (i = 0; *(p_vbase) && i < MAX_VBASES; i++) {
        p_vbase = read_vbase(p_vbase, '_', &vbases[i], bmap, NAME); //read vbase to delimiter '_'
        if (*p_vbase == '_') {      //reached delimiter, vbase was not the last vbase
            p_vbase = &p_vbase[1];  //set vbase pointer to the start of the next vbase
        }
    }
    crnt_vbase_cnt = i; //set crnt_vbase_cnt to the number of vbases in the current variant id
    
    if (crnt_vbase_cnt == MAX_VBASES && *(p_vbase)) { //check that MAX_VBASES is not exceeded by the variant ID
        printf("test_variant_match: error - maximum number of vbases in variant name (%d) was exceeded. aborting...\n", MAX_VBASES);
        abort();
    }
    
    //if input data is standard TECdisplay format, validate data-containing segment of line
    if (!nonstandard) {
        validate_data_line(*p_id, *p_vals, tot_vals);
    }
    
    //test if variant is a match to a constraint
    return fnd_match = test_bases(*p_id, &vbases[0], crnt_vbase_cnt, cons, cons_cnt, bmap, &match[0]);
    
}

/* validate_data_line: validate data line from merged values file */
void validate_data_line(char * p_id, char * p_vals, int tot_vals)
{
    int i = 0;               //general purpose index
    int crrnt_tot_vals = 0;  //number of value fields observed in the current data line

    /* test that line begins with non-tab character */
    if (p_vals[0] != '\t') { //first char of data line cannot be a tab
        crrnt_tot_vals = 1;  //initialize field value counter to 1, since first char is not a tab
    } else {
        printf("validate_data_line: error - values string begins with a tab character. aborting...\n");
        abort();
    }
    
    /* test that line contains only valid characters and
     count the number of tab-separated columns in the line */
    for (i = 0; p_vals[i] && i < MAX_LINE; i++) { //test every character of the line
        if (p_vals[i] == '\t' && p_vals[i+1]) {   //tab followed by non-null indicates entering next value field
            crrnt_tot_vals++;                     //increment current total values counter
            
        } else if (!isdigit(p_vals[i]) &&  //if not tab, char must be a digit,
                   p_vals[i] != '.'    &&  //or a '.',
                   p_vals[i] != '-') {     //or a '-' //TODO: not sure negatives are actually allowed?
            
            printf("validate_data_line: error - values string contains an unexpected '%c' (%d) character. aborting...\n", p_vals[i], p_vals[i]);
        }
    }
    
    /* check observed number of value fields matches expected number */
    if (crrnt_tot_vals != tot_vals) {
        printf("validate_data_line: error - variant %s contains %d value fields, which does not match the expected number of value fields determined by the column headers (%d). aborting...\n", p_id, crrnt_tot_vals, tot_vals);
        abort();
    }
    
    return;
}

/* test_bases: test if variant matches any input constraints */
int test_bases(char * p_id, base_params * vbases, int vbase_cnt, constraints * cons, int cons_cnt,  basemap * bmap, int * match)
{
    int c = 0;    //index of constraint
    int i = 0;    //index for base/pair settings of constraint
    int j = 0;    //index for member bases of a pair
            
    int match_cnt = 0;  //number of constraint matches observed for the variant
    int fnd_vbMtch = 0; //flag that match was found for current variable base
    int vb_mtchCnt = 0; //counts total number of variable base matches
    
    int uniq_match_tbl[MAX_VBASES] = {0}; //table for tracking unique base constraint matches

    int  mate_valid[2] = {0}; //flag for testing validity of each mate of a pair
    char base[2] = {0};       //array to store pair bases, used to handle variable and constant bases simultaneously
    int  pair_type = 0;       //integer value representing pair type, stores output of ispair function
    
    base_params * tmp_vb = NULL; //used to store search_for_vbase_match return value
    
    /**** test variant to determine if it is a match to each constraint ****/
    for (c = 0; c < cons_cnt; c++) { //for every constraint

        for (j = 0; j < MAX_CONSTRAINTS; j++) {uniq_match_tbl[j] = 0;} //zero unique match table
        
        /**** test if vbases of input variant match base identity constraints ****/
        for (i = 0, vb_mtchCnt = 0; i < vbase_cnt; i++) {
            
            //for every variable base in the variant id, check
            //for a match to a constrained variable base
            
            for (j = 0, fnd_vbMtch = 0; !fnd_vbMtch && j < cons[c].bcnt && j < MAX_VBASES; j++) {
                
                if (vbases[i].pos == cons[c].base[j].pos &&              //test position match
                    vbases[i].ins == cons[c].base[j].ins &&              //test insertion match
                    vbases[i].typ == cons[c].base[j].typ &&              //test type match
                    is_dgnrt_mtch(vbases[i].seq, cons[c].base[j].seq)) { //test seq match
                    
                    if (!uniq_match_tbl[j]) { //first match to current base constraint
                        fnd_vbMtch = 1;       //set found vbase match to TRUE
                        vb_mtchCnt++;         //increment vbase match counter
                        uniq_match_tbl[j]++;  //increment unique match tracker
                        
                    } else { //not the first match to the current constraint. throw error and abort
                        printf("test_bases: error - variant id contains duplicate base entry. aborting...\n");
                        abort();
                    }
                }
            }
        }
        
        if (vb_mtchCnt != vbase_cnt) { //test for complete vbase match
            match[c] = 0;              //if not complete match, set match[c] to FALSE
        }
        
        /**** test if vbases of input variant match pair constraints ****/
        for (i = 0; i < cons[c].pcnt; i++) {
            
            mate_valid[0] = mate_valid[1] = 0; //zero mate validity flags
            tmp_vb = NULL;                     //set temp vb pointer to NULL
            base[0] = base[1] = 0;             //zero base array
            pair_type = PAIR_TYPE_INIT;        //initialize pair_type to PAIR_TYPE_INIT
            
            /**** test if pair constraint bases match base constraints/reference sequence ****/
            for (j = 0; j < 2; j++) { //for each member of the base pair
                
                if (cons[c].pair[i].bs[j].typ == VARIABLE) { //if pair constraint is a variable base
                    
                    //check that vbase is a match to a base constraint
                    if ((tmp_vb = search_4_vbase_match(vbases, &cons[c].pair[i].bs[j], cons[c].bcnt, DGNRT_VBASE_MATCH)) != NULL) {
                        mate_valid[j] = 1;     //set flag that mate is valid
                        base[j] = tmp_vb->seq; //set base to base params seq value
                        
                    } else {
                        match[c] = 0; //mate not valid, set match to zero
                    }
                    
                } else if (cons[c].pair[i].bs[j].typ == CONSTANT) { //if pair constraint is a constant base
                    
                    //check that vbase matches the reference sequence
                    if (cons[c].pair[i].bs[j].seq == bmap->lkp0[cons[c].pair[i].bs[j].pos][cons[c].pair[i].bs[j].ins]) {
                        
                        mate_valid[j] = 1;                   //set flag that mate is valid
                        base[j] = cons[c].pair[i].bs[j].seq; //set base to constant seq val
                        
                    } else {
                        match[c] = 0; //mate not valid, set match to zero
                    }
                    
                } else { //pair constraint base type must be VARIABLE or CONSTANT
                    printf("filter_values: error - invalid base type for pair constraint. aborting...");
                    abort();
                }
            }
            
            /**** test if bases are a match to the constrained pair type ****/
            if (mate_valid[0] && mate_valid[1]) {     //both mates match the current base constraints/ref seq
                
                pair_type = ispair(base[0], base[1]); //set integer value of the base pair using ispair function
                 
                switch (cons[c].pair[i].att) {
                    case ANY_PAIR: if (pair_type < 1)                  {match[c] = 0;} break;
                    case WC_PAIR:  if (pair_type < 2)                  {match[c] = 0;} break;
                    case STRONG:   if (pair_type < 3)                  {match[c] = 0;} break;
                    case WEAK:     if (pair_type < 1 || pair_type > 2) {match[c] = 0;} break;
                    case WEAK_AT:  if (pair_type != 2)                 {match[c] = 0;} break;
                    case WEAK_AU:  if (pair_type != 2)                 {match[c] = 0;} break;
                    case WEAK_GT:  if (pair_type != 1)                 {match[c] = 0;} break;
                    case WEAK_GU:  if (pair_type != 1)                 {match[c] = 0;} break;
                    case MISMATCH: if (pair_type)                      {match[c] = 0;} break;
                    case NO_CONSTRAINT: break;
                    default:
                        break;
                }
            } 
        }
        
        if (match[c]) {          //variant id matches the current constraint
            match_cnt++;         //increment match count
            if (match_cnt > 1) { //if variant matches multiple constraints, print warning
                printf("test_bases: WARNING - variant %s matches more than one constraint\n", p_id);
                //TODO: print list of multi-constraint matches?
            }
        }
    }
    
    return match_cnt;
}

/* get_sample_info: parse column headers to obtain sample names and validate formatting */
int get_sample_info(char * line, sample_values * smpl_vals, int * tot_vals)
{
    extern const char TECdsply_clmn_hdrs[4][32]; //column headers from TECdisplay_output_column_headers.c
    
    int i = 0;  //general purpose index
    int j = 0;  //general purpose index
    
    char * p_id = NULL;            //pointer to id column string
    char * p_nm = NULL;            //pointer to column name string
    char tmp_nm[MAX_NAME+1] = {0}; //array for storing sample_name from current column name
    char * tmp_typ = NULL;         //array for storing field type from current column name
    char * last_uscore = NULL;     //pointer to last '_' in column name, used for string splitting
    
    int proceed = 1;    //flag to proceed to the next iteration of the loop
    int sample = 0;     //index of smpl_vals array, number of samples
    int val_count = 1;  //counts the number of value fields for a given sample
    int field1 = 1;     //flag that current field is the first field of a sample and sample name should be stored
    int disordered = 0; //flag that field types for current sample are not in the expected order
    
    /* set start of sample value names string*/
    for (i = 0; line[i] != '\t' && line[i] && i < MAX_LINE; i++) {;} //iterate to first tab
    if (line[i] == '\t') {
        line[i] = '\0';
        p_id = &line[0];   //set pointer to column name of id column
        p_nm = &line[i+1]; //set values pointer to the index after the first tab
    } else {
        printf("filter_values: unexpected format for header line in merged values file. aborting...\n");
        abort();
    }
        
    //column name of id column contain the substring "id" //TODO: make this test look at end of substring?
    if (strstr(p_id, TECdsply_clmn_hdrs[TDSPLY_VID_HDR]) == NULL) {
        printf("get_sample_info: error - unexpected format for column header line. the name of the first column is expected to contain the substring 'id' aborting...\n");
        abort();
    }
    
    /* read column names in groups of VALUE_FIELDS (currently 3), which is the number of value
     fields expected for each sample. the sample name from each column name is assessed for
     concordance across the sample. the field type for each sampe name is assessed to determine
     that the field types have the correct name, and are in the correct order. */
    for (sample = 0, i = 0; proceed && sample < MAX_SAMPLE; sample++) {
        
        field1 = 1; //field1 flag signals to get sample name from first value field of a sample
                    //sample names from all subsequente fields are compared to the first sample
                    //name to confirm concordance
        
        /* process the colulmn names for each sample to identify sample name and field type. each
         sample has three value fields: bound (int), unbound (int), and frac_bound (double). this
         loop also determines the number of value fields that will be exected for each data line
         using the tot_vals variable */
        for (val_count = 1; proceed && val_count <= VALUE_FIELDS; val_count++, (*tot_vals)++) {
            
            //copy value field to tmp_val until the tab delimiter or terminal null is reached
            for (j = 0, last_uscore = NULL; p_nm[i] != '\t' && p_nm[i] && j < MAX_NAME; i++, j++) {
                tmp_nm[j] = p_nm[i];
                if (tmp_nm[j] == '_') {
                    last_uscore = &tmp_nm[j]; //set pointer to last '_' to split the column name into two
                                              //strings later. first string is sample name, second is the
                                              //field name
                }
            }
            tmp_nm[j] = '\0';
            
            if (p_nm[i] != '\t' && p_nm[i] != '\0') { //check
                printf("get_sample_info: error - column name in merged values file is too long. aborting...\n");
                abort();
            }
            
            if (last_uscore != NULL) {     //column name must contain '_', otherwise formatting is invalid
                *last_uscore = '\0';       //set last '_' to null to split string
                tmp_typ = &last_uscore[1]; //set pointer to the field type string
                
                if (field1) {
                    //the value field is the first field of a sample. copy the string to the sample_values
                    //struct so that the two subsequent field names can be validated against the first
                    strcpy(smpl_vals[sample].nm, tmp_nm); //copy sample name to sample_values struct
                    field1 = 0;                           //turn off field 1 flag
                } else {
                    //the value field is not the first field of a sample. compare the sample name of the
                    //current field to the sample name that was obtained from the first column name
                    if (strcmp(tmp_nm, smpl_vals[sample].nm)) {
                        printf("get_sample_info: error - column names for sample %s are discordant. aborting...\n", smpl_vals[sample].nm);
                        abort();
                    }
                }
                
                disordered = 0; //zero disordered flag for current analysis
                
                //test the field type string from the current column against the field type
                //string that is expected to be in this column location. 'bound' comes first,
                //then 'unbound', then 'fracBound'.
                switch (val_count) {
                    case TDSPLY_BND_HDR: disordered = strcmp(tmp_typ, TECdsply_clmn_hdrs[TDSPLY_BND_HDR]); break;
                    case TDSPLY_UNB_HDR: disordered = strcmp(tmp_typ, TECdsply_clmn_hdrs[TDSPLY_UNB_HDR]); break;
                    case TDSPLY_FRC_HDR: disordered = strcmp(tmp_typ, TECdsply_clmn_hdrs[TDSPLY_FRC_HDR]); break;
                    default:
                        //note: this should be unreachable because val_count
                        //should never exceecd the VALUE_FIELDS constant
                        printf("get_sample_info: error - val_count exceeded the VALUE_FIELDS constant. this should never happen! aborting...\n");
                        abort();
                        break;
                }
                
                if (disordered) {
                    //if the disordered flag is on, the column names were read in an unexpected order
                    printf("get_sample_info: error - column names for sample %s are disordered. aborting...\n", smpl_vals[sample].nm);
                    abort();
                }
                
            } else {
                //all column names should contain an underscore that separates the sample name from the field type string
                printf("get_sample_info: error - unexpected format for column header name. column header string did not contain an underscore that separates the sample name from the field type. aborting...\n");
                abort();
            }
            
            /* after processing the last column name, check that the number of observed sample
             columns is divisible by the expected number of VALUE_FIELDS per sample column. */
            if (!p_nm[i]) {
                proceed = 0; //set proceed to off to exit outermost loop
                if (((val_count) % VALUE_FIELDS) != 0) {
                    printf("filter_values: error - unanticipated number of value fields for sample. expected %d value fields. aborting...\n", VALUE_FIELDS);
                }
            } else {
                i++; //increment i to the index of the next samples first column name
            }
        }
    }
    
    return sample; //return number of samples that were processed
}

//
//  merge_reactivity_profiles.c
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#include <stdio.h>
#include <stdlib.h>

#include "../global/global_defs.h"
#include "../mkmtrx/mkmtrx_defs.h"

#include "../cotrans_preprocessor/cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"

#include "./merge_TECprobeVL_replicates_defs.h"
#include "./merge_TECprobeVL_replicates_structs.h"

#include "merge_reactivity_profiles.h"

/* merge_profiles: combine data from input reactivity profiles to generate merged reactivity profile files */
int merge_profiles(SM2_analysis_directory * an_dir, int dir_count, output_files * outfiles)
{
    //the op array contains instructions for how to combine
    //the input data in the merged reactivity profile
    
    //the following op codes are used
    //COPY 0  directly copy the data from the first input profile to the merged profile
    //RSUM 1  sum the input file read counts
    //RATE 2  calculate a new mutation rate from the summed read counts
    //BSUB 3  calculate a new bgrnd-subtracted mut rate from the new mut rates
    //MASK 4  mask the data in the merged profile (used to mask normalized reactivity and stddev fields)
    
    int op[PRFL_CLMNS] = {
        COPY,   // 0 Nucleotide
        COPY,   // 1 Sequence
        RSUM,   // 2 Modified_mutations
        RSUM,   // 3 Modified_read_depth
        RSUM,   // 4 Modified_effective_depth
        RATE,   // 5 Modified_rate
        RSUM,   // 6 Modified_off_target_mapped_depth
        RSUM,   // 7 Modified_low_mapq_mapped_depth
        RSUM,   // 8 Modified_mapped_depth
        RSUM,   // 9 Untreated_mutations
        RSUM,   //10 Untreated_read_depth
        RSUM,   //11 Untreated_effective_depth
        RATE,   //12 Untreated_rate
        RSUM,   //13 Untreated_off_target_mapped_depth
        RSUM,   //14 Untreated_low_mapq_mapped_depth
        RSUM,   //15 Untreated_mapped_depth
        RSUM,   //16 Denatured_mutations
        RSUM,   //17 Denatured_read_depth
        RSUM,   //18 Denatured_effective_depth
        RATE,   //19 Denatured_rate
        RSUM,   //20 Denatured_off_target_mapped_depth
        RSUM,   //21 Denatured_low_mapq_mapped_depth
        RSUM,   //22 Denatured_mapped_depth
        BSUB,   //23 Reactivity_profile
        MASK,   //24 Std_err
        MASK,   //25 HQ_profile
        MASK,   //26 HQ_stderr
        MASK,   //27 Norm_profile
        MASK,   //28 Norm_stderr
    };
    
    int tl = 0; //transcript length index
    int d = 0;  //directory index
    int i = 0;  //general purpose index
    int j = 0;  //general purpose index
    int k = 0;  //general purpose index
    
    int dir0_line = 0; //used to test success of get_line_local across all directories
    int EOF_rchd = 0;  //flag that the end of file was reached
    
    char line[MAX_RUNS][MAX_LINE] = {{'\0'}}; //array to store lines from each input file
    
    int col = 0;           //column index
    int hanging_tab = 0;   //flag that line contains a hangingtab
    
    double * nmrtr = NULL; //pointer to numerator for mut rate calculations
    double * denom = NULL; //pointer to denominator for mut rate calculations
    
    char in_vals[MAX_RUNS][PRFL_CLMNS+1][MAX_FIELD] = {{{'\0'}}}; //+1 handles hanging tabs
    char sequence[MAX_RUNS][END_MAX] = {{0}};  //arrays to store profile sequences
    int sq[MAX_RUNS] = {0};            //sequence indices
    double out_vals[PRFL_CLMNS] = {0}; //array to store calculated output values when needed
    
    char prev_seq[END_MAX] = {0};
    char * p_prev = {0};
        
    for (tl = an_dir[0].min_tl; tl <= an_dir[0].max_tl; tl++) {  //for every transcript length
        
        for (j = 0; j < MAX_RUNS; j++) {
            sequence[j][0] = '\0'; //zero first index of sequence storage
            sq[j] = 0;             //zero sequence array storage
        }
        
        for (i = 0, EOF_rchd = 0; !EOF_rchd; i++) {              //for every line of the profile file
            for (j = 0; j < MAX_RUNS; j++) {         //for each directory
                
                line[j][0] = '\0';                   //zero first index of line storage
                for (k = 0; k < PRFL_CLMNS+1; k++) { //for SM2 profile column
                    in_vals[j][k][0] = '\0';         //zero the input values strings
                }
            }
            
            for (j = 0; j < PRFL_CLMNS; j++) { //for each SM2 profile column
                out_vals[j] = 0;               //zero the output values array
            }
            
            for (d = 0; d < dir_count; d++) {  //for every input directory
                if (!d) {                      //if reading the first input directory
                    //record get_line_local success so that success for all other dirs can be checked below
                    dir0_line = get_line_local(line[d], an_dir[d].prf[tl]) > 0 ? 1 : 0;
                    
                //if reading non-first input directory, check get_line_local success
                } else if ((get_line_local(line[d], an_dir[d].prf[tl]) > 0 ? 1 : 0) != dir0_line) {
                    printf("merge_profiles: profiles for transcript length %d end asynchronously. aborting...\n", tl);
                    abort();
                }
            }
            
            if (!dir0_line) { //if get_line_local failed
                EOF_rchd = 1; //set EOF_rchd to true
            }
            
            if (!i) { //reading header line
                for (d = 0; d < dir_count; d++) {
                    validate_header(line[d]);
                }
                fprintf(outfiles->ofp[tl], "%s\n", line[0]); //print header line to output file
                
            //TODO: check conditional below - not sure if second test is needed
            } else if (i || (EOF_rchd && line[0][0])) { //internal line or EOF and last line contained info
                for (d = 0; d < dir_count; d++) {  //for each directory
                    hanging_tab = 0;               //zero the hanging tab flag
                    
                    //parse the line for value fields
                    for (j = 0, k = 0, col = 0; line[d][j] != '\r' && line[d][j] != '\0'; j++) {
                        if (line[d][j] == '\t') {        //reached the end of a value field
                            in_vals[d][col][k] = '\0';   //terminate the current input values string
                            col++;                       //increment the column index
                            k = 0;                       //reset k to 0
                            
                            if (line[d][j+1] == '\n' ||  //adjust column count when
                                line[d][j+1] == '\r' ||  //there is a hanging tab
                                line[d][j+1] == '\0' ) {
                                hanging_tab = 1;
                            }
                            
                        } else {
                            in_vals[d][col][k] = line[d][j]; //copy current char to current input vals string
                            k++;                             //increment k
                        }
                    }
                    
                    if (!hanging_tab) {             //if the line did not contain a hanging tab
                        in_vals[d][col][k] = '\0';  //terminate the input values string
                        col++;                      //increment the column count
                    }
                    
                    if (col != PRFL_CLMNS) { //test whether column count matches expected number
                        printf("merge_profiles: error - unexpected column number (%d). aborting...\n", col);
                        abort();
                    }
                }
            } else {
                printf("this should not be reachable. aborting...\n");
                abort();
            }
            
            //TODO: add profile validation to ensure matched profiles
            //TODO: this should involve comparing nucleotide sequence
            if (i) { //if reading a data line, perform column-specific operations
                for (col = 0; col < PRFL_CLMNS; col++) { //for each column
                    
                    if (col == SEQ_COL) {
                        for (d = 0; d < dir_count; d++) {
                            if (strlen(in_vals[d][col]) == 1) {
                                sequence[d][sq[d]++] = in_vals[d][col][0];
                                if (sq[d] == END_MAX) {
                                    printf("merge_profiles: error - sequence exceeded array bounds. aborting...\n");
                                    abort();
                                }
                                
                            } else {
                                printf("%s\n", in_vals[d][col]);
                                
                                printf("merge_profiles: error - expected a single character in Nucleotide column value field. aborting...\n");
                                abort();
                            }
                        }
                    }
                    
                    //if op=COPY, copy value from first file to merged out
                    if (op[col] == COPY) {
                        if (col == 0) { //only print tab if non-first col
                            fprintf(outfiles->ofp[tl], "%s", in_vals[0][col]);
                        } else {
                            fprintf(outfiles->ofp[tl], "\t%s", in_vals[0][col]);
                        }
                    
                    //if op=RSUM, sum read count from each profile
                    } else if (op[col] == RSUM) {
                        for (d = 0; d < dir_count; d++) {
                            out_vals[col] += strtof(in_vals[d][col], NULL);
                        }
                        fprintf(outfiles->ofp[tl], "\t%.0f", out_vals[col]);
                        
                    //if op=RATE calculate new mutation rate from new read counts
                    } else if (op[col] == RATE) {

                        switch (col) {
                            case MOD_RATE: //calculating modified mutation rate
                                nmrtr = &out_vals[MOD_MUT];
                                denom = &out_vals[MOD_EFF_DEP];
                                break;
                                
                            case UNT_RATE: //calculating untreated mutation rate
                                nmrtr = &out_vals[UNT_MUT];
                                denom = &out_vals[UNT_EFF_DEP];
                                break;
                            
                            case DEN_RATE: //calculating denatured mutation rate
                                nmrtr = &out_vals[DEN_MUT];
                                denom = &out_vals[DEN_EFF_DEP];
                                break;
                                
                            default:
                                break;
                        }
                        
                        out_vals[col] = (*nmrtr)/(*denom); //calculate new mutation rate
                        fprintf(outfiles->ofp[tl], "\t%.6f", out_vals[col]);
                        
                    //if op=BSUB, calculate new background-subtracted mutation rate
                    } else if (op[col] == BSUB) {
                        
                        //if no denatured reads are detected, perform MOD_RATE - UNT_RATE
                        if (((int)out_vals[DEN_READ_DEP]) == 0) {
                            out_vals[col] = out_vals[MOD_RATE] - out_vals[UNT_RATE];
                        
                        //if denatured reads are detected, perform (MOD_RATE-UNT_RATE)/DEN_RATE
                        } else {
                            out_vals[col] = (out_vals[MOD_RATE] - out_vals[UNT_RATE])/out_vals[DEN_RATE];
                        }
                        fprintf(outfiles->ofp[tl], "\t%.6f", out_vals[col]);
                        
                    //if op=MASK, mask output data
                    } else if (op[col] == MASK) {
                        //TODO: what about cases where SC1 leader is not present, should detect?
                        if (i <= 16) { //if in SC1 leader, mask with nan
                            fprintf(outfiles->ofp[tl], "\tnan");
                        } else {       //if beyond SC1 leader, mask with 0.000000
                            fprintf(outfiles->ofp[tl], "\t0.000000");
                        }
                    }
                }
                
                if (!EOF_rchd) { //print newline if EOF was not yet reached
                    fprintf(outfiles->ofp[tl], "\n");
                }
            }
        }
        
        for (d = 0; d < dir_count; d++) { //for each directory
            sequence[d][sq[d]++] = '\0';  //terminate the sequence string
            if (d && strcmp(sequence[d], sequence[0])) {
                printf("merge_profiles: discordant sequences detected in input directory. aborting...\n");
                abort();
            }
        }
        
        //check that previous sequence is a substring of the current sequence
        if ((p_prev = strstr(sequence[0], prev_seq)) == NULL) {
            printf("merge_profiles: error - sequence of transcript %d is not a substring of transcript %d. aborting...\n", tl-1, tl);
            abort();
        
        //if processing non-min transcript length, check that the length of the
        //current sequence string is one nt less than the length of the previous
        //seuence string
        } else if (tl > an_dir[0].min_tl ) {
            
            if ((uint64_t)(p_prev) != (uint64_t)(&sequence[0][0])) {
                printf("merge_profiles: error - sequence of transcripts %d and %d do not start at the same position. aborting...\n", tl-1, tl);
                abort();
                
            } else if (strlen(prev_seq) != (strlen(sequence[0])-1)) {
                printf("merge_profiles: error - sequence of transcript %d is not 1 nt shorter than sequence of transcript %d. aborting...\n", tl-1, tl);
                abort();
            }
        }
        strcpy(prev_seq, sequence[0]);
    }
    return 1;
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
            line[++i] = '\0'; //append terminating null character
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

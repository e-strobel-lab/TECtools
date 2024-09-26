//
//  mk_rdat.c
//  
//
//  Created by Eric Strobel on 2/9/23.
//

#include "mk_rdat.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../utils/io_management.h"

#include "./cotrans_mtrx.h"
#include "./mkmtrx_defs.h"
#include "./mkmtrx_structs.h"

#include "../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"

/* mk_rdat: manages rdata config parsing/checking and
 construction of output rdat file*/
int mk_rdat(FILE * fp_config, cotrans_matrix * mtrx, int mode, int reps)
{
    rdat_metadata rdat_meta = {{0}}; //struct for rdat metadata
    rdat_meta.rtype = -1;            //intialize rtype to -1 because possible values are >= 0
    rdat_meta.smooth = -1;           //initialize smooth to -1 because true/false is 1/0
    
    parse_rdat_config(fp_config, &rdat_meta); //parse rdat config file
    check_rdat_config(&rdat_meta, mtrx);      //check that config contains required info
    print_rdat(&rdat_meta, mtrx, mode, reps); //print rdat file
    return 1;
}

/* parde_rdat_config: read rdat config file and store values in rdat_metadata struct */
int parse_rdat_config(FILE * ifp, rdat_metadata * rdat_meta)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    char line[MAX_LINE] = {0};       //array for storing config file line
    char setting[MAX_LINE] = {0};    //array for storing config file setting name
    char val[MAX_LINE] = {0};        //array for storing config file setting value
    
    char * setting_start = {NULL};   //pointer to start of config file setting in input line
        
    //in each iteration, get a line from the config file, check if it contains
    //a setting, and set the corresponding value in the config_MLT struct
    while (get_line(line, ifp)) {
        
        //reset setting and val variables
        setting[0] = '\0';
        val[0] = '\0';
        
        //all setting lines begin with an alphabetic characters
        if (isalpha(line[0])) {
            
            //copy setting name to setting array
            setting_start = &line[0]; //set pointer to start of setting name
            for (i = 0; setting_start[i] != '=' && setting_start[i] && i < MAX_LINE; i++) {
                setting[i] = setting_start[i];
            }
            setting[i] = '\0';
            
            if (setting_start[i] != '=') { //check that loop ended on '=' character
                printf("parse_rdat_config: error unexpected setting name format for setting %s. aborting... \n", setting);
                abort();
            }
            
            //copy setting value to val array
            for (i++, j = 0; setting_start[i] && i < MAX_LINE; i++, j++) {
                val[j] = setting_start[i];
            }
            val[j] = '\0';
            
            if (setting_start[i]) { //check that loop ended on null character
                printf("parse_rdat_config: error unexpected setting value format. aborting... \n");
                abort();
            }
            
            //check for setting value errors
            //required settings have a default value of "FILL_IN"
            //optional settings have a default value of "NULL"
            //and no setting should have an empty value
            if (!strcmp(val, "FILL_IN")) { //check if required setting value has been omitted
                printf("parse_rdat_config: error - missing required value for %s setting. aborting...\n", setting);
                abort();
            } else if (strstr(val, "FILL_IN") != NULL) {   //check if required setting value has been omitted
                printf("parse_rdat_config: error - missing required value for %s %s setting. aborting...\n", setting, val);
                abort();
            } else if (val[0] == '\0') {   //check if setting value is empty in config file
                printf("parse_rdat_config: error - value for %s setting is empty. optional settings without a value must be set to NULL. aborting...\n", setting);
                abort();
            } else if (!strcmp(val, "NULL")) {
                //nothing to do here, just forcing NULL values away from the if-else statement below
            } else {
                
                //set output rdat file name
                if (!strcmp(setting, "OUTPUT_FILE_NAME")) {
                    strcpy(rdat_meta->fn, val);
                    
                //set sample name
                } else if (!strcmp(setting, "NAME")) {
                    strcpy(rdat_meta->nm, val);
                    
                //TODO: NOT CURRENTLY USED, structure is set to a string of '.'
                } else if (!strcmp(setting, "STRUCTURE")) {
                    strcpy(rdat_meta->str, val);
                    
                //set offset
                } else if (!strcmp(setting, "OFFSET")) {
                    for (i = 0; val[i]; i++) {  //check that offset value is a string of digits
                        if (!isdigit(val[i])) {
                            printf("parse_rdat_config: error - offset value must be a string of digits. aborting...\n");
                            abort();
                        }
                    }
                    rdat_meta->offset = atoi(val);
                
                //set reactivity type
                } else if (!strcmp(setting, "REACTIVITY_TYPE")) {
                    if (!strcmp(val, "RAW")) {
                        rdat_meta->rtype = RAW;
                    } else if (!strcmp(val, "NORMALIZED")) {
                        rdat_meta->rtype = NORMALIZED;;
                    } else {
                        printf("parse_rdat_config: error - unrecognized reactivity type. reactiviity type must be 'RAW' or 'NORMALIZED'. aborting...");
                        abort();
                    }
                    
                //set smooth
                } else if (!strcmp(setting, "SMOOTHED")) {
                    set_TF_value(val, setting, &(rdat_meta->smooth));
                    
                //set annotation
                } else if (!strcmp(setting, "ANNOTATION")) {
                    if (strstr(val, "experimentType:") != NULL) {
                        strcpy(rdat_meta->xtype, val);
                        
                    } else if (strstr(val, "experiment:") != NULL) {
                        strcpy(rdat_meta->xmeth, val);
                        
                    } else if (strstr(val, "modifier:") != NULL) {
                        strcpy(rdat_meta->probe, val);
                        
                    } else if (strstr(val, "chemical:Tris-HCl:") != NULL) {
                        strcpy(rdat_meta->tris, val);
                        
                    } else if (strstr(val, "chemical:KCl:") != NULL) {
                        strcpy(rdat_meta->KCl, val);
                        
                    } else if (strstr(val, "chemical:EDTA:") != NULL) {
                        strcpy(rdat_meta->EDTA, val);
                        
                    } else if (strstr(val, "chemical:DTT:") != NULL) {
                        strcpy(rdat_meta->DTT, val);
                        
                    } else if (strstr(val, "chemical:MgCl2:") != NULL) {
                        strcpy(rdat_meta->MgCl2, val);
                        
                    } else if (strstr(val, "chemical:NTPs:") != NULL) {
                        strcpy(rdat_meta->NTPs, val);
                        
                    } else if (strstr(val, "chemical:BSA:") != NULL) {
                        strcpy(rdat_meta->BSA, val);
                        
                    } else if (strstr(val, "chemical:") != NULL) {
                        if (rdat_meta->oc_cnt < MAX_OTHER_CHEM) {
                            strcpy(rdat_meta->other_chem[rdat_meta->oc_cnt++], val);
                        } else {
                            printf("parse_rdat_config: error - number of other chemicals exceeds the maximum (%d). aborting...\n", MAX_OTHER_CHEM);
                            abort();
                        }
                        
                    } else if (strstr(val, "temperature:") != NULL) {
                        strcpy(rdat_meta->temp, val);
                        
                    } else {
                        printf("parse_rdat_config: error - unrecognized annotation %s. aborting...\n", val);
                        abort();
                    }
                    
                //set user supplied comments
                } else if (strstr(setting, "COMMENT")) {
                    if (rdat_meta->cmnt_cnt < MAX_COMMENTS) {
                        strcpy(rdat_meta->comment[rdat_meta->cmnt_cnt++], val);
                    } else {
                        printf("parse_rdat_config: error - number of user-supplied comments exceeds the maximum (%d). aborting...\n", MAX_COMMENTS);
                        abort();
                    }
                    
                    
                //throw error for unexpected setting name
                } else {
                    printf("parse_rdat_config: error - unrecognized setting (%s). aborting...\n", setting);
                    abort();
                }
            }
        } else if (line[0] != '#') {  //# can be used for comments in rdat_config
            printf("parse_rdat_config: error - unexpected line initator (%c). lines must start with # (comment) or an alphabetic character (setting). aborting...\n", line[0]);
            abort();
        }
    }    
    return 1;
}

/* check_rdat_config: check that required rdat metadata has been supplied
 and is consistent with the shapemapper2 output */
int check_rdat_config(rdat_metadata * rdat_meta, cotrans_matrix * mtrx)
{
    int i = 0; //general purpose index
     
    FILE * out_fp = NULL;              //parsed config file pointer
    char out_nm[MAX_LINE+23+1] = {0};  //parsed file name
    
    //generate output parsed config file
    if (!rdat_meta->fn[0]) {
        printf("check_rdat_config: error - output filename missing. aborting...\n");
        abort();
    } else {
        sprintf(out_nm, "parsed_rdat_config_%s.txt", rdat_meta->fn);
        if ((out_fp = fopen(out_nm, "w")) == NULL) {
            printf("check_rdat_config: ERROR - could not generate parsed_rdat_config file. Aborting program...\n");
            abort();
        }
    }
     
    //check name variable
    printf("parsed rdat config:\n");
    fprintf(out_fp, "parsed rdat config:\n");
    if (!rdat_meta->nm[0]) {
        printf("check_rdat_config: error - name is not set. aborting\n");
        abort();
    } else {
        printf("sample name\t%s\n", rdat_meta->nm);
        fprintf(out_fp, "sample name\t%s\n", rdat_meta->nm);
    }
        
    //check structure variable
    //TODO: NOT CURRENTLY USING STRUCTURE VALUE
    /*if (!rdat_meta->str[0]) {
        printf("check_rdat_config: error - structure is not set. aborting\n");
        abort();
    } else {
        printf("structure\t%s\n", rdat_meta->str);
        fprintf(out_fp, "structure\t%s\n", rdat_meta->str);
    }*/
    
    //TODO: if using structure variable, will check if it matches mtrx->seq-1 len
    //check that length of sequence and structure match
    /*if (strlen(rdat_meta->seq) != strlen(rdat_meta->str)) {
        printf("check_rdat_config: error - length of sequence and structure do not match. aborting...\n");
        fprintf(out_fp, "check_rdat_config: error - length of sequence and structure do not match. aborting...\n");
        abort();
    }*/
    
    //check that offset is not less than zero or >= sequence length
    //less than zero should be impossible, since input must be a string of digits, check here anyway
    //strlen(mtrx->sq)-1) is used because index zero of mtrx->sq is '>'
    if (rdat_meta->offset < 0 || rdat_meta->offset >= (strlen(mtrx->sq)-1)) {
        printf("check_rdat_config: error - invalid value for offset (%d). aborting\n", rdat_meta->offset);
        abort();
    } else {
        printf("offset\t\t%d\n", rdat_meta->offset);
        fprintf(out_fp, "offset\t\t%d\n", rdat_meta->offset);
    }
    
    //check that reactivity type was set
    if (rdat_meta->rtype < 0 || rdat_meta->rtype > 1) {
        printf("check_rdat_config: error - invalid value for reactivity type. aborting\n");
        abort();
    } else {
        if (rdat_meta->rtype == RAW) {
            printf("reactivity\tRAW\n");
            fprintf(out_fp, "reactivity\tRAW\n");
        } else if (rdat_meta->rtype == NORMALIZED) {
            printf("reactivity\tNORMALIZED\n");
            fprintf(out_fp, "reactivity\tNORMALIZED\n");
        }
    }
    
    //check that smoothing flag was set to OFF (0) or ON (1)
    if (rdat_meta->smooth != 0 && rdat_meta->smooth != 1) {
        printf("check_rdat_config: error - invalid value for smoothed. aborting\n");
        abort();
    } else {
        if (rdat_meta->smooth) {
            printf("smoothed\tTRUE\n");
            fprintf(out_fp, "smoothed\tTRUE\n");
        } else {
            printf("smoothed\tFALSE\n");
            fprintf(out_fp, "smoothed\tFALSE\n");
        }
    }
    
    //check experimentType variable
    if (!rdat_meta->xtype[0]) {
        printf("check_rdat_config: error - experimentType is not set. aborting\n");
        abort();
    } else {
        printf("experimentType\t%s\n", rdat_meta->xtype);
        fprintf(out_fp, "experimentType\t%s\n", rdat_meta->xtype);
    }
    
    //check experiment variable
    if (!rdat_meta->xmeth[0]) {
        printf("check_rdat_config: error - experiment is not set. aborting\n");
        abort();
    } else {
        printf("experiment\t%s\n", rdat_meta->xmeth);
        fprintf(out_fp, "experiment\t%s\n", rdat_meta->xmeth);
    }
    
    //check modifier variable
    if (!rdat_meta->probe[0]) {
        printf("check_rdat_config: error - modifier is not set. aborting\n");
        abort();
    } else {
        printf("modifier\t%s\n", rdat_meta->probe);
        fprintf(out_fp, "modifier\t%s\n", rdat_meta->probe);
    }
    
    //check tris variable
    if (!rdat_meta->tris[0]) {
        printf("check_rdat_config: warning - tris is not set.\n");
    } else {
        printf("\ntris\t%s\n", rdat_meta->tris);
        fprintf(out_fp, "\ntris\t%s\n", rdat_meta->tris);
    }
    
    //check KCl variable
    if (!rdat_meta->KCl[0]) {
        printf("check_rdat_config: warning - KCl is not set.\n");
    } else {
        printf("KCl\t%s\n", rdat_meta->KCl);
        fprintf(out_fp, "KCl\t%s\n", rdat_meta->KCl);
    }
    
    //check EDTA variable
    if (!rdat_meta->EDTA[0]) {
        printf("check_rdat_config: WARNING - EDTA is not set.\n");
    } else {
        printf("EDTA\t%s\n", rdat_meta->EDTA);
        fprintf(out_fp, "EDTA\t%s\n", rdat_meta->EDTA);
    }
    
    //check DTT variable
    if (!rdat_meta->DTT[0]) {
        printf("check_rdat_config: warning - DTT is not set.\n");
    } else {
        printf("DTT\t%s\n", rdat_meta->DTT);
        fprintf(out_fp, "DTT\t%s\n", rdat_meta->DTT);
    }
    
    //check MgCl2 variable
    if (!rdat_meta->MgCl2[0]) {
        printf("check_rdat_config: warning - MgCl2 is not set.\n");
    } else {
        printf("MgCl2\t%s\n", rdat_meta->MgCl2);
        fprintf(out_fp, "MgCl2\t%s\n", rdat_meta->MgCl2);
    }
    
    //check NTPs variable
    if (!rdat_meta->NTPs[0]) {
        printf("check_rdat_config: warning - NTPs is not set.\n");
    } else {
        printf("NTPs\t%s\n", rdat_meta->NTPs);
        fprintf(out_fp, "NTPs\t%s\n", rdat_meta->NTPs);
    }
    
    //check BSA variable
    if (!rdat_meta->BSA[0]) {
        printf("check_rdat_config: error - BSA is not set.\n");
    } else {
        printf("BSA\t%s\n", rdat_meta->BSA);
        fprintf(out_fp, "BSA\t%s\n", rdat_meta->BSA);
    }
    
    //print other chemical settings
    for (i = 0; i < rdat_meta->oc_cnt; i++) {
        printf("other\t%s\n", rdat_meta->other_chem[i]);
        fprintf(out_fp, "other\t%s\n", rdat_meta->other_chem[i]);
    }
    
    //check temperature variable
    if (!rdat_meta->temp[0]) {
        printf("check_rdat_config: error - temperature is not set. aborting\n");
        abort();
    } else {
        printf("temp\t%s\n", rdat_meta->temp);
        fprintf(out_fp, "temp\t%s\n", rdat_meta->temp);
    }
    
    //print user supplied comments
    if (rdat_meta->cmnt_cnt) {
        printf("\nuser-supplied comments:\n");
        fprintf(out_fp, "\nuser-supplied comments:\n");
        for (i = 0; i < rdat_meta->cmnt_cnt; i++) {
            printf("comment\t%s\n", rdat_meta->comment[i]);
            fprintf(out_fp, "comment\t%s\n", rdat_meta->comment[i]);
        }
    }

    //close output file
    if ((fclose(out_fp)) == EOF) {
        printf("check_rdat_config: ERROR - error occurred when closing parsed_rdat_config  file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

/* print_rdat: print rdat file */
int print_rdat(rdat_metadata * rdat_meta, cotrans_matrix * mtrx, int mode, int reps)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    FILE * out_fp = NULL;            //output file pointer
    char out_nm[MAX_LINE+5+1] = {0}; //output file name
    
    //generate output parsed config file
    sprintf(out_nm, "%s.rdat", rdat_meta->fn);
    if ((out_fp = fopen(out_nm, "w")) == NULL) {
        printf("print_rdat: ERROR - could not generate rdat file. Aborting program...\n");
        abort();
    }
     
    fprintf(out_fp, "RDAT_VERSION\t0.4\n");         //print rdat version line
    fprintf(out_fp, "NAME\t%s\n", rdat_meta->nm);   //print name line
    
    //i starts at 1 because mtrx->sq[0] is '>'
    fprintf(out_fp, "SEQUENCE\t");                  //print sequence line
    for (i = 1; mtrx->sq[i]; i++) {                 //using sequence obtained from
        fprintf(out_fp, "%c", mtrx->sq[i]);         //shapemapper 2 output files
    }
    fprintf(out_fp, "\n");
    
    //i starts at 1 because mtrx->sq[0] is '>'
    fprintf(out_fp, "STRUCTURE\t");                 //print structure line
    for (i = 1; mtrx->sq[i]; i++) {                 //TODO: currently this prints a string of '.'
        fprintf(out_fp, ".");                       //that matches sequence length
    }
    fprintf(out_fp, "\n");
    
    fprintf(out_fp, "OFFSET\t%d\n", rdat_meta->offset); //print offset value
    fprintf(out_fp, "SEQPOS");                          //print seqpos line
    for (i = 1; mtrx->sq[i]; i++) { //i starts at 1 because mtrx->sq[0] is '>'
        fprintf(out_fp, "\t%c%d", mtrx->sq[i], i-rdat_meta->offset);
    }
    fprintf(out_fp, "\n\n");
    
    //printf annotations
    fprintf(out_fp, "ANNOTATION");
    fprintf(out_fp, "\t%s", rdat_meta->xtype);
    fprintf(out_fp, "\t%s", rdat_meta->xmeth);
    fprintf(out_fp, "\t%s", rdat_meta->probe);
    if (rdat_meta->tris[0]) {
        fprintf(out_fp, "\t%s", rdat_meta->tris);
    }
    if (rdat_meta->KCl[0]) {
        fprintf(out_fp, "\t%s", rdat_meta->KCl);
    }
    if (rdat_meta->EDTA[0]) {
        fprintf(out_fp, "\t%s", rdat_meta->EDTA);
    }
    if (rdat_meta->DTT[0]) {
        fprintf(out_fp, "\t%s", rdat_meta->DTT);
    }
    if (rdat_meta->MgCl2[0]) {
        fprintf(out_fp, "\t%s", rdat_meta->MgCl2);
    }
    if (rdat_meta->NTPs[0]) {
        fprintf(out_fp, "\t%s", rdat_meta->NTPs);
    }
    if (rdat_meta->BSA[0]) {
        fprintf(out_fp, "\t%s", rdat_meta->BSA);
    }
    for (i = 0; i < rdat_meta->oc_cnt; i++) {
        fprintf(out_fp, "\t%s", rdat_meta->other_chem[i]);
    }
    fprintf(out_fp, "\t%s", rdat_meta->temp);
    fprintf(out_fp, "\n\n");
    
    //print comments
    fprintf(out_fp, "COMMENT\tSequencing reads were demultiplexed by transcript length and channel using the TECtools cotrans_preprocessor script\n");
    if (rdat_meta->smooth) {
        fprintf(out_fp, "COMMENT\tNeighboring transcript smoothing was applied\n");
    }
    
    fprintf(out_fp, "COMMENT\tMutation mapping and reactivity calculation was performed using shapemapper2 (Busan and Weeks,  DOI: 10.1261/rna.061945.117)\n");
    
    if (rdat_meta->rtype == RAW) {
        fprintf(out_fp, "COMMENT\tReactivity values are raw backround-subtracted mutations rates\n");
    } else if (rdat_meta->rtype == NORMALIZED) {
        fprintf(out_fp, "COMMENT\tReactivity values were normalized as described in section 3.2.1 of Low and Weeks, 2010, DOI: 10.1016/j.ymeth.2010.06.007\n");
    }
    
    fprintf(out_fp,"COMMENT\tThe leading SC1 hairpin sequence (5-ATGGCCTTCGGGCCAA), which is masked by a primer, was trimmed\n");
    if (mode == SINGLE) {
        fprintf(out_fp, "COMMENT\tTrailing Illumina adapter sequence (GATCGTCGGACTGTAGAACTCTGAAC-3), which is masked by a primer, was trimmed\n");
    }
    if (mode == MULTI) {
        fprintf(out_fp, "COMMENT\tTranscripts from %d to %d and after %d were not enriched by template DNA strand biotin-streptavidin roadblocks and may be poor quality\n", mtrx->last_iStl+1, mtrx->frst_tStl-1, mtrx->last_tStl);
    }
    
    //print user supplied comments
    if (rdat_meta->cmnt_cnt) {
        for (i = 0; i < rdat_meta->cmnt_cnt; i++) {
            fprintf(out_fp, "COMMENT\t%s\n", rdat_meta->comment[i]);
        }
    }
    
    fprintf(out_fp, "\n");
    
    //print data annotations
    int a_indx = 0;  //annotation index
    int col_cnt = 0; //number of columns in current row
    
    int cRep = 0;    //current replicate index
    
    for (i = mtrx->tl[MIN], a_indx = 1; i <= mtrx->tl[MAX]; i++) {
        
        //print REACTIVITY ANNOTATION
        col_cnt = atoi(mtrx->vals[i][0]); //set column count, which is stored at index zero of the row
        
        //print a data annotation line for each replicate
        for (cRep = 0; cRep < reps; cRep++) {
            fprintf(out_fp, "DATA_ANNOTATION:%d\tsequence:", a_indx++);
            for (j = mtrx[cRep].nt[MIN]; j < mtrx[cRep].nt[MIN]+col_cnt; j++) {
                fprintf(out_fp, "%c", mtrx[cRep].sq[j]);
            }
            fprintf(out_fp, "\tdatatype:REACTIVITY\t%s\tID:Length%d", rdat_meta->probe, col_cnt);
            
            if (mode == LM_RDAT) {
                fprintf(out_fp, "_Sample%d", i);
            }
            
            if (reps > 1) {
                fprintf(out_fp, "_Rep%d", cRep+1);
            }
            
            fprintf(out_fp, "\n");
        }
        
        
    }
    
    for (i = mtrx->tl[MIN], a_indx = 1; i <= mtrx->tl[MAX]; i++) {
        
        //print REACTIVITY data
        col_cnt = atoi(mtrx->vals[i][0]); //set column count, which is stored at index zero of the row
        
        //print a data line for each replicate
        for (cRep = 0; cRep < reps; cRep++) {
            fprintf(out_fp, "DATA:%d", a_indx++);
            for (j = mtrx[cRep].nt[MIN]; j < mtrx[cRep].nt[MIN]+col_cnt; j++) {
                fprintf(out_fp, "\t%s", mtrx[cRep].vals[i][j]);
            }
            fprintf(out_fp, "\n");
        }
    }
    
    //close output file
    if ((fclose(out_fp)) == EOF) {
        printf("check_rdat_config: ERROR - error occurred when closing rdat file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

//
//  parse_VL_sample_name.c
//  
//
//  Created by Eric Strobel on 2/5/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../utils/debug.h"
#include "../../utils/io_management.h"

#include "../../cotrans_preprocessor/run_script_gen/UNV/config_struct.h"
#include "../../cotrans_preprocessor/run_script_gen/UNV/mk_run_nm.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"

#include "parse_VL_sample_name.h"

/* parse_VL_sample_name: parse sample name for attributes that were specified in the TECprobe analysis config */
void parse_VL_sample_name(char * ipt_nm, tprobe_configuration * cfg)
{
    extern int debug;
    
    int i = 0; //general purpose index
    int u = 0; //counts underscores during runID parsing
    
    char snm[MAX_LINE] = {0};            //copy of the sample name to be parsed
    parsed_sample_name prsd_sn = {NULL}; //storage for parsed sample name sub-strings
    
    char CoTxn[8] = {"_CoTxn_"}; //cotranscriptional folding type indicator, used as parsing anchor
    char Equil[8] = {"_Equil_"}; //equilibrium folding type indicator, used as parsing anchor
    char * scnd_srch = NULL;     //pointer used to perform a search for a second anchor match
    char * fld_typ = NULL;       //pointer to fold folding type indicator (CoTxn or Equil variables)
    char * othr_typ = NULL;      //pointer to the folding type indicator that was not identified as an anchor
    
    char * p_SM_flag = NULL;     //pointer to the smoothing flag
    
    char *tmp_ptr = NULL;        //temporary pointer
    
    int parse_complete = 0;      //flag that entire string was parsed
    
    char str1[2] = {"1"};        //string for setting runCT pointer when in non-concatenated data sets
    
    int found_conc_start = 0;    //flag that start of concentration string was found
    
    int ret = 0; //variable for storing snprintf return value
    
    //copy sample name to snm array for parsing into sub-strings
    ret = snprintf(snm, MAX_LINE, "%s", ipt_nm);
    if (ret >= MAX_LINE || ret < 0) {
        printf("parse_sample_name: error when storing copy of sample name. aborting...\n");
        abort();
    }
    
    /*** parse sample name to identify substring attributes ***/
        
    //search for the '_CoTxn' and '_Equil' anchor sequences
    if ((prsd_sn.anchr = strstr(snm, CoTxn)) != NULL) {        //if string contains '_CoTxn'
        fld_typ = &CoTxn[0];                                   //set fold type pointer to '_CoTxn'
        othr_typ = &Equil[0];                                  //and other type pointer to '_Equil'
    } else if ((prsd_sn.anchr = strstr(snm, Equil)) != NULL) { //if string contains '_Equil'
        fld_typ = &Equil[0];                                   //set fold type pointer to '_Equil'
        othr_typ = &CoTxn[0];                                  //and other type pointer to '_CoTxn'
    } else { //throw error and abort if no anchor was detected
        printf("parse_sample_name: error - could not detect folding type. aborting...\n");
        abort();
    }
    
    //search for a second anchor sequence. if a second anchor sequence is found,
    //throw error and abort. otherwise, proceed with sample name parsing
    if ((scnd_srch = strstr(&prsd_sn.anchr[strlen(fld_typ)], fld_typ)) != NULL) {
        printf("parse_sample_name: error - sample name contains two instances of the fold type string (%s). aborting...\n", fld_typ);
        
    } else if ((scnd_srch = strstr(&prsd_sn.anchr[strlen(fld_typ)], othr_typ)) != NULL) {
        printf("parse_sample_name: error - sample name contains both %s and %s fold type strings. aborting...\n", fld_typ, othr_typ);
        
    } else { //proceed with sample name parsing
            
        //set folding type substring
        prsd_sn.fld = &prsd_sn.anchr[1];         //set fld substring pointer to start of fold type string
        prsd_sn.anchr[strlen(fld_typ)-1] = '\0'; //truncate the fold substring by replacing uscore with '\0'
        
        //set RNA name substring
        prsd_sn.rna = &snm[0]; //set rna name pointer to start of snm
        if (((uint64_t)(prsd_sn.rna)) < ((uint64_t)(prsd_sn.anchr))) { //if chars precede the anchor
            prsd_sn.anchr[0] = '\0'; //set the uscore that precedes fold type to '\0' to terminate rna name
            
        } else { //if chars do not precede the anchor, throw error and abort
            printf("parse_sample_name: error - expected characters to precede the fold type string in sample name\n");
            abort();
        }
        
        //set probe type substring
        prsd_sn.prb = &prsd_sn.anchr[strlen(fld_typ)]; //set prb substring pointer to char after anchor
        for (i = 0; prsd_sn.prb[i] && prsd_sn.prb[i] != '_'; i++) {;} //find the index of uscore after prb type
        if (prsd_sn.prb[i]) {                //if the uscore was found
            prsd_sn.prb[i] = '\0';           //replace it with a '\0' to terminate the probe type substring
            prsd_sn.rst = &prsd_sn.prb[i+1]; //point the 'rest' pointer to the char after the uscore
            
        } else { //if no uscore was detected, throw error and abort
            printf("parse_sample_name: error - expected additional information beyond probe type string. aborting...");
            abort();
        }
        
        //set smoothing flag/substring
        if ((p_SM_flag = strstr(prsd_sn.rst, "_SM"))) { //search for smoothing flag at end of sample name
            if (!p_SM_flag[strlen(p_SM_flag)]) {        //confirm smoothing flag is at end of sample name
                prsd_sn.SM = &p_SM_flag[1];             //set SM flag substring pointer to SM flag
                p_SM_flag[0] = '\0';                    //change preceding uscore to '\0' to truncate SM flag
            }
        }
        
        //set ligand name/conc substrings
        //perform search for concentration unit substrings
        if ((prsd_sn.conc = strstr(prsd_sn.rst, "mM_")) != NULL) {
            ;
        } else if ((prsd_sn.conc = strstr(prsd_sn.rst, "uM_")) != NULL) {
            ;
        } else if ((prsd_sn.conc = strstr(prsd_sn.rst, "nM_")) != NULL) {
            ;
        } else if ((prsd_sn.conc = strstr(prsd_sn.rst, "pM_")) != NULL) {
            ;
        } else if ((prsd_sn.conc = strstr(prsd_sn.rst, "fM_")) != NULL) {
            ;
        } else if ((prsd_sn.conc = strstr(prsd_sn.rst, "aM_")) != NULL) {
            ;
        }
        
        //if concentration unit substring was found and is preceeded by characters
        if (prsd_sn.conc != NULL && (((uint64_t)(prsd_sn.conc)) - ((uint64_t)(prsd_sn.rst))) > 0) {
            
            for (i = -1, found_conc_start = 0; ((uint64_t)(&prsd_sn.conc[i])) != ((uint64_t)(prsd_sn.rst)) && !found_conc_start; i--) {
                
                if (!isdigit(prsd_sn.conc[i]) && prsd_sn.conc[i] != '_') {
                    printf("parse_VL_sample_name: error - found non-digit character %c in ligand concentration string. aborting...\n", prsd_sn.conc[i]);
                    abort();
                } else if (prsd_sn.conc[i] == '_') {
                    
                    prsd_sn.lig = prsd_sn.rst;  //set ligand ptr to the start of the 'rest' ptr
                    
                    if (prsd_sn.conc[3]) {      //if characters follow the conc unit substring
                        prsd_sn.rst = &prsd_sn.conc[3]; //set 'rest' ptr to start of unparsed string
                        
                    } else { //missing information in sample name, abort
                        printf("parse_sample_name: error - expected run ID to follow ligand concentration in sample name. aborting...\n");
                        abort();
                    }
                    
                    prsd_sn.conc[2] = '\0'; //terminate conc substring by replacing uscore w/ '/0'
                    prsd_sn.conc[i] = '\0'; //terminate ligand substring by replacing uscore w/ '/0'
                    tmp_ptr = prsd_sn.conc;       //set tmp ptr to conc pointer
                    prsd_sn.conc = &tmp_ptr[i+1]; //set conc pointer to start of conc substring
                    found_conc_start = 1;         //set flag that concentratio nstart was found
                }
            }
            /*
            if (prsd_sn.conc[-4] == '_'   && //check that the concentration unit
                isdigit(prsd_sn.conc[-3]) && //substring is preceded by a three-digit,
                isdigit(prsd_sn.conc[-2]) && //which is preceded by an underscore
                isdigit(prsd_sn.conc[-1]) ){
                
                prsd_sn.lig = prsd_sn.rst;   //set ligand pointer to the start of the 'rest' pointer
                prsd_sn.conc[2] = '\0';      //terminate the conc substring by replacing uscore w/ '/0'
                prsd_sn.conc[-4] = '\0';     //terminate the ligand substring by replacing uscore w/ '/0'
                tmp_ptr = prsd_sn.conc;      //set tmp ptr to conc pointer
                prsd_sn.conc = &tmp_ptr[-3]; //set conc pointer to start of conc substring
                
                if (tmp_ptr[3]) {              //if characters follow the conc unit substring
                    prsd_sn.rst = &tmp_ptr[3]; //set 'rest' pointer to start of remaining unparsed string
                    
                } else { //missing information in sample name, abort
                    printf("parse_sample_name: error - expected run ID to follow ligand concentration in sample name. aborting...\n");
                    abort();
                }
                
            } else { //unexpected formating for concentration substring
                printf("parse_sample_name: error - expected string of three digits at start of ligand concentration. aborting...\n");
                abort();
            } */
        }
        
        
        //parse run ID
        if (prsd_sn.rst[0] == 'C' && prsd_sn.rst[1] == 'A' && prsd_sn.rst[2] == 'T') { //test for CAT string
            
            if (isdigit(prsd_sn.rst[3]) && prsd_sn.rst[4] == '_') { //CAT has new formatting
                prsd_sn.run_cnt_S = &prsd_sn.rst[3];         //set run_cnt_S pointer to run count string
                prsd_sn.rst[4] = '\0';                       //terminate run count string
                prsd_sn.run_cnt_I = atoi(prsd_sn.run_cnt_S); //convert run count to int
                
                if (prsd_sn.rst[5]) {          //if there are characters downstream of the run count string
                    tmp_ptr = prsd_sn.rst;     //set tmp_ptr to the start of the 'rest' string
                    prsd_sn.rst = &tmp_ptr[5]; //set the 'rest string to the start of the unparsed string
                    
                } else { //error - sample name must contain run IDs
                    printf("parse_sample_name: error - expected additional information beyond CAT string. aborting...");
                    abort();
                }
                
            } else if (prsd_sn.rst[3] == '_') { //CAT has old formatting
                printf("parse_sample_name: error - old format CAT string detected in sample name '%s'. current format cat strings contain a digit >= 2 and <=%d following 'CAT' to indicate how many samples were concatenated. accurate TECprobe sample name parsing can only be performed on sample names of the current format. it is recommended to reprocess the input data sets so that sample names are up-to-date. aborting...\n", ipt_nm, MAX_RUNS);
                abort();
            } else { //CAT has unknown formatting
                printf("parse_sample_name: error - unknown CAT string format detected in sample name '%s'. aborting...", ipt_nm);
                abort();
            }
            
        } else { //sample name contains single run ID
            prsd_sn.run_cnt_S = &str1[0]; //set run_cnt_S pointer to "1" string
            prsd_sn.run_cnt_I = 1;        //set run_cnt_I value to 1
        }
                        
        //parse runID strings
        for (u = 0, parse_complete = 0; u < prsd_sn.run_cnt_I && !parse_complete; u++) {
            
            //find end of run ID by searching for underscore character
            for (i = 0; prsd_sn.rst[i] && prsd_sn.rst[i] != '_'; i++) { ;}
            
            prsd_sn.runID[u] = &prsd_sn.rst[0];              //set runID pointer to start of 'rest' string
            if (prsd_sn.rst[i] == '_' && prsd_sn.rst[i+1]) { //if an underscore was found w/ chars downstream
                prsd_sn.rst[i] = '\0';                       //replace uscore with '\0' to terminate runID
                tmp_ptr = &prsd_sn.rst[0];                   //set temp pointer to start of 'rest' string
                prsd_sn.rst = &tmp_ptr[i+1];                 //set 'rest' pointer to char after uscore
                
            } else {
                prsd_sn.rst = NULL;  //otherwise, set 'rest' pointer to NULL
                parse_complete = 1;  //set flag that there is nothing left to parse
            }
        }
                
        if (u != prsd_sn.run_cnt_I) { //parsed expected number of runIDs
            printf("parse_VL_sample_name: error - expected %d run IDs but detected %d. aborting...\n", prsd_sn.run_cnt_I, u);
            abort();
        }
                    
        //parse custom value fields
        if (!parse_complete) { //if parse is incomplete, the remaining string should be custom value fields
            
            for (u = 0; u < MAX_FIELDS && !parse_complete; u++) {
                
                //find end of custom field by searching for underscore character
                for (i = 0; prsd_sn.rst[i] && prsd_sn.rst[i] != '_'; i++) { ;}
                
                prsd_sn.field[u] = &prsd_sn.rst[0];              //set field pointer to start of 'rest' string
                if (prsd_sn.rst[i] == '_' && prsd_sn.rst[i+1]) { //if an underscore was found w/ chars downstream
                    prsd_sn.rst[i] = '\0';                       //replace uscore with '\0' to terminate field
                    tmp_ptr = &prsd_sn.rst[0];                   //set temp pointer to start of 'rest' string
                    prsd_sn.rst = &tmp_ptr[i+1];                 //set 'rest' pointer to char after uscore
                    
                } else {
                    prsd_sn.rst = NULL;  //otherwise, set 'rest' pointer to NULL
                    parse_complete = 1;  //set flag that there is nothing left to parse
                }
            }
            
            prsd_sn.field_cnt = u; //set custom value field count
            
            if (!parse_complete) { //entire sample name should be parsed at this point
                printf("parse_sample_name: error - unparsable string (%s) in sample name. aborting...\n", prsd_sn.rst);
                abort();
            }
        }
                
        /*** record parsed substring attributes in config structure ***/
        //store RNA name in config
        set_cfg_string(&cfg->input_name, prsd_sn.rna, 0);
                
        //store folding condition in config
        if (!strcmp(prsd_sn.fld, "CoTxn")) {
            set_TF_value("TRUE", "cotranscriptional", &cfg->cotranscriptional);
        } else if (!strcmp(prsd_sn.fld, "Equil")) {
            set_TF_value("FALSE", "cotranscriptional", &cfg->cotranscriptional);
        } else {
            printf("parse_sample_name: error - how did we get here?\n");
            abort();
        }
                
        //store chemical probe in config
        set_cfg_string(&cfg->chemical_probe, prsd_sn.prb, 0);
                
        //store run ID in config
        if (prsd_sn.run_cnt_I == 1) {       //one runID, set flag that sample is not concatenated
            set_TF_value("FALSE", "concatenated", &cfg->concatenated);
        } else if (prsd_sn.run_cnt_I > 1) { //more than one runID, set flag that sample is concatenated
            set_TF_value("TRUE", "concatenated", &cfg->concatenated);
        } else {                            //should be unreachable
            printf("parse_VL_sample_name: run ID count is less than 1. this should not be possible. aborting...\n");
        }
                
        cfg->run_count = prsd_sn.run_cnt_I; //set run count
                
        for (i = 0; i < prsd_sn.run_cnt_I; i++) { //copy runIDs to config struct
            set_cfg_string(&cfg->runID[i], prsd_sn.runID[i], 0);
        }
                
        //store smoothing flag in config
        if (prsd_sn.SM != NULL) {
            set_TF_value("TRUE", "smoothing", &cfg->smoothing);
        } else {
            set_TF_value("FALSE", "smoothing", &cfg->smoothing);
        }
                    
        //sample name contains ligand name/conc
        if (prsd_sn.lig != NULL && prsd_sn.conc != NULL) {
            
            //store ligand name in config
            set_cfg_string(&cfg->ligand_name, prsd_sn.lig, 0);
                          
            //store ligand concentration in config
            set_cfg_string(&cfg->ligand_conc, prsd_sn.conc, 0);
        }
                        
        //store custom fields in config
        cfg->field_count = prsd_sn.field_cnt;
        for (i = 0; i < prsd_sn.field_cnt; i++) {
            set_cfg_string(&cfg->field[i], prsd_sn.field[i], 0);
        }
        
        if (debug) {
            print_parsed_fields(ipt_nm, cfg); //print parsed fields from config values
        }
    }
}

/* remove_out_suffix: remove "_out" suffix from sample name */
void remove_out_suffix (char * snm)
{
    char out_sffx[5] = {"_out"}; //output directory suffix string
    char * p_out_sffx = NULL;    //pointer to the output directory suffix string
    
    //remove '_out' suffix and preceding transcript length value from sample name
    if ((p_out_sffx = strstr(snm, out_sffx)) != NULL) { //search sample name for '_out' suffix
        
        if ((((uint64_t)(p_out_sffx)) - ((uint64_t)(snm))) > 4) { //check that characters precede suffix
            if (p_out_sffx[-4] == '_'  &&  //check for three-digit transcript length
                isdigit(p_out_sffx[-3]) && //value and for an underscore that
                isdigit(p_out_sffx[-2]) && //precedes the transcript length value
                isdigit(p_out_sffx[-1])) {
                
                p_out_sffx[-4] = '\0'; //truncatetranscript length value and '_out' suffix
                
            } else { //'_out' suffix must be preceded by'_<3-digit transcript length>'
                printf("parse_sample_name: error - unexpected format for out suffix. aborting...\n");
                abort();
            }
            
        } else { //'_out' suffix must be preceded by'_<3-digit transcript length>'
            printf("parse_sample_name: error - expected >4 characters to preceded '_out' string in directory name. aborting...\n");
            abort();
        }
        
    } else { //suffix not found, throw error and abort
        printf("parse_sample_name: error - could not detect out suffix. aborting...\n");
        abort();
    }
}

/* remove_simple_suffix: remove a string of characters from the end of a TECprobe sample name */
void remove_simple_suffix(char * sn, char * str2rmv)
{
    char * sffx_ptr = NULL;
    
    if ((sffx_ptr = strstr(sn, str2rmv)) != NULL) {
        if (!sffx_ptr[strlen(str2rmv)]) {
            sn[(uint64_t)(sffx_ptr) - (uint64_t)(sn)] = '\0';
        } else {
            printf("remove_simple_suffix: error - suffix string was found but was not at the end of the string. aborting...\n");
            abort();
        }
    } else {
        printf("remove_simple_suffix: error - input string did not contain suffix. aborting...\n");
        abort();
    }
}

/* print_parsed_fields: print parsed fields from config values */
void print_parsed_fields(char * ipt_nm, tprobe_configuration * cfg)
{
    int i = 0; //general purpose index
    
    //print RNA name
    printf("\nRNA name:      %s\n", cfg->input_name);
    
    //print folding strategy
    printf("folding:       %s\n", (cfg->cotranscriptional) ? "cotranscriptional" : "equilibrium");
    
    //print probe type
    printf("probe:         %s\n", cfg->chemical_probe);
    
    //if present, print ligand info
    if (cfg->ligand_name != NULL && cfg->ligand_conc != NULL) {
        printf("ligand name:   %s\n", cfg->ligand_name);
        printf("ligand conc:   %s\n", cfg->ligand_conc);
    }
    
    //print concatenated status, run count, and run IDs
    printf("concatenated:  %s\n", (cfg->concatenated) ? "TRUE" : "FALSE");
    printf("run IDs:       %d\n", cfg->run_count);
    for (i = 0; i < cfg->run_count && i < MAX_RUNS; i++) {
        printf("ID%d:           %s\n", i+1, cfg->runID[i]);
    }
    
    //print smoothing status
    printf("smoothing:     %s\n", (cfg->smoothing) ? "TRUE" : "FALSE");
    
    //printf value fields
    printf("custom fields: %d\n", cfg->field_count);
    for (i = 0; i < cfg->field_count && i < MAX_FIELDS; i++) {
        printf("FIELD%d:        %s\n", i+1, cfg->field[i]);
    }
    
    //reconstruct sample name and compare to input sample name
    char reconstructed[MAX_LINE] = {0};
    mk_run_nm(&reconstructed[0], cfg); //reconstruct the sample name
    printf("\ninput name:    %s\n", ipt_nm);
    printf("reconstructed: %s\n\n", reconstructed);

    
}

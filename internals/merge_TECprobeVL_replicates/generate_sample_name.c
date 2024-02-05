//
//  generate_sample_name.c
//  
//
//  Created by Eric Strobel on 1/26/24.
//

#include <stdio.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../utils/io_management.h"

#include "../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"
#include "../cotrans_preprocessor/run_script_gen/MLT/mk_MLT_run_nm.h"

#include "../mkmtrx/cotrans_mtrx.h"
#include "../mkmtrx/mkmtrx_defs.h"

#include "./merge_TECprobeVL_replicates_defs.h"
#include "./merge_TECprobeVL_replicates_structs.h"

#include "generate_sample_name.h"

/* generate_sample_name: manages sample name parsing and generation*/
void generate_sample_name (sample_names * sn)
{
    int i = 0; //general purpose index
    
    //parse each sample name
    for (i = 0; i < sn->cnt; i++) {
        parse_sample_name(&sn->ipt[i][0], &sn->cfg[i]);
    }
    
    //merge the parsed sample names
    merge_sample_names(sn);
}

/* parse_sample_name: parse sample name for attributes that were specified in the TECprobe analysis config */
void parse_sample_name(char * ipt_nm, configuration_MLT * cfg)
{
    int i = 0;
    
    char snm[MAX_LINE] = {0};            //copy of the sample name to be parsed
    parsed_sample_name prsd_sn = {NULL}; //storage for parsed sample name sub-strings
    
    char out_sffx[5] = {"_out"}; //pointer to output directory suffix string
    char * p_out_sffx = NULL;    //pointer to the output directory suffix string
    
    char CoTxn[8] = {"_CoTxn_"}; //cotranscriptional folding type indicator, used as parsing anchor
    char Equil[8] = {"_Equil_"}; //equilibrium folding type indicator, used as parsing anchor
    char * scnd_srch = NULL;     //pointer used to perform a search for a second anchor match
    char * fld_typ = NULL;       //pointer to fold folding type indicator (CoTxn or Equil variables)
    char * othr_typ = NULL;      //pointer to the folding type indicator that was not identified as an anchor
    
    char * p_SM_flag = NULL;     //pointer to the smoothing flag
    
    char *tmp_ptr = NULL;        //temporary pointer
    
    //copy sample name to snm array for parsing into sub-strings
    if (snprintf(snm, MAX_LINE, "%s", ipt_nm) >= MAX_LINE) {
        printf("parse_sample_name: sample name exceeded buffer. aborting...\n");
        abort();
    }
    
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
        if (prsd_sn.conc != NULL && (((uint64_t)(prsd_sn.conc)) - ((uint64_t)(prsd_sn.rst))) > 4) {
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
            }
        } else {
            printf("parse_sample_name: error - expected string of three digits at start of ligand concentration. aborting...\n");
            abort();
        }
        
        //parse run ID
        if (prsd_sn.rst[0] == 'C' && //if input data was concatenated throw error and abort. this script
            prsd_sn.rst[1] == 'A' && //is intended to be used in place of concatenation, and should not be
            prsd_sn.rst[2] == 'T') { //with concatenated data sets
            printf("parse_sample_name: error - sample name contains the string 'CAT', indicating that reads from multiple samples were concatenated prior to analysis. merge_replicate_SM2_out should only be run on non-concatenated data sets. aborting...\n");
            abort();
            
        } else { //sample name contains single run ID
            
            //find end of run ID by searching for underscore character
            for (i = 0; prsd_sn.rst[i] && prsd_sn.rst[i] != '_'; i++) { ;}
            
            prsd_sn.runID = &prsd_sn.rst[0];                 //set runID pointer to start of 'rest' string
            if (prsd_sn.rst[i] == '_' && prsd_sn.rst[i+1]) { //if an underscore was found above
                prsd_sn.rst[i] = '\0';                       //replace uscore with '\0' to terminate runID
                tmp_ptr = &prsd_sn.rst[0];                   //set temp pointer to start of 'rest' string
                prsd_sn.rst = &tmp_ptr[i+1];                 //set 'rest' pointer to char after uscore
            } else {
                prsd_sn.rst = NULL;                          //otherwise, set 'rest' pointer to NULL
            }
        }
        
        /*** record parsed substring attributes in config structure ***/
        
        //store RNA name in config
        if (snprintf(cfg->input_name, MAX_LINE, "%s", prsd_sn.rna) >= MAX_LINE) {
            printf("parse_sample_name: error - RNA name exceeded buffer. aborting...");
            abort();
        }
        
        //store folding condition in config
        if (!strcmp(prsd_sn.fld, "CoTxn")) {
            cfg->cotranscriptional = 1;
        } else if (!strcmp(prsd_sn.fld, "Equil")) {
            cfg->cotranscriptional = 0;
        } else {
            printf("parse_sample_name: error - how did we get here?\n");
            abort();
        }
        
        //store chemical probe in config
        if (snprintf(cfg->chemical_probe, MAX_LINE, "%s", prsd_sn.prb) >= MAX_LINE) {
            printf("parse_sample_name: error - chemical probe name exceeded buffer. aborting...");
            abort();
        }
        
        //store run ID in config
        cfg->concatenated = 0;
        cfg->run_count = 1;
        if (snprintf(cfg->runID[0], MAX_LINE, "%s", prsd_sn.runID) >= MAX_LINE) {
            printf("parse_sample_name: error - run_id exceeded buffer. aborting...");
            abort();
        }
        
        //store smoothing flag in config
        if (prsd_sn.SM != NULL) {
            cfg->smoothing = 1;
        } else {
            cfg->smoothing = 0;
        }
        
        //store ligand name in config
        if (snprintf(cfg->ligand_name, MAX_LINE, "%s", prsd_sn.lig) >= MAX_LINE) {
            printf("parse_sample_name: error - ligand name exceeded buffer. aborting...");
            abort();
        }
        
        //store ligand concentration in config
        if (snprintf(cfg->ligand_conc, MAX_LINE, "%s", prsd_sn.conc) >= MAX_LINE) {
            printf("parse_sample_name: error - ligand concentration exceeded buffer. aborting...");
            abort();
        }
    }
}

/* merge_sample_names: confirm that sample name attributes match and generate merged sample name */
void merge_sample_names(sample_names *sn)
{
    int i = 0; //general purpose index
    
    //set merged sample name config attributes. in most cases,
    //attributes are set using the config from the first input
    //directory. comparisons with configs from all other input
    //directories are then performed below
    
    strcpy(sn->mrgd_cfg.input_name, sn->cfg[0].input_name);         //set RNA name
    sn->mrgd_cfg.cotranscriptional = sn->cfg[0].cotranscriptional;  //set folding type
    strcpy(sn->mrgd_cfg.chemical_probe, sn->cfg[0].chemical_probe); //set chemical probe
    sn->mrgd_cfg.concatenated = 1;                                  //set concatenated flag to TRUE
    sn->mrgd_cfg.run_count = sn->cfg[0].run_count;                  //set run count (will be incremented below)
    strcpy(sn->mrgd_cfg.runID[0], sn->cfg[0].runID[0]);             //set first run ID
    sn->mrgd_cfg.smoothing = sn->cfg[0].smoothing;                  //set smoothing flag
    strcpy(sn->mrgd_cfg.ligand_name, sn->cfg[0].ligand_name);       //set ligand name
    strcpy(sn->mrgd_cfg.ligand_conc, sn->cfg[0].ligand_conc);       //set ligand concentration
    
    //compare attributes of the first input directory (now in the merged config)
    //to attributes of all other input directories. during this process, the
    //run ID of each additional directory is added to the merged config and
    //the run_count variable with the config struct is incremented
    
    for (i = 1; i < sn->cnt; i++) { //for each input directory
        
        //check that RNA names are identical
        if (strcmp(sn->cfg[i].input_name, sn->mrgd_cfg.input_name)) {
            printf("merge_sample_names: error - samples contain discordant input names. aborting...\n");
            abort();
        }
        
        //check that folding type is identical
        if (sn->cfg[i].cotranscriptional != sn->mrgd_cfg.cotranscriptional) {
            printf("merge_sample_names: error - samples contain discordant folding type. aborting...\n");
            abort();
        }
        
        //check that chemical probe is identical
        if (strcmp(sn->cfg[i].chemical_probe, sn->mrgd_cfg.chemical_probe)) {
            printf("merge_sample_names: error - samples did not use the same chemical probe. aborting...\n");
            abort();
        }
        
        //if MAX_RUNS hasn't been reached, copy the run ID of the current input
        //directory to the merged config and increment the run_count variable
        if (sn->mrgd_cfg.run_count == MAX_RUNS) {
            printf("merge_sample_names: error - input directory count exceeds the maximum allowed number (%d). aborting...\n", MAX_RUNS);
            abort();
        } else {
            strcpy(sn->mrgd_cfg.runID[sn->mrgd_cfg.run_count], sn->cfg[i].runID[0]);
            sn->mrgd_cfg.run_count++;
            
        }
        
        //check that the smoothing flag is identical
        if (sn->cfg[i].smoothing != sn->mrgd_cfg.smoothing) {
            printf("merge_sample_names: error - application of neighboring transcript smoothing is not uniform across all samples. aborting...\n");
            abort();
        }

        //if either sample contained a ligand name, check that
        //ligand names and concentrations are identical
        if (sn->mrgd_cfg.ligand_name[0] || sn->cfg[i].ligand_name[0]) {
            if (strcmp(sn->cfg[i].ligand_name, sn->mrgd_cfg.ligand_name)) {
                printf("merge_sample_names: error - samples contain discordant ligand names. aborting...\n");
                abort();
            }
            
            if (strcmp(sn->cfg[i].ligand_conc, sn->mrgd_cfg.ligand_conc)) {
                printf("merge_sample_names: error - samples contain discordant ligand names. aborting...\n");
                abort();
            }
        }
    }
    
    mk_MLT_run_nm(&sn->mrg[0], &sn->mrgd_cfg); //generate a run name using the merged config
}

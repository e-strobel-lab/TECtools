//
//  parse_3pEnd_trgts.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/io_management.h"
#include "../../../seq_utils/revcomp.h"
#include "../../../seq_utils/seq2bin_hash.h"

#include "parse_3pEnd_trgts.h"

/* parse_3pEnd_trgts: parse targets file to obtain target ids, sequences, attributes,
 and min/max transcript lengths */
void parse_3pEnd_trgts(FILE * ifp, target * trgts, opt_3pEnd * end3p, target3p_params * trg_prms)
{
    extern int debug;				//flag to run debug mode
    
    int i = 0;
    int j = 0;
    
    char line[MAX_LINE];			//array for input line
    char *trgt_id = NULL;			//pointer to target name in target id
    char *trgt_sq = NULL;			//pointer to target sequence in target id
    char rcsq[MAX_TRGT_SQ] = {0};	//reverse complement of target sequence
    char len[LEN_CODE] = {0};		//array for storing end target length code
    char mut[MUT_CODE] = {0};		//array for storing end target mutation code
    
    char *cnstnt = NULL;			//pointer to start of constant region in target name
    int last = 0;					//index of last character in target id
    
    opt_3pEnd * crnt_3pOpt = NULL;	//pointer to opt_3pEnd struct for dereferencing void pointer in targets struct
    
    while (get_line(line, ifp) && trg_prms->cnt < BLOCK_SIZE) {
        if (debug) { printf("\nipt: %s\n", line);}	//when running debug mode, print input line
        
        //split input line into target name and sequence
        for (i = 0; line[i] != '\t' && line[i] && i < MAX_LINE; i++) { ;} //iterate to first tab
        if (!line[i] || i == MAX_LINE) { //check loop exit condition for unexpected line format
            printf("parse_3pEnd_trgts: error - unexpected end target line format (missing tab delimiter). aborting...\n");
            abort();
        } else {
            line[i] = '\0';			//replace tab with null to split line
            trgt_id = &line[0];		//set pointer to target id (starts at index 0)
            trgt_sq = &line[i+1];	//set pointer to target sequence (starts at index i+1)
            
            //check that target sequence length does not exceed seq2bin limits. in this case,
            //the target sequence input is directly used to generate the key and SEQ2BIN_MAX_KEY
            //is the tighter limit. however, code below throws an error in either case
            if (strlen(trgt_sq) > SEQ2BIN_MAX_IPT || strlen(trgt_sq) > SEQ2BIN_MAX_KEY) {
                printf("parse_3pEnd_trgts: error - input target sequence length (%lu) exceeds maximum allowed input (%d) or key (%d). aborting...\n", strlen(trgt_sq), SEQ2BIN_MAX_IPT, SEQ2BIN_MAX_KEY);
                abort();
            }
        }
        
        //check that target id and sequence don't exceed array bounds
        if (strlen(trgt_id) > (MAX_TRGT_ID-1) || strlen(trgt_sq) > (MAX_TRGT_SQ-1)) {
            printf("parse_3pEnd_trgts: error - target id (%lu chars) and sequence (%lu chars) must be <%d and <%d chars long, respectively. aborting...\n", strlen(trgt_id), strlen(trgt_sq), MAX_TRGT_ID, MAX_TRGT_SQ);
            abort();
        }
        
        //reverse complement target sequence. end target mapping is performed using
        //read 1, which will always be the reverse complement of the target sequence due
        //to sequencing library structure. it is therefore necessary to reverse complement
        //the 3' end target sequences. this step also checks that end target sequences do
        //not contain any unexpected bases.
        
        //reverse_complement returns 1 if sequence contains an N/n base.
        //the reverse_complement function itself will throw an error if
        //the sequence contains a non-ATGCN/atgcn character
        if (reverse_complement(rcsq, trgt_sq, REVCOMP)) {
            printf("parse_3pEnd_trgts: error - target sequence %s contains unexpected character. aborting...\n", trgt_sq);
        }
        
        //set trgts structure values
        strcpy(trgts[trg_prms->cnt].id, trgt_id);	//copy target name to target entry
        strcpy(trgts[trg_prms->cnt].sq, trgt_sq);	//copy target sequence to target entry
        strcpy(trgts[trg_prms->cnt].rc, rcsq);		//copy revcomp of target sequence to target entry
        trgts[trg_prms->cnt].key = &trgts[trg_prms->cnt].rc[0];
        
        //when running debug mode, print parsed target and revcomp
        if (debug) {
            printf(" id: %s\n", trgts[trg_prms->cnt].id);
            printf(" sq: %s\n", trgts[trg_prms->cnt].sq);
            printf(" rc: %s\n", trgts[trg_prms->cnt].rc);
            printf("key: %s\n", trgts[trg_prms->cnt].rc);
        }
        
        //set pointer to start of constant region in target id
        if ((last = strlen(trgt_id)) < CNSTNT_STRT_NAT) { //check that target id is an expected length
            printf("parse_3pEnd_trgts: error - unexpected short target id. aborting...\n");
            abort();
        }

        if (trgt_id[last-2] == 'A') {		 //mut code is NAT
            cnstnt = &trgt_id[last-CNSTNT_STRT_NAT]; //start of constant region is last-7
        } else if (trgt_id[last-2] == '_') { //mut code is not NAT
            cnstnt = &trgt_id[last-CNSTNT_STRT_MUT]; //start of constant region is last-11
        } else { //all target ids should have an 'A' or '_' at index last-2
            printf("parse_3pEnd_trgts: error - unexpected target id format. aborting...\n");
            abort();
        }
        
        //when running debug modee, print constant region of target id
        if (debug) {printf("con: %s\n", cnstnt);}
        
        //copy length characters from target name
        //length code should always stop at an underscore
        for (i = 0, j = 0; cnstnt[i] != '_' && j < LEN_CODE-1; i++, j++) {
            if (isdigit(cnstnt[i])) {	//length code should always be composed of digits
                len[j] = cnstnt[i];
            } else {
                printf("parse_3pEnd_trgts: error - target length code contains a non-digit character. aborting...\n");
                abort();
            }
        }
        len[j] = '\0';
        
        //throw error if loop did not end on '_' character
        if (cnstnt[i] != '_') {
            printf("parse_3pEnd_trgts: error - target id contains unexpected length code format. aborting...\n");
            abort();
        }
        
        //throw error if length code is an unexpected string length
        if (strlen(len) != LEN_CODE-1) {
            printf("parse_3pEnd_trgts: error - target id contains unexpected length code format. aborting...\n");
            abort();
        }
        
        //copy mutation type from target name
        //mut code should always stop at null or digit character for native and mutant targets, respectively
        for (i++, j = 0; cnstnt[i] && !isdigit(cnstnt[i]) && j < MUT_CODE-1; i++, j++) {
            if (isalpha(cnstnt[i])) { //mutation code should only contain alphabetic characters
                mut[j] = cnstnt[i];
            } else {
                printf("parse_3pEnd_trgts: error - target mutation code contains non-alphabetic character. aborting...\n");
                abort();
            }
        }
        mut[j] = '\0';
        
        //throw error if loop did not end on terminating null or digit character
        if (cnstnt[i] != '\0' && !isdigit(cnstnt[i])) {
            printf("parse_3pEnd_trgts: error - target id contains unexpected mutation code format. aborting...\n");
            abort();
        }
        
        //throw error if mut code is an unexpected string length
        if (strlen(mut) != MUT_CODE-1) {
            printf("parse_3pEnd_trgts: error - target id contains unexpected mutation code format. aborting...\n");
            abort();
        }
        
        //if parsing first target, set minimum transcript length
        //and length of 3' end target seed sequences
        if (trg_prms->cnt == 0) {
            trg_prms->min = atoi(len);			//set minimum transcript length
            trg_prms->sdLen = strlen(trgt_sq);	//set target sequence length
        }
        
        //set opt (void pointer) to point to opt_3pEnd struct (contains 3' end values)
        //for 3' end targets, opt will contain the transcript length value and the type
        //of end (NAT, SUB, INS, DEL). in the code below, opt (void) is dereferenced to
		//the crnt_3pOpt pointer for the purpose of readability
        trgts[trg_prms->cnt].opt = &end3p[trg_prms->cnt]; 	//point opt to corresponding 3' end values
        crnt_3pOpt = (opt_3pEnd *)trgts[trg_prms->cnt].opt; //dereference opt to crnt_3pOpt pointer
        
        //copy length to target entry
        crnt_3pOpt->len = atoi(len);

        //set mutation type value in target entry
        if (!(strcmp(mut, "NAT"))) {
            crnt_3pOpt->typ = NAT;
        } else if (!(strcmp(mut, "SUB"))) {
            crnt_3pOpt->typ = SUB;
        } else if (!(strcmp(mut, "DEL"))) {
            crnt_3pOpt->typ = DEL;
        } else if (!(strcmp(mut, "INS"))) {
            crnt_3pOpt->typ = INS;
        } else {
            printf("parse_3pEnd_trgts: error - unexpected code %s in target id. aborting...\n", mut);
            abort();
        }
        
        //when running debug mode, print length and mutation codes
        if (debug) {printf("len: %03d\nmut: %d\n", crnt_3pOpt->len, crnt_3pOpt->typ);}
        
        trg_prms->cnt++; //increment parsed target count
    }
    
    //TODO: could manage rare use cases by allowing the targets block to be expanded or
    //by allowing users to manually change allocated memory with an option
    //for now, throw error if trg_prms->cnt reaches BLOCK_SIZE
    if (trg_prms->cnt == BLOCK_SIZE) {
        printf("parse_3pEnd_trgts: error - unanticipated large number of targets (>=%d). check that end target sequences are not excessively long. ideally, end target sequences should be <=17 nts long. aborting...\n", BLOCK_SIZE);
        abort();
    }
    
    trg_prms->max = atoi(len); //set max_target_len to the length of the last target entry
    
    //check that trg_prms->max is <1000. a maximum transcript length of 999 should exceed
    //the needs of all current use cases
    if (trg_prms->max > END_MAX) {
        printf("parse_3pEnd_trgts: error - max target length (%d) is expected to be <=999. aborting...\n", trg_prms->max);
        abort();
    }
    
    //print target parameters to screen
    printf("\n\n********************  3' end target parsing  ********************\n");
    printf("the 3' end targets file contained %d sequences\n", trg_prms->cnt);
    printf("min transcript length: %3d\n", trg_prms->min);
    printf("max transcript length: %3d\n", trg_prms->max);
    printf("    end target length: %3d\n\n", trg_prms->sdLen);
}

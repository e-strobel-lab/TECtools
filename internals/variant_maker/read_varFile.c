//
//  read_varFile.c
//  
//
//  Created by Eric Strobel on 8/3/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../utils/io_management.h"
#include "../utils/gen_utils.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/basemap.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"

#include "read_varFile.h"

/* read_varFile: open file input and assign to file pointer */
int read_varFile(FILE * ifp, wt_source * wt, basemap * bmap, int mode)
{
    extern FILE * prcs_ofp;            //output file pointer for processing messages
    extern char prcs_out_nm[MAX_LINE]; //name of processing message output file
    extern char out_msg[MAX_LINE];     //output message
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    char line1[MAX_LINE] = {0}; //array to store input line, which contains sequence name
    char line2[MAX_LINE] = {0}; //array to store second input line(s), which contains pairing constraints
    char * p_name = NULL;       //pointer to start of sequence name
    char * p_wtsq = NULL;       //pointer to start of source sequence
    
    int mask_prs = 0; //variable indicating which mode set_basemap_pairs should be run in
    
    if (mode == MAKE_VARIANTS) {       //in MAKE_VARIANTS mode
        mask_prs = NO_MASK;            //do not mask pairs
    } else if (mode == MAKE_BARCODES){ //in MAKE_BARCODES mode
        mask_prs = MASK_PAIRS;         //mask pairs
    } else {                           //unrecognized mode
        printf("read_varFile: error - unrecognized mode. aborting...\n");
        abort();
    }
    
    /* get sequence name and source sequence. the first line of the variants file
     contains the name of the sequence from which the variant templates are derived.
     this line beings with the string '/name'. the second line of the variants file
     contains the sequence from which the variant templates are derived. this line
     begins with the string '/seq. */
    
    if (get_line(line1, ifp) && get_line(line2, ifp)) { //get first and second lines line, which
                                                        //contain the sequence name and wt sequence
        
        if (!memcmp(line1, "/name", 5) && line1[5] == '\t' && line1[6] && //1st line begins with '/name' and a tab
            !memcmp(line2, "/wtsq", 5) && line2[5] == '\t' && line2[6]) { //2nd line begins with '/wtsq' and a tab
            
            p_name = &line1[6]; //set pointer to 1st line value, which follows the tab
            p_wtsq = &line2[6]; //set pointer to 2nd line value, which follows the tab
            
            set_wt_seq(wt, p_name, p_wtsq); //initialize and set the wt sequence fasta structure
        
        } else {
            //1st and 2nd lines are formatted incorrectly
            printf("%s\n%s\n", line1, line2);
            printf("read_varFile: error - unexpected format for variants file. the first line of the file should contain the sequence name and the second line should contain the wt source sequence. aborting...\n");
            abort();
        }
    } else {
        //line was blank
        printf("read_varFile: error - could not get name and source sequence lines from input file. aborting...\n");
        abort();
    }
    
    
    /* get variant templates. each variant template comprises a sequence line and one or more pair
     lines. the code below checks that the sequence and pair lines have the same identifier and
     that the length of the sequence and pair lines matches that of the source sequence. if the
     variant template passes these checks, variant_template struct is initialized and the sequence
     and pair lines are stored. */
    
    char *p_id_s = NULL;           //pointer to sequence line id
    char *p_id_p = NULL;           //pointer to pair line id
    char *p_sq = NULL;             //pointer to sequence line value string
    char *p_pr = NULL;             //pointer to pair line value string
    int id_lkup[MAXREF] = {0};     //array to track used id numbers
    
    int wtsq_len = strlen(wt->sq); //length of source sequence
    int end_of_vTmp = 0;           //flag that the end of a variant template was found
    
    line1[0] = line2[0] = 0;       //zero line arrays
    
    for (i = 0; get_line(line1, ifp); i++) { //continue loop until failure to get new variant sequence line
        
        if (i == MAXREF) { //check that processing current variant template will not exceed MAXREF
            printf("read_varFile: error - number of input variant templates exceeds the maximum (%d). aborting...\n", MAXREF);
            abort();
        }
        
        if (!memcmp(line1, "/s", 2) && line1[5] == '\t' && line1[6]) {  //sequence line begins with '/s___' and a tab
                        
            line1[5] = '\0';    //split line 1 at the tab
            p_id_s = &line1[2]; //set pointer to seq line id
            p_sq = &line1[6];   //set pointer to sequence line value string
            
            validate_seq_line(p_id_s, p_sq, &id_lkup[0], wtsq_len); //validate line 1 as a sequence line
            set_basemap_seq(&bmap[i], p_id_s, p_sq, wt);            //set basemap sequence values
            
            for (j = 0, end_of_vTmp = 0; !end_of_vTmp; j++) { //run loop until end of variant template is found (# character)
                
                if (!get_line(line2, ifp)) { //check that file does not end prematurely
                    printf("read_varFile: error - expected every variant template to be terminated with a '#' character. aborting...\n");
                    abort();
                }
                
                if (line2[0] == '#') {  //found variant template termination character (#)
                    end_of_vTmp = 1;    //set end of variant template flag
                    set_basemap_pairs(&bmap[i], mask_prs); //set base pair parameters in basemap struct
                    //print_basemap(&bmap[i]);
                    
                //found sec struct line. sec struct begins with '/p___' and a tab
                } else if (!memcmp(line2, "/p", 2) && line2[5] == '\t' && line2[6]) {
                    line2[5] = '\0';    //split line 2 at the tab
                    p_id_p = &line2[2]; //set pointer to pair line id
                    p_pr = &line2[6];   //set pointer to pair line value string
                    
                    validate_pair_line(p_id_p, p_id_s, p_pr, wtsq_len); //validate line 2 as a pair line
                    add_pairs_to_basemap(&bmap[i], p_pr, p_id_p);       //set variant template pair constraint
                    
                } else {
                    //variant template file does not conform to the expected format
                    printf("%s\n", line2);
                    printf("read_varFile: error - unexpected format for variant template line\n");
                    abort();
                }
            }
            
            if (!j) { //no line after sequence line
                printf("read_varFile: error - variant templates must be associated with at least one secondary structure aborting...\n");
                abort();
            }
            
            line1[0] = line2[0] = 0; //zero line arrays
        
        } else {
            //variant template file does not conform to the expected format
            printf("%s\n", line1);
            printf("read_varFile: error - unexpected format for variant template line\n");
            abort();
        }
    }
    
    return i; //return number of variant templates
}

/* validate_seq_line: check that sequence line is valid */
int validate_seq_line(char * id, char * seq, int * id_lkup, int wtsq_len)
{
    int i = 0; //general purpose index
    
    //check whether id was used previously
    if (!id_lkup[atoi(id)]) {
        id_lkup[atoi(id)] = 1; //if not, mark that id was used
    } else {                   //if so, print error and abort
        printf("read_varFile: error - multiple variant template entries for id %s\n", id);
        abort();
    }
    
    //check that variant sequence length matches source sequence length
    if (strlen(seq) != wtsq_len) {
        printf("read_varFile: error - variant template %s length (%lu) does not match wt sequence length (%d). aborting...", id, strlen(seq), wtsq_len);
        abort();
    }
    
    //check that sequence does not exceed maximum length
    if (strlen(seq) > MAXLEN) {
        printf("read_varFILE: error - sequence for variant template %s (%lu nts long) exceeds the sequence length limit (%d nts). aborting...\n", id, strlen(seq), MAXLEN);
        abort();
    }
    
    //check that sequence string contains valid characters
    for (i = 0; seq[i] && i < (MAX_LINE-6); i++) { //MAX_LINE-6 due to lack of 5 char id + tab
        if (!isIUPACbase(seq[i]) && //char must be IUPAC base
            seq[i] != '.'        && //or a '.'
            seq[i] != '-') {        //or a '-'
            printf("read_varFile: error - variant template %s sequence contains unrecognized character %c. aborting...\n", id, seq[i]);
            abort();
        }
    }
    
    return 1;
}

/* validate_pair_line: check that pair line is valid */
int validate_pair_line(char * id_p, char * id_s, char * prs, int wtsq_len)
{
    int i = 0; //general purpose index
    
    //check that pair line id matches the associated sequence line id
    if (strcmp(id_p, id_s)) {
        printf("read_varFile: error - variant template sequence (%s) and pair (%s) line ids are associated but do not match. aborting...\n", id_s, id_p);
        abort();
    }

    //check that the secondary structure length matches source sequence length
    if (strlen(prs) != wtsq_len) {
        printf("read_varFile: error - variant template %s length (%lu) does not match wt sequence length (%d). aborting...", id_p, strlen(prs), wtsq_len);
        abort();
    }
    
    //check that sequence and pair strings contain only valid characters
    for (i = 0; prs[i] && i < (MAX_LINE-6); i++) { //MAX_LINE-6 due to lack of 5 char id + tab
        
        //check that pair constraint string char is valid
        if (prs[i] != '.' &&       //char must be a '.'
            prs[i] != '-' &&       //or a '-'
            prs[i] != '(' &&       //or a '('
            prs[i] != ')' ){       //or a ')'
            printf("read_varFile: error - variant template pair constraint containts unrecognized character %c. aborting...\n", prs[i]);
            abort();
        }
    }
 
    return 1;
}

//
//  parse_constraints.c
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"
#include "../seq_utils/test_possible_pairs.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/basemap.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"
#include "./read_vbase.h"
#include "./search_4_vbase_match.h"


#include "parse_constraints.h"

/* parse_constraints: read constraints file and store the settings for each constraint in a constraints struct */
int parse_constraints(FILE * ipt, constraints * cons, basemap * bmap, char * cnstnt_indels)
{
    extern int debug;      //flag to turn on debug mode
    
    int i = 0;    //index for constraint groups
    int bi = 0;   //index for constrained bases within a constraint group
    int pi = 0;   //index for constrained pairs within a constraint group
    
    int j = 0;    //index for input line
    int k = 0;    //index for code and tmp arrays
    
    char line[MAX_LINE+1] = {0};         //array to store input line
    char code[MAX_CODE+1] = {0};         //array to store constraint code
    char tmp_pos[MAX_POS_ARRAY+1] = {0}; //array to store position value characters
    char tmp_ins[MAX_POS_ARRAY+1] = {0}; //array to store insertion value characters
    
    char *p_base_prms = NULL; //pointer to start of base parameters in line array
        
    int new_constraint = 1;   //flag to indicate new constraint
    
    while (get_line(line, ipt)) {
        
        if (debug) { printf("\nipt = %s\n", line); } //in debug mode, print input line
        
        if (new_constraint) {
            /***** parsing the first line of a new constraint, store constraint name *****/
            
            //check that number of constraints doesn't exceed MAX_CONSTRAINTS
            if (i >= MAX_CONSTRAINTS) {
                printf("parse_constraints: error - number of supplied constraints exceeds the expected maximum (%d). aborting...\n", MAX_CONSTRAINTS);
                abort();
            }
            
            //copy name to constraints struct, replace any spaces with underscores
            //newlines are trimmed by get_line, so only internal spaces should remain
            //replacing spaces with underscores is necessary for datagraph compatibility
            for (j = 0; line[j] && j < MAX_LINE; j++) {
                cons[i].nm[j] = (isspace(line[j])) ? '_' : line[j];
            }
            cons[i].nm[j] = '\0';
            
            new_constraint = 0; //turn off new constraint flag
            
        } else {
            /********* parsing body of a constraint *********/
            
            code[0] = tmp_pos[0] = tmp_ins[0] = '\0'; //zero arrays

            /* get constraint code
             base = constrained base
             pair = constrained pair
             # = end of constraint*/
            for (j = 0, k = 0; !isspace(line[j]) && line[j] && k < MAX_CODE; j++, k++) {
                code[k] = line[j];
            }
            code[k] = '\0';
            
            if (!isspace(line[j]) && code[0] != '#') {
                printf("parse_constraints: error - unexpected format for constraint line.\n%s\naborting...\n", line);
                abort();
            }
            
            if (debug) { printf("code= %s\n", code); } //in debug mode, print constraint code
            
            /********* process line based on constraint code *********/
            if (code[0] == '#') { //reached end of constraint
                
                //add constant indels
                if (strcmp(cnstnt_indels, "none")) {                              //if constant indels string contains values
                    set_cnstnt_indel_constraints(cons, cnstnt_indels, &bi, bmap); //add constant indels to base constraints
                }
                
                //check that constraint contained the expected number of vbases
                if (bi != (bmap->vb_cnt + bmap->ci_cnt + bmap->d_cnt)) {
                    printf("parse_constraints: error - all constraints must contain the same number of variable bases. aborting...\n");
                    abort();
                }
                
                cons[i].bcnt = bi; //set base count
                cons[i].pcnt = pi; //set pair count
                
                print_base_constraints(&cons[i]);
                
                /* during constraint parsing, read_vbase checks that:
                 1. vbase matches expected format
                 2. vbase in constraint matches vbase reference seq
                 3. constant in constraint matches reference seq
                 */
                
                chk_vbase_complete(&cons[i], bmap); //check that there is a base constraint for every reference vbase
                chk_vbase_duplicates(&cons[i]);     //check that there are no constrained base duplicates
                chk_pair_vbases(&cons[i], bmap);    //check for conflicting vbase/pair constraints
                chk_pair_duplicates(&cons[i]);      //check that there are no constrained pair duplicates
                validate_pairs(&cons[i], bmap);     //check that constrained pairs are possible
                //chk_pair_overlap(&cons[i]);       //check for conflicting pair attributes in overlapping pairs
                                                    //^^TODO: add pair overlap test later

                bi = pi = 0;        //zero base and pair indices
                i++;                //increment constraint index
                new_constraint = 1; //turn on new constraint flag
                
                if (debug) {
                    printf("\n");
                }
                
            } else if (!strcmp(code, "base")) { //line specifies constrained base
                
                if (bi == MAX_VBASES) { //check that maximum number of base constraints will not be exceeded
                    printf("parse_constraints: error - number of base constraints (including constant indels) exceeds the maximum (%d). aborting...\n", MAX_VBASES);
                    abort();
                }
                
                while (isspace(line[j])) {j++;} //iterate to next non-space char
                p_base_prms = &line[j];         //set pointer to base parameters
                read_vbase(p_base_prms, '\0', &cons[i].base[bi], bmap, BASE_CONSTRAINT); //parse base parameters
                
                if (debug) {
                    printf("%d typ = %d\n", bi, cons[i].base[bi].typ);
                    printf("%d seq = %c\n", bi, cons[i].base[bi].seq);
                    printf("%d pos = %d\n", bi, cons[i].base[bi].pos);
                    printf("%d ins = %d\n", bi, cons[i].base[bi].ins);
                }
        
                bi++; //increment base index
                
            } else if (!strcmp(code, "pair")) { //line specifies constrained pair
                
                if (pi == MAX_PAIRS) { //check that maximum number of pair constraints will not be exceeded
                    printf("parse_constraints: error - number of pair constraints exceeds the maximum (%d). aborting...\n", MAX_PAIRS);
                    abort();
                }
                
                while (isspace(line[j])) {j++;} //iterate to next non-space char
                p_base_prms = &line[j];         //set pointer to first base parameters
                
                //parse base parameters of first pair member
                if (*(p_base_prms = read_vbase(p_base_prms, ',', &cons[i].pair[pi].bs[POS1], bmap, PAIR_CONSTRAINT)) == ',') {
                    p_base_prms = &p_base_prms[1]; //skip comma, set pointer to second base parameters
                } else {
                    printf("parse_constraints: error - unexpected format for constrained pair. use the format <base1>,<base2><tab><pair_attribute>. aborting...\n");
                    abort();
                }
                
                //parse base parameters of second pair member
                if (*(p_base_prms = read_vbase(p_base_prms, '\t', &cons[i].pair[pi].bs[POS2], bmap, PAIR_CONSTRAINT)) == '\t') {
                    p_base_prms = &p_base_prms[1]; //skip comma, set pointer to pair attribute string
                } else {
                    printf("parse_constraints: error - unexpected format for constrained pair. use the format <base1>,<base2><tab><pair_attribute>. aborting...\n");
                    abort();
                }
                
                //check that pair constraint does not specify a self-pair
                if (cons[i].pair[pi].bs[POS1].pos == cons[i].pair[pi].bs[POS2].pos &&
                    cons[i].pair[pi].bs[POS1].ins == cons[i].pair[pi].bs[POS2].ins) {
                    printf("parse_constraints: error - pair constraint specifies a self-pair for ");
                    print_vbase(&cons[i].pair[pi].bs[POS1]);
                    printf(". aborting...\n");
                    abort();
                }

                //set base pair constraint
                if (!strcmp(p_base_prms, "ANY_PAIR")) {
                    cons[i].pair[pi].att = ANY_PAIR;
                    
                } else if (!strcmp(p_base_prms, "WC_PAIR")) {
                    cons[i].pair[pi].att = WC_PAIR;
                    
                } else if (!strcmp(p_base_prms, "STRONG")) {
                    cons[i].pair[pi].att = STRONG;
                    
                } else if (!strcmp(p_base_prms, "WEAK")) {
                    cons[i].pair[pi].att = WEAK;
                    
                } else if (!strcmp(p_base_prms, "WEAK_AT")) {
                    cons[i].pair[pi].att = WEAK_AT;
                    
                } else if (!strcmp(p_base_prms, "WEAK_AU")) {
                    cons[i].pair[pi].att = WEAK_AU;
                    
                } else if (!strcmp(p_base_prms, "WEAK_GT")) {
                    cons[i].pair[pi].att = WEAK_GT;
                    
                } else if (!strcmp(p_base_prms, "WEAK_GU")) {
                    cons[i].pair[pi].att = WEAK_GU;
                    
                } else if (!strcmp(p_base_prms, "MISMATCH")) {
                    cons[i].pair[pi].att = MISMATCH;
                    
                } else if (!strcmp(p_base_prms, "NO_CONSTRAINT")) {
                    cons[i].pair[pi].att = NO_CONSTRAINT;
                    
                } else {
                    printf("parse_constraints: error - unexpected pair constraint (%s). use the constratint ANY_PAIR, WC_PAIR, STRONG, WEAK, WEAK_AU/WEAK_AT, WEAK_GU/WEAK_GT, MISMATCH, or NO_CONSTRAINT. aborting...\n", p_base_prms);
                    abort();
                }
                
                if (debug) {
                    printf("ps1 = %d %d %d %c\n", cons[i].pair[pi].bs[POS1].typ, cons[i].pair[pi].bs[POS1].pos, cons[i].pair[pi].bs[POS1].ins, cons[i].pair[pi].bs[POS1].seq);
                    printf("ps2 = %d %d %d %c\n", cons[i].pair[pi].bs[POS2].typ, cons[i].pair[pi].bs[POS2].pos, cons[i].pair[pi].bs[POS2].ins, cons[i].pair[pi].bs[POS2].seq);
                    printf("att = %d\n", cons[i].pair[pi].att);
                }
                
                pi++; //increment pair index
                
            } else { //unrecognized constraint code
                printf("parse_constraints: unrecognized constraint code \"%s\". aborting...", code);
                abort();
            }
        }
    }
    
    //printf("exiting parse_constraints\n");
    return i;
}

/* set_cnstnt_indel_constraints: add constant indel constraints to the constraints structure base constraints */
void set_cnstnt_indel_constraints(constraints * cons, char * cnstnt_indels, int * bi, basemap * bmap)
{
    char * p_base_prms = &cnstnt_indels[0]; //set base params pointer to start of constant indels string

    while (*p_base_prms != '\0') { //until the end of the constant indels string is reached
        
        //read the next vbase
        p_base_prms = read_vbase(p_base_prms, '_', &cons->base[(*bi)], bmap, CNST_CONSTRAINT);
    
        if (p_base_prms[0] == '_') {       //if the return value of read_vbase is a '_' delimiter
            p_base_prms = &p_base_prms[1]; //set the base params pointer to the star of the next vbase
        }
        
        (*bi)++; //increment the base index
    }
    
    return;
}

/* chk_vbase_complete: check that there is a base constraint for every reference vbase */
void chk_vbase_complete(constraints * cons, basemap * bmap)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int mtch_fnd = 0; //found match to reference vbase
    int error = 0;    //flag that match to reference vbase was not found
    
    char vb_nm[MAX_NAME+1] = {0}; //array for constructing variable base name from bmap
    
    for (i = 0; i < bmap->vb_cnt; i++) {                                       //for every reference vbase
        for (j = 0, mtch_fnd = 0; j < cons->bcnt && !mtch_fnd; j++) {          //check every base constraint for a match
            
            if (cons->base[j].typ == bmap->typ[bmap->vb_ix[i]] &&              //test type match
                cons->base[j].pos == bmap->pos[bmap->vb_ix[i]] &&              //test position match
                cons->base[j].ins == bmap->v_ins[bmap->vb_ix[i]].ins_contig && //test insertion match
                is_dgnrt_mtch(cons->base[j].seq, bmap->rS[bmap->vb_ix[i]])) {  //test degenerate sequence match
                
                mtch_fnd = 1; //set flag that match was found
            }
        }
        
        if (!mtch_fnd) { //if no match was found
            
            mk_vbase_nm(bmap, bmap->vb_ix[i], vb_nm, MAX_NAME, '\0');                          //generate variable base name
            printf("chk_vbase_complete: no match for vbase %s in base constraints.\n", vb_nm); //print error message
            error = 1;                                                                         //set error flag
        }
    }
    
    if (error) {                 //if there was an error
        printf("aborting...\n"); //print abort message
        abort();                 //and abort
    }
    
    return;
}

/* chk_vbase_duplicates: check that constraints do not contain duplicate vbase entries */
void chk_vbase_duplicates(constraints * cons)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    for (i = 0; i < cons->bcnt; i++) {                     //for every variable base
        for (j = 0; j < cons->bcnt; j++) {                 //make comparison with every other variable base
            if (cons->base[i].pos == cons->base[j].pos &&  //test for identical position number
                cons->base[i].ins == cons->base[j].ins &&  //test for identical insertion number
                i != j) {                                  //ignore self-comparison
                
                printf("chk_vbase_duplicates: error - duplicate base constraint entry. bases ");
                print_vbase(&cons->base[i]);
                printf(" and ");
                print_vbase(&cons->base[j]);
                printf(" have identical position and insertion values. aborting...\n");
                abort();
            }
        }
    }
    
    return;
}

/* chk_pair_vbases: check that bases in pair constraints match vbases
 in base constraints or constant bases in the reference sequence*/
void chk_pair_vbases(constraints * cons, basemap * bmap)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int error = 0; //flag that error was found
    
    for (i = 0; i < cons->pcnt; i++) { //for ever pair
        for (j = 0; j < 2; j++) {      //for each member of the pair
            
            if (cons->pair[i].bs[j].typ == VARIABLE) {
                
                //if the base is VARIABLE test, for an
                //exact match to a base constraint vbase
                
                if (search_4_vbase_match(&cons->base[0],  &cons->pair[i].bs[j],  cons->bcnt, EXACT_VBASE_MATCH) == NULL) {
                    
                    //no match, print error message and set error flag
                    printf("chk_pair_vbases: error - pair constraint vbase ");
                    print_vbase(&cons->pair[i].bs[j]);
                    printf(" does not match any variable base constraint.\n");
                    error = 1;
                }
                
            } else if (cons->pair[i].bs[j].typ == CONSTANT) {
                
                //if the base is CONSTANT, test for an
                //exact match to the reference sequence
                
                if (cons->pair[i].bs[j].seq != bmap->lkp0[cons->pair[i].bs[j].pos][cons->pair[i].bs[j].ins]) {
                    
                    //no match, print error message and set error flag
                    printf("chk_pair_vbases: error - pair constraint constant base ");
                    print_vbase(&cons->pair[i].bs[j]);
                    printf(" does not match the reference sequence.\n");
                    error = 1;
                }
            }
        }
    }
    
    if (error) {
        printf("aborting...\n");
        abort();
    }
    
    return;
}


/* chk_pair_duplicates: check that constraints do not contain duplicate pair entries */
void chk_pair_duplicates(constraints * cons)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    for (i = 0; i < cons->pcnt; i++) {                                       //for every variable pair
        for (j = 0; j < cons->pcnt; j++) {                                   //make comparison with every other variable pair
            if (cons->pair[i].bs[POS1].pos == cons->pair[j].bs[POS1].pos &&  //test base 1 for identical position number
                cons->pair[i].bs[POS1].ins == cons->pair[j].bs[POS1].ins &&  //test base 1 for identical insertion number
                cons->pair[i].bs[POS2].pos == cons->pair[j].bs[POS2].pos &&  //test base 2 for identical position number
                cons->pair[i].bs[POS2].ins == cons->pair[j].bs[POS2].ins &&  //test base 2 for identical insertion number
                i != j) {                                                    //ignore self-comparison
                
                printf("chk_pair_duplicates: error - duplicate pair entry. pairs ");
                print_vbase(&cons->pair[i].bs[POS1]);
                printf(",");
                print_vbase(&cons->pair[i].bs[POS2]);
                printf(" and ");
                print_vbase(&cons->pair[j].bs[POS1]);
                printf(",");
                print_vbase(&cons->pair[j].bs[POS2]);
                printf(" have identical position and insertion values. aborting...\n");
                abort();
            }
        }
    }
    
    return;
}

/* check for conflicting pair attributes in overlapping pairs */
void chk_pair_overlap(constraints * cons)
{
    //TODO: Add pair overlap test later
}



/* validate_pairs: validate base pair constraints */
void validate_pairs(constraints * cons, basemap * bmap)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int e = 0; //error index
    
    char base[2] = {0}; //array to store base character for testing
    int  pos[2]  = {0}; //array to store base position for error messages
    
    int fnd_error = 0;        //flag that error was found
    char * error[2] = {NULL}; //pointer to error message NOTE: code is written so array index 2 is never exceeded
    
    char any_nv_err[128]   = "ANY_PAIR but have been set as specific bases.";
    char any_mm_err[128]   = "ANY_PAIR but are unable to pair.";
    char wc_pair_err[128]  = "a WC_PAIR but are unable to form a Watson-Crick pair.";
    char strong_err[128]   = "a STRONG pair but are unable to form a STRONG pair.";
    char weak_err[128]     = "a WEAK pair but are unable to form a WEAK pair.";
    char weak_au_err[128]  = "a WEAK_AU pair but are unable to form a WEAK_AU pair.";
    char weak_gu_err[128]  = "a WEAK_GU pair but are unable to form a WEAK_GU pair.";
    char mismatch_err[128] = "a MISMATCH but always form a pair.";
        
    for (i = 0; i < cons->pcnt; i++) { //test every pair constraint
        
        e = 0;                                   //set error index to zero
        error[0] = error[1] = NULL;              //set error message pointers to null
        base[POS1] = cons->pair[i].bs[POS1].seq; //set base 1
        base[POS2] = cons->pair[i].bs[POS2].seq; //set base 2
        pos[POS1]  = cons->pair[i].bs[POS1].pos; //set position 1
        pos[POS2]  = cons->pair[i].bs[POS2].pos; //set position 2
        
        //NOTE: currently only ANY_PAIR requires variable bases
        
        switch (cons->pair[i].att) {
            case ANY_PAIR: //test that the bases can form a base pair and that >= 1 is a variable base
                if (isDNAbase(base[POS1]) && isDNAbase(base[POS2])) {error[e++] = &any_nv_err[0];}
                if (!test_possible_pairs(base[POS1], base[POS2], TEST_BASEPAIR, NO_TYPE_TEST)) {error[e++] = &any_mm_err[0];}
                break;
                
            case WC_PAIR:  //test that the bases can form a Watson-Crick base pair
                if (!test_possible_pairs(base[POS1], base[POS2], TEST_PAIRTYPE, TEST_AT_GC)) {error[e++] = &wc_pair_err[0];}
                break;
                
            case STRONG:   //test that the bases can form a strong base pair
                if (!test_possible_pairs(base[POS1], base[POS2], TEST_PAIRTYPE, TEST_GC)) {error[e++] = &strong_err[0];}
                break;
                
            case WEAK:     //test that the bases can form a weak pair
                if (!test_possible_pairs(base[POS1], base[POS2], TEST_PAIRTYPE, TEST_GT_AT)) {error[e++] = &weak_err[0];}
                break;
                
            case WEAK_AU:  //test that the bases can form a AU pair
                if (!test_possible_pairs(base[POS1], base[POS2], TEST_PAIRTYPE, TEST_AT)) {error[e++] = &weak_au_err[0];}
                break;
                
            case WEAK_GU:  //test that the bases can form a GU pair
                if (!test_possible_pairs(base[POS1], base[POS2], TEST_PAIRTYPE, TEST_GT)) {error[e++] = &weak_gu_err[0];}
                break;
                
            case MISMATCH: //test that the bases can be a mismatch
                if (!test_possible_pairs(base[POS1], base[POS2], TEST_MISMATCH, NO_TYPE_TEST)) {error[e++] = &mismatch_err[0];}
                break;
                
            case NO_CONSTRAINT: //no need to do anything for NO_CONSTRAINT
                break;
            
            default: //unrecognized pair attribute
                printf("validate_pairs: error - unrecognized pair attribute. aborting...\n");
                abort();
                break;
        }
        
        if (e) {           //if an error was found
            fnd_error = 1; //set found error flag
        }

        
        for (j = 0; j < e; j++) { //print error messages
            printf("validate_pairs: ERROR - nucleotides %d%c and %d%c are constrained as %s",
                   pos[POS1], base[POS1], pos[POS2], base[POS2], error[j]);
        }
    }
    
    if (fnd_error) { //if an error was found, abort
        printf("aborting...\n");
        abort();
    }
    
    return ;
}

/* print_vbase: print formatted vbase */
void print_vbase(base_params * vbase)
{
    switch (vbase->typ) {
        case CONSTANT: printf("c"); break;  //lead CONSTANT bases with 'c'
        case DELETION: printf("d"); break;  //lead DELETION bases with 'd'
        default:
            break;
    }
    printf("%d", vbase->pos);        //print the base position
    if (vbase->ins) {                //if the base is an insertion
        printf("i%d", vbase->ins);   //print 'i' and the insertion contig value
    }
    printf("%c", vbase->seq);        //print the base sequence
    
    return;
}

/* print_base_constraints: print list of all base constraints */
//NOTE: currently for debugging only, might expand in future
void print_base_constraints(constraints * cons)
{
    int i = 0; //general purpose index
    
    for (i = 0; i < cons->bcnt; i++) {          //for every base constraint
        printf("%d  ", i);                      //print the current index
        print_vbase(&cons->base[i]);            //print the current vbase
        printf("\n");                           //print newline
    }
    printf("\n");
    
    return;
}

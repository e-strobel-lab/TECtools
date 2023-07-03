//
//  read_vbase.c
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/is_dgnrt_mtch.h"
#include "../seq_utils/basemap.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"

#include "read_vbase.h"


/* read_vbase: parse variable base information from an input string and store in a base_params structure.
 read_vbase performs extensive checks to confirm the validity of the parsed variable base with respect to
 the reference sequence and variable base reference sequence.
 
 format:
   variable base, no insertion: <pos><seq>          e.g., 19N
   variable base, w/ insertion: <pos>i<ins><seq>    e.g., 19i1N
   constant base:               c<pos><seq>         e.g., c19A
   deletion site:               d<pos><seq>         e.g., d19
 */

char * read_vbase(char * vbase, char delimiter, base_params * prsd_vbase, basemap * bmap, int mode)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    char tmp_pos[MAX_POS_ARRAY+1] = {0}; //array for storing vbase position
    char tmp_ins[MAX_POS_ARRAY+1] = {0}; //array for storing insertion number
    
    uint64_t bi = 0;   //reference sequence base index of vbase that was read
    
    int typ_error = 0; //flag for vbase type error
    int pos_error = 0; //flag for vbase position error
    int ins_error = 0; //flag for vbase insertion contig error
    int seq_error = 0; //flag for vbase sequence error
    
    while (isspace(vbase[i])) {i++;} //bypass whitespace
    
    /***** set flag to indicate type of vbase (variable, deletion, constant) *****/
    //all variable bases must lead with a digit or '-'. deletion and constant
    //bases lead with 'd' and 'c', respectively, followed by a digit or '-'
    
    if (isdigit(vbase[i]) || vbase[i] == '-') { //variable bases lead with a digit or negative sign in position string
        if (mode == CNST_CONSTRAINT) {
            printf("read_vbase: error - base in constants string %s was supplied as a variable base. constant base constraints cannot be declared as variable bases. aborting...\n", vbase);
            abort();
        }
        
        prsd_vbase->typ = VARIABLE;             //set VARIABLE flag
        
    } else if (isalpha(vbase[i])) { //deletion and constant bases lead with an alphabetic char
        
        if (mode == BASE_CONSTRAINT) {  //vbases in base constraints cannot be constants or deletions
            printf("read_vbase: error - variable base in base constraint %s was supplied as a constant or deletion. variable bases used as base constraints cannot be declared as constants or deletions. aborting...\n", vbase);
            abort();
        }
        
        /* determine if base is a deletion or a constant */
        switch (vbase[i]) {
            case 'd': //deletion bases lead with 'd'
                if (mode == PAIR_CONSTRAINT) {  //vbases in pair constraints cannot be deletions
                    printf("read_vbase: error - vbases in pair constraints cannot be declared as deletions. aborting...\n");
                    abort();
                }
                prsd_vbase->typ = DELETION; //set DELETION flag
                break;
            
            case 'c': //constant bases lead with 'c'
                prsd_vbase->typ = CONSTANT; //set CONSTANT flag
                break;
                
            default: //unrecognized vbase format
                printf("read_vbase: error - unexpected vbase format in string %s. leading alphabetic character must be 'c' to indicate constant base or 'd' to indicate deletion. aborting...\n", vbase);
                abort();
                break;
        }
        
        i++; //increment i to reach the start of the position string
        if (!isdigit(vbase[i]) && vbase[i] != '-') { //check that the position string starts with a digit or '-'
            printf("read_vbase: error - unexpected format for deletion/constant base in string %s. leading 'd' (deletion) and 'c' (constant) must be followed by digit(s) (or '-' if negative base position) that specify base position. aborting...\n", vbase);
            abort();
        }
        
    } else { //vbase contains unexpected leading character
        printf("read_vbase: error - vbase in string %s contains unexpected leading character %c. aborting...\n", vbase, vbase[i]);
        abort();
    }
    
    
    /***** read and store vbase position *****/
    //vbase position is a string of digits that precedes either an IUPAC base
    //char (vbase sequence) or 'i' (which indicates the vbase is an insertion)
    for (j = 0; (isdigit(vbase[i]) || vbase[i] == '-') && vbase[i] && j < MAX_POS_ARRAY; i++, j++) {
        tmp_pos[j] = vbase[i];
    }
    tmp_pos[j] = '\0';
        
    //if MAX_POS_ARRAY was reached, check that the next character is alphabetic.
    //the next character should always be either the base sequence or an 'i' that
    //indicates that the base is an insertion
    if (j == MAX_POS_ARRAY && !isalpha(vbase[i])) { //check that complete position string was read
        printf("read_vbase: error - variable base position string (%s) exceeds maximum allowed length (%d). aborting...\n", vbase, MAX_POS_ARRAY);
        abort();
    
    //if the end of the string was reached, throw error
    } else if (!vbase[i]) {
        printf("read_vbase: error - reached the end of vbase string (%s) before processing base identity. aborting...\n", vbase);
        abort();
    }
    
    prsd_vbase->pos = atoi(tmp_pos); //set vbase position in base_params struct
        
    /***** parse insertion value *****/
    //vbases that are insertions are specified by an 'i' following the vbase position.
    //'i' precedes a string of digits that specifies the position of the insertion
    //within a series of insertions (first insertion is '1', second is '2', etc.).
    if (vbase[i] == 'i') { //vbase is an insertion
        
        if (prsd_vbase->typ == VARIABLE || prsd_vbase->typ == CONSTANT) { //only variable/constant bases can be insertions
            
            if (isdigit(vbase[++i])) { //if base that follows 'i' is a digit
                
                //read and store insertion value in base_params struct
                for (j = 0; isdigit(vbase[i]) && vbase[i] && j < MAX_POS_ARRAY; i++, j++) {
                    tmp_ins[j] = vbase[i];
                }
                tmp_pos[j] = '\0';
                
                //if MAX_POS_ARRAY was reached, check that the next character is alphabetic.
                //the next character should always be the base sequence
                if (j == MAX_POS_ARRAY && !isalpha(vbase[i])) { //check that complete insertion contig was read
                    printf("read_vbase: error - variable base insertion contig string (%s) exceeds maximum allowed length (%d). aborting...\n", vbase, MAX_POS_ARRAY);
                    abort();
                
                } else if (!vbase[i]) { //if the end of the string was reached, throw error
                    printf("read_vbase: error - reached the end of vbase string (%s) before processing base identity. aborting...\n", vbase);
                    abort();
                }
                
                prsd_vbase->ins = atoi(tmp_ins); //set insertion number in base_params struct
                
            } else {
                printf("read_vbase: error - unexpected vbase format in string %s. the insertion position that follows 'i' must be digit(s). aborting...\n", vbase);
                abort();
            }
            
        } else if (prsd_vbase->typ == DELETION) { //deletions cannot be insertions
            printf("read_vbase: error - unexpected vbase format in string %s. deletions cannot be insertions. aborting...\n", vbase);
            abort();
            
        } else { //unrecognized type, this should never happen
            printf("read_vbase: error - unexpected type for string %s. this error suggests an internal problem aborting...\n", vbase);
            abort();
        }
        
    } else {
        prsd_vbase->ins = 0; //value zero indicates not an insertion
    }

    
    /***** read and store vbase sequence in base_params struct *****/
    //check that vbase[i] is an IUPAC base and that the next character terminates the vbase string
    //permissible terminating characters are the specified delimiter, space characters, or terminating null
    if (isIUPACbase(vbase[i]) && (vbase[i+1] == delimiter || isspace(vbase[i+1]) || !vbase[i+1])) {
        prsd_vbase->seq = vbase[i];
    } else {
        printf("read_vbase: error - unexpected vbase format in string %s. expected vbase to terminate with an IUPAC base followed by a '%c', space character, or null character. terminating char: %c %d. aborting...\n", vbase, delimiter, vbase[i+1], vbase[i+1]);
        abort();
    }
    
    /***** test that deletion/constant vbases were not declared with a degenerate base sequence*/
    if ((prsd_vbase->typ == DELETION || prsd_vbase->typ == CONSTANT) && //vbase was declared as DELETION or CONSTANT
        !isDNAbase(prsd_vbase->seq) && isIUPACbase(prsd_vbase->seq)) {  //vbase is a degenerate IUPAC base
        printf("read_vbase: error - unexpected vbase format in string %s. deletion/constant vbases cannot also be degenerate bases.\n", vbase);
        abort();
    }
    
    i++; //increment to delimiter
        
    /***** validate vbase by comparison to reference sequence *****/
    //NOTE: position and insertion errors should not be possible since these
    //values are used to look up the variable base in the reference sequence
    //however, I am leaving these tests in for now.
    //NOTE: it is necessary to test the sequence of deletions against the wt
    //source sequence because the variable base reference sequence contains
    //a '-' at the deletion site.
    
    if (mode == BASE_CONSTRAINT || mode == PAIR_CONSTRAINT || mode == CNST_CONSTRAINT) {
        
        //get base index by subtracting the starting address of the
        //bmap lkp array from the address of the current variable base.
        bi = &(bmap->lkp0[prsd_vbase->pos][prsd_vbase->ins]) - *bmap->lkp_start; //use pointer subtraction to get base index
        
        typ_error = (prsd_vbase->typ != bmap->typ[bi]) ? 1 : 0;                  //test type match
        pos_error = (prsd_vbase->pos != bmap->pos[bi]) ? 1 : 0;                  //test position match
        ins_error = (prsd_vbase->ins != bmap->v_ins[bi].ins_contig) ? 1 : 0;     //test insertion match
        
        if (bmap->typ[bi] == VARIABLE || bmap->typ[bi] == CONSTANT ) {           //if variable or constant base...
            seq_error = (!is_dgnrt_mtch(prsd_vbase->seq, bmap->rS[bi])) ? 1 : 0; //test degenerate seq match using vbs ref
            
        } else if (bmap->typ[bi] == DELETION) {                                  //if deletion...
            seq_error = (prsd_vbase->seq != bmap->wt->sq[bi]) ? 1 : 0;           //test exact seq match using wt source
        }

        //if any errors were detected, print error message and abort
        if (typ_error || pos_error || ins_error || seq_error) {
            print_read_vbase_error(vbase, prsd_vbase, bmap, bi, typ_error, pos_error, ins_error, seq_error);
            abort();
        }
    }
    
    return &vbase[i]; //return pointer to delimiter/terminating character
}

/* print_read_vbase_error: print error message failure to match vbase to variable base reference sequence */
void print_read_vbase_error(char * vbase, base_params * prsd_vbase, basemap * bmap, uint64_t bi, int typ_error, int pos_error, int ins_error, int seq_error)
{
    printf("read_vbase: error - vbase in string \"%s\" does not match the reference sequence.\n", vbase);
    
    if (typ_error) {
        printf("type error:\n");
        printf("  ref=");
        switch (bmap->typ[bi]) {
            case VARIABLE: printf("VARIABLE\n"); break;
            case DELETION: printf("DELETION\n"); break;
            case CONSTANT: printf("CONSTANT\n"); break;
            case SPACER:   printf("SPACER\n"); break;
            default: printf("UNRECOGNIZED\n");break;
        }
        printf("  ipt=");
        switch (prsd_vbase->typ) {
            case VARIABLE: printf("VARIABLE\n"); break;
            case DELETION: printf("DELETION\n"); break;
            case CONSTANT: printf("CONSTANT\n"); break;
            case SPACER:   printf("SPACER\n"); break;
            default: printf("UNRECOGNIZED\n");break;
        }
    }
    
    if (pos_error) {
        printf("position error:\n");
        printf("  ref=%d\n", bmap->pos[bi]);
        printf("  ipt=%d\n", prsd_vbase->pos);
    }
    
    if (ins_error) {
        printf("insertion contig error:\n");
        printf("  ref=%d\n", bmap->v_ins[bi].ins_contig);
        printf("  ipt=%d\n", prsd_vbase->ins);
    }
    
    if (seq_error) {
        printf("sequence error:\n");
        
        if (prsd_vbase->typ == VARIABLE ||
            prsd_vbase->typ == CONSTANT) {
            printf("  ref=%c\n", bmap->rS[bi]);
            
        } else if (prsd_vbase->typ == DELETION) {
            printf("  ref=%c\n", bmap->wt->sq[bi]);
            
        } else {
            printf("it should not be possible to get here\n");
            abort();
        }
        printf("  ipt=%c\n", prsd_vbase->seq);
    }
    
    printf("aborting...\n");
    
    return;
}

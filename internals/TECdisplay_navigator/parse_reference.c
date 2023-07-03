//
//  parse_reference.c
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/basemap.h"

#include "../variant_maker/variant_maker_defs.h"
#include "../variant_maker/variant_maker_structs.h"
#include "../variant_maker/read_varFile.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"

#include "parse_reference.h"

/* parse_reference: read reference lines in constraints file and store the reference sequence, variable base
 reference sequence, and parsed vbases in wt_source and basemap structs. parsing vbases from the variable base
 reference sequence establishes the reference that vbases in each constraint and variant names will be checked
 against to establish their validity. */
int parse_reference(FILE * ipt, basemap * bmap, wt_source * wt, char ** cnstnt_indels)
{
        
    int i = 0;  //tracks absolute array position
    int j = 0;  //tracks ref seq position number, only incremented at non-insertion chars in ref seq
    
    char seq_code[5] = "/seq";                 //code indicating line is wt source sequence line
    char vbs_code[5] = "/vbs";                 //code indicating line is variable base sequence line
    char indels_code[17] = "/constant_indels"; //code indicating line is constant indels line
    
    char seq[MAX_LINE] = {0};    //array to store input reference sequence line
    char vbs[MAX_LINE] = {0};    //array to store input variable bases sequence line
    char indels[MAX_LINE] = {0}; //array to store constant indels line
    
    char *p_sq = NULL;  //pointer to the start of the reference sequence
    char *p_vb = NULL;  //pointer to the start of the variable bases sequence
    char *p_cid = NULL; //pointer to the start of the constant indels string
    
    /* copy first (ref seq), second (variable bases seq), and third
     (indels) lines from input file into seq, vbs, and indels arrays */
    if (!get_line(seq, ipt)) {
        printf("parse_reference: error - failed to get reference sequence. aborting...\n");
        abort();
    }
    
    if (!get_line(vbs, ipt)) {
        printf("parse_reference: error - failed to get reference vbases. aborting...\n");
        abort();
    }
    
    if (!get_line(indels, ipt)) {
        printf("parse_reference: error - failed to get constant indels. aborting...\n");
        abort();
    }
    
    /* test that ref seq line begins with /seq, set start of the reference sequence */
    if (!memcmp(seq, seq_code, strlen(seq_code)) && seq[strlen(seq_code)] == '\t') {
        p_sq = &seq[strlen(seq_code)+1]; //set pointer to the start of the reference sequence
    } else {
        printf("parse_reference: error - unexpected format for constraints file. reference sequence is not preceeded by code /seq<tab>. aborting...\n");
        abort();
    }
    
    /* test that variable bases seq line begins with /vbs, set start of variable bases sequence */
    if (!memcmp(vbs, vbs_code, strlen(vbs_code)) && vbs[strlen(vbs_code)] == '\t') {
        p_vb = &vbs[strlen(vbs_code)+1]; //set pointer to the start of the variable bases sequence
    } else {
        printf("parse_reference: error - unexpected format for constraints file. reference vbase sequence is not preceeded by code /vbs<tab>. aborting...\n");
        abort();
    }
    
    /* test that indels line begins with "/constant_indels", set start of constant indels sequence */
    if (!memcmp(indels, indels_code, strlen(indels_code)) && indels[strlen(indels_code)] == ':') {
        p_cid = &indels[strlen(indels_code)+1]; //set pointer to the start of the constant indels sequence
    } else {
        printf("parse_reference: error - unexpected format for constraints file. constant indels string is not preceeded by string \"/constant indels:\". aborting...\n");
        abort();
    }
    
    printf("\nseq: %s\nvbs: %s\nconstant_indels: %s\n\n", p_sq, p_vb, p_cid);
    
    //allocate memory for constant indels string
    if (((*cnstnt_indels) = malloc((strlen(p_cid)+1) * sizeof(*(*cnstnt_indels)))) == NULL) {
        printf("parse_reference: error - memory allocation for constant indels string failed. aborting...\n");
        abort();
    }
    strcpy(*cnstnt_indels, p_cid); //store constant indels string
    
    /* construct basemap */
    set_wt_seq(wt, "seq", p_sq);            //set wt source sequence values
    set_basemap_seq(bmap, "vbs", p_vb, wt); //set basemap sequence values
    print_vbases(bmap);                     //print parsed references vbases
    
    return 1;
}

//
//  id2variant.c
//  
//
//  Created by Eric Strobel on 9/2/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"
#include "../seq_utils/basemap.h"
#include "../seq_utils/isDNAbase.h"

#include "../TECdisplay_navigator/TECdisplay_navigator_defs.h"
#include "../TECdisplay_navigator/TECdisplay_navigator_structs.h"
#include "../TECdisplay_navigator/parse_reference.h"
#include "../TECdisplay_navigator/read_vbase.h"

#define MAX_VARIANT_IDS 32768   //maximum number of variant ids that can be processed

typedef struct variant_input {  //storage for input file attributes
    char nm[MAX_NAME+1];        //array to store sample name
    char fn[MAX_LINE+1];        //array to store file name
    FILE * fp;                  //pointer to values file
} variant_input;

typedef struct variant_id {     //storage for variant id and alias
    char * id;                  //pointer for allocating memory for variant id
    char * alias;               //pointer for allocating memory for variant alias
} variant_id;


/* parse_variant_ids: read and store variant ids and aliases */
int parse_variant_ids(FILE * ipt, variant_id * varID);

/* print_variants: reconstruct variant sequences from variant ids and print to file */
void print_variants(variant_id * varID, int id_cnt, basemap * bmap, char * out_dir_nm, int max_nt);

int main(int argc, char *argv[])
{
    int i = 0;                        //general purpose index
    
    variant_input ipt = {{0}};        //input file attributes
    
    int ipt_provided = 0;             //number of input files provided
    int out_nm_provided = 0;          //number of output directory names provided
    
    char out_dir_nm[MAX_NAME] = {0};  //array to store output directory name
    
    int max_nt = 0;                   //maximum nucleotide to print
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"variants",       required_argument,  0,  'v'},  //variants file input
            {"out-name",       required_argument,  0,  'o'},  //set output directory name
            {"max-nucleotide", required_argument,  0,  'x'},  //set maximum nucleotide to print
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "v:o:x:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
                
                       case 'v': //get values file
                if (!ipt_provided) {    //check that input file was not yet provided
                    strcpy(ipt.fn, argv[optind-1]);          //store file name
                    get_sample_name(argv[optind-1], ipt.nm); //get sample name from file name
                    get_file(&(ipt.fp), argv[optind-1]);     //set file pointer to variants file
                    ipt_provided++;                          //count values files provided
                } else {
                    printf("id2variant: error - more than one input file was provided. aborting...\n");
                    abort();
                }
                break;
                
            case 'o': //set output_directory name
                if (!out_nm_provided) { //check that output directory was not previously provided
                    if (strlen(argv[optind-1]) < MAX_NAME) { //check output directory name length
                        strcpy(out_dir_nm, argv[optind-1]);  //store output directory name
                    } else {
                        printf("id2variant: error - output directory name is longer than the maximum length (%d). setting output directory to default name 'out'.\n", MAX_NAME);
                    }
                    out_nm_provided++; //increment output name counter
                    
                } else {
                    printf("id2variant: error - more than one output directory name was provided. aborting\n");
                    abort();
                }
                break;
                
            case 'x': //set maximum nucleotide limit
                for (i = 0; argv[optind-1][i]; i++) {  //check that max nucleotide limit string is digits
                    if (!isdigit(argv[optind-1][i])) {
                        printf("id2variant: error - maximum nucleotide limit must be composed of digits. aborting...\n");
                        abort();
                    }
                }
                max_nt = atoi(argv[optind-1]); //set max nucleotide limit
                break;

            default:
                printf("error: unrecognized option. Aborting program...\n");
                abort();
        }
    }
    
    if (optind < argc) {
        printf("\nnon-option ARGV-elements:\n");
        while (optind < argc)
            printf("\n%s \n", argv[optind++]);
        putchar('\n');
        printf("Aborting program.\n\n");
        abort();
    }
    /*********** end of option parsing ***********/

    /* make output directory */
    if (!out_dir_nm[0]) {          //if no output directory name was provided
        strcpy(out_dir_nm, "out"); //if not, set output directory name to "out"
    }
    mk_out_dir(out_dir_nm);        //make output directory
    
    
    wt_source wt = {0};           //wt source sequence
    basemap bmap = {0};           //basemap for storing variable base reference sequence
    char * cnstnt_indels = NULL;  //storage for constant indels string
    parse_reference(ipt.fp, &bmap, &wt, &cnstnt_indels);  //construct basemap from reference sequence
    

    int id_cnt = 0;               //number of variant ids in input file
    variant_id * varID = NULL;    //storage for variant ids and aliases
    
    //allocate memory for variant ids
    if ((varID = calloc(MAX_VARIANT_IDS, sizeof(*(varID)))) == NULL) {
        printf("id2variant: error - memory allocation for variant ids failed. aborting...\n");
        abort();
    }
    id_cnt = parse_variant_ids(ipt.fp, varID); //parse variant id input file
    
    print_variants(varID, id_cnt, &bmap, out_dir_nm, max_nt); //reconstruct variant sequences and print to file
}

/* parse_variant_ids: read and store variant ids and aliases */
int parse_variant_ids(FILE * ipt, variant_id * varID)
{
    int i = 0;                    //general purpose index
    
    int id_cnt = 0;               //number of variant ids in input file
    char line[MAX_LINE] = {0};    //array to store lines of input file
    
    char * p_varID = NULL;        //pointer to variant id string
    char * p_alias = NULL;        //pointer to alias string
    
    while (get_line(line, ipt)) { //while there are still lines left to parse
        p_varID = &line[0];       //set variant id pointer to point to start of line
        p_alias = NULL;           //set alias pointer to NULL
        
        //iterate to the end of the line string or until a tab is reached
        //while checking that the variant id contains only permissible
        //characters. if the line contains a tab, it should also contain
        //an alias. if there is no tab, the line only contains a variant id
        
        for (i = 0; line[i] && line[i] != '\t' && line[i] != '\r'; i++) { //until a null, tab, or CR is reached
            if (!isdigit(line[i])   &&  //char must be a digit
                !isDNAbase(line[i]) &&  //or a DNA base
                line[i] != '-'      &&  //or a dash
                line[i] != 'i'      &&  //or an 'i'
                line[i] != '_') {       //or an underscore
                printf("parse_variant_ids: variant id contains unexpected character %c (ascii: %d). aborting...\n", line[i], line[i]);
                abort();
            }
        }
        
        //remove carriage return
        if (line[i] == '\r') {
            line[i] = '\0';
        }
        
        //if a tab was found, process the alias
        if (line[i] == '\t') {  //if a tab was found
            line[i] = '\0';     //set the tab to a term null
            
            if (line[i+1] && line[i+1] != '\r') { //if there is a char after the tab
                p_alias = &line[i+1]; //set the alias pointer to 1 char after the tab
                
                for (i = 0; p_alias[i] && p_alias[i] != '\r'; i++) { //check that chars are permissible
                    if (!isalnum(p_alias[i]) &&  //char must be alphanumeric
                        p_alias[i] != '_'    &&  //or an underscore
                        p_alias[i] != '-'    &&  //or a dash
                        p_alias[i] != '.') {     //or a period
                        printf("parse_variant_ids: variant alias contains unexpected character %c (ascii: %d). aliases should only contain alphanumeric, underscore, dash, and period characters. aborting...\n", p_alias[i], p_alias[i]);
                        abort();
                    }
                }
                
                //remove carriage return
                if (p_alias[i] == '\r') {
                    p_alias[i] = '\0';
                }
                
                if (i == 0) {       //if i == 0, there was no alias string
                    p_alias = NULL; //set p_alias to NULL;
                }
            }
        }
        
        //store variant id
        if (strlen(p_varID) > 0) { //check that variant id contains characters
    
            //allocate memory for variant id
            if ((varID[id_cnt].id = malloc((strlen(p_varID)+1) * sizeof(*(varID[id_cnt].id)))) == NULL) {
                printf("parse_variant_ids: error - memory allocation for variant id failed. aborting...\n");
                abort();
            }
            strcpy(varID[id_cnt].id, p_varID); //store variant id
            
        } else {
            printf("parse_variant_ids: error - variant id string contains no characters. aborting...\n");
            abort();
        }
        
        //store alias
        if (p_alias != NULL) { //if the current variant id is associated with an alias
            
            //allocate memory for alias
            if ((varID[id_cnt].alias = malloc((strlen(p_alias)+1) * sizeof(*(varID[id_cnt].alias)))) == NULL) {
                printf("parse_variant_ids: error - memory allocation for variant alias failed. aborting...\n");
                abort();
            }
            strcpy(varID[id_cnt].alias, p_alias); //store alias
        }
        
        //printf("%s\n%s\n\n", varID[id_cnt].id, varID[id_cnt].alias);
        
        id_cnt++;  //increment the variant id count
    }
    
    return id_cnt; //return the variant id count
}

/* print_variants: reconstruct variant sequences from variant ids and print to file */
void print_variants(variant_id * varID, int id_cnt, basemap * bmap, char * out_dir_nm, int max_nt)
{
    FILE * lfp = NULL;            //list file pointer
    FILE * ofp = NULL;            //individual output file pointer
    char out_nm[MAX_LINE] = {0};
    
    int i = 0;  //general purpose index
    int j = 0;  //general purpose index
    
    //during variant reconstruction variable bases within the reference sequence
    //are changed to the bases specified by the variant id. it is therefore
    //necessary to keep a copy of source reference sequence and restore the
    //reference sequence in the basemap structure before reconstructing each
    //variant sequence
    
    char ref_src[MAX_LINE] = {0};   //storage for source reference sequence
    strcpy(ref_src, bmap->rS);      //store reference sequence
    
    char * p_vbase = NULL;          //pointer to variable base in variant id string
    base_params vbPrms = {0};       //storage for parsed variable base parameters
    
    sprintf(out_nm, "./%s/%s_list.txt", out_dir_nm, out_dir_nm);
    
    //open list output file
    if ((lfp = fopen(out_nm, "w")) == NULL) {
        printf("print_variants: error - could not open output file. Aborting program...\n");
        abort();
    }
    
    for (i = 0; i < id_cnt; i++) {  //for every variant id
        
        strcpy(bmap->rS, ref_src);  //restore the base map reference sequence
        p_vbase = &varID[i].id[0];  //set pointer to start of the variant id
        
        sprintf(out_nm, "./%s/%s.txt", out_dir_nm, (varID[i].alias != NULL) ? varID[i].alias : varID[i].id);
        
        //open output file for current variant
        if ((ofp = fopen(out_nm, "w")) == NULL) {
            printf("print_variants: error - could not open output file. Aborting program...\n");
            abort();
        }
        
        while (*p_vbase) { //while there is another variable base to parse
            
            p_vbase = read_vbase(p_vbase, '_', &vbPrms, bmap, NAME); //parse the variable base
            bmap->lkp0[vbPrms.pos][vbPrms.ins] = vbPrms.seq; //edit the ref seq to match the variant id
                        
            if (*p_vbase == '_') {     //if an underscore delimiter was found, there is another vbase
                p_vbase = &p_vbase[1]; //set the vbase pointer to the start of the next vbase
            }
        }
        
        //check that ref seq contains only non-degenerate DNA bases
        for (j = 0; bmap->rS[j]; j++) {
            if (!isDNAbase(bmap->rS[j])) {
                printf("print_variants: error - incomplete variant reconstruction for reference sequence %s.\nreconstructed sequence contains a %c base:\n%s\naborting...\n", varID[i].id, bmap->rS[j], bmap->rS);
                abort();
            }
        }
        
        //print contents of output file
        if (varID[i].alias != NULL) {
            fprintf(lfp, ">%s_%s\n", varID[i].alias, varID[i].id);
            fprintf(ofp, ">%s_%s\n", varID[i].alias, varID[i].id);
        } else {
            fprintf(lfp, ">%s\n", varID[i].id);
            fprintf(ofp, ">%s\n", varID[i].id);
        }
        
        if (max_nt) {
            for (j = 0; j < max_nt && bmap->rS[j]; j++) {
                fputc(bmap->rS[j], lfp);
                fputc(bmap->rS[j], ofp);
            }
            fputc('\n', lfp);
            fputc('\n', ofp);
            
        } else {
            fprintf(lfp, "%s\n", bmap->rS);
            fprintf(ofp, "%s\n", bmap->rS);
        }
        
        /* close current variant output file */
        if (fclose(ofp) == EOF) {
            printf("print_variants: error - error occurred when closing output file. Aborting program...\n");
            abort();
        }
    }
    
    /* close list output file */
    if (fclose(lfp) == EOF) {
        printf("print_variants: error - error occurred when closing output file. Aborting program...\n");
        abort();
    }
}

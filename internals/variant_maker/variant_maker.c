//
//  make_variants.c
//  
//
//  Created by Eric Strobel on 4/15/20.
//

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../utils/io_management.h"
#include "../utils/gen_utils.h"
#include "../seq_utils/ispair.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/t2u.h" //may not need this
#include "../seq_utils/basemap.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"
#include "./read_varFile.h"
#include "./count_variants.h"
#include "./expand_variant_template.h"
#include "./filter_variants.h"

#include "make_barcodes.h"

void check_input(names * nm, int varFile_supplied, int brcdFile_supplied, int append_barcode, char * out_dir);
void print_output(names * nm, basemap * bmap, int vTmpCnt, int varCnt, char * out_dir);

//variant storage
fasta *vrnts = NULL; //pointer to fasta structures used when generating variants
int v_indx = 0;      //index of current variant

FILE * prcs_ofp = NULL;           //output file pointer for processing messages
char prcs_out_nm[MAX_LINE] = {0}; //name of processing message output file
char out_msg[MAX_LINE] = {0};     //output message

/* end of global variables */

int main(int argc, char *argv[])
{
    extern fasta *vrnts;       //pointer to fasta structures used when generating variants
    extern int v_indx;         //index of current variant
    
    extern FILE * prcs_ofp;            //output file pointer for processing messages
    extern char prcs_out_nm[MAX_LINE]; //name of processing message output file
    extern char out_msg[MAX_LINE];     //output message
    
    FILE * fp_var = NULL;      //pointer to variant parameters file
    FILE * fp_brcd = NULL;     //pointer to barcode sequences file
    
    struct names nm = {{0}};   //input names storage
    
    int mode = 0;              //run mode
    int varFile_supplied = 0;  //flag that variant template file was provided
    int brcdFile_supplied = 0; //flag that barcode file was provided
    int brcds2mk = 0;          //number of barcodes to make
    int append_barcode = 0;    //flag to append barcodes to variants
        
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    /* ******* parse options using getopt_long ******* */
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"mk-variants",     required_argument,  0,  'v'},  //variant template input file, sets MAKE_VARIANTS mode
            {"mk-barcodes",     required_argument,  0,  'b'},  //flag to make barcode file
            {"append_barcode",  required_argument,  0,  'a'},  //flag to append barcodes to variants
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "v:b:a:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            
            //variant template input file
            case 'v':
                mode = MAKE_VARIANTS;                       //set mode to MAKE_VARIANTS
                get_file(&fp_var, argv[optind-1]);          //open variant template file
                strcpy(nm.iptV, argv[optind-1]);            //copy variant template file name to names struct
                get_sample_name(argv[optind-1], nm.vTmp);   //get variant template sample name
                varFile_supplied++;                         //increment number of variant template files supplied
                break;
            
            //run barcode generation
            case 'b':
                mode = MAKE_BARCODES;                       //set mode to MAKE_BARCODES
                brcds2mk = atoi(argv[optind-1]);            //get number of barcodes to make
                
                //check that brcds2mk does not exceed barcode generation cap
                if (brcds2mk > MAX_BRCDS_2_MK) {
                    printf("main: error - the maximum number of barcodes to make is currently %d. aborting...\n", MAX_BRCDS_2_MK);
                    abort();
                }
                
                mk_brcds(brcds2mk);                         //generate barcodes 
                return 1;                                   //end program after barcodes are made
                break;                                      //leaving this in case I decide to remove the return later
            
            //append pre-generated barcodes to variants
            case 'a':
                append_barcode = 1;                         //set flag to append barcode,
                get_file(&fp_brcd, argv[optind-1]);         //open barcode file
                strcpy(nm.iptB, argv[optind-1]);            //copy barcode file name to names struct
                get_sample_name(argv[optind-1], nm.brcd);   //get barcode file name
                brcdFile_supplied++;                        //increment number of barcode files supplied
                break;
                
            default: printf("error: unrecognized option. Aborting program...\n"); abort();
        }
    }
    /* ******* end of options parsing ******* */
    
    /* ******* make output directory ******* */
    char out_dir[512] = {0};             //array to store output directory name
    sprintf(out_dir, "%s_out", nm.vTmp); //construct output directory name
    mk_out_dir(out_dir);                 //make output directory
    
    //generate input details file
    sprintf(prcs_out_nm, "./%s/%s_processing.txt", out_dir, nm.vTmp);
    if ((prcs_ofp = fopen(prcs_out_nm, "w")) == NULL) {
        printf("main: ERROR - could not open processing messages file. Aborting program...\n");
        abort();
    }
    
    check_input(&nm, varFile_supplied, brcdFile_supplied, append_barcode, out_dir); //check and print input files
    
    /* ******* make variants ******* */
    struct wt_source wt  =  {0};         //storage for wt name and sequence
    struct basemap bmap[MAXREF] = {{0}}; //basemap structures for storing variant template information
    int vTmpCnt = 0;                     //number of variant templates
    
    vTmpCnt = read_varFile(fp_var, &wt, &bmap[0], MAKE_VARIANTS); //parse input file
    
    //allocate memory for the variants encoded by the supplied variant templates
    //memory is allocated for the calculated number of variants encoded by all
    //variant templates plus extra space for each variant template reference
    int varCnt = count_variants(&bmap[0], vTmpCnt, MAKE_VARIANTS);
    
    if ((vrnts = calloc(varCnt + vTmpCnt + 1, sizeof(*vrnts))) == NULL) {
        printf("main: error - variant memory allocation failed\n");
        return 0;
    }
    
    char var_outpt[MAXLEN+1] = {0};     //array to store output sequence during variant generation
    struct trgt_vb var_lcl_bases = {0}; //variable bases used to construct local variant output sequence
    char ref_name[MAX_LINE+1] = {0};
    
    sprintf(out_msg, "\nmaking variants...\n\n");
    printf2_scrn_n_fl(prcs_ofp, out_msg);
    
    for (i = 0; i < vTmpCnt; i++) {     //for every variant template
        sprintf(out_msg, "expanding ref\t%s\n", bmap[i].rS);
        printf2_scrn_n_fl(prcs_ofp, out_msg);
        
        sprintf(ref_name, "REF:%s_%s", bmap[i].wt->nm, bmap[i].nm); //construct reference name
        
        vrnts[v_indx].nm = malloc(((strlen(ref_name)) * sizeof(vrnts[v_indx].nm))+1);   //allocate mem for reference name
        vrnts[v_indx].sq = malloc(((strlen(bmap[i].rS)) * sizeof(vrnts[v_indx].sq))+1); //allocate mem for reference seq
        
        if (vrnts[v_indx].nm == NULL || vrnts[v_indx].sq == NULL) {                     //check memory allocation success
            printf("main: reference sequence memory allocation failed. aborting...\n");
            abort();
        }
        
        strcpy(vrnts[v_indx].nm, ref_name);   //set reference name as "REF"
        strcpy(vrnts[v_indx].sq, bmap[i].rS); //store reference sequence
        v_indx++;                             //increment variant index
        
        //expand variant template
        expand_variant_template(&bmap[i], 0, var_outpt, var_lcl_bases, MAKE_VARIANTS);
        
        if (bmap[i].cnt[EXPANDED] == bmap[i].cnt[CALCULATED]) { //check success of variant template expansion
            
            //print the number of variants after expansion and filtering
            sprintf(out_msg, "variants made:\t%d\nvariants kept:\t%d\n\n", bmap[i].cnt[EXPANDED], bmap[i].cnt[FILTERED]);
            printf2_scrn_n_fl(prcs_ofp, out_msg);

        } else {
            printf("variant_maker: error - expected variant template %d to expand into %d variants before filtering but expansion generated %d variants. aborting...\n", i+1, bmap[i].cnt[CALCULATED], bmap[i].cnt[EXPANDED]);
            abort();
        }
    }
    
    print_output(&nm, &bmap[0], vTmpCnt, varCnt, out_dir);
    sprintf(out_msg, "\n%d variants were generated from %d variant template(s)\n", v_indx-vTmpCnt, vTmpCnt);
    printf2_scrn_n_fl(prcs_ofp, out_msg);
    /* ******* end of variant generation ******* */
    
    //close output file
    if (fclose(prcs_ofp) == EOF) {
        printf("main: error - error occurred when closing processing messages output file. Aborting program...\n");
        abort();
    }
}

/* check_input: check that correct input files were supplied for
 variant generation and print input file information to file*/
void check_input(names * nm, int varFile_supplied, int brcdFile_supplied, int append_barcode, char * out_dir)
{
    FILE * out_fp = NULL;                                       //output file pointer
    char out_nm[512] = {0};                                     //output file name
    sprintf(out_nm, "./%s/%s_input.txt", out_dir, nm->vTmp);    //construct input details file name
    
    //generate input details file
    if ((out_fp = fopen(out_nm, "w")) == NULL) {
        printf("check_input: ERROR - could not generate input details file. Aborting program...\n");
        abort();
    }
    
    //check that one variant templates file was supplied
    if (varFile_supplied == 1) {
        fprintf(out_fp, "variants: %s\n", nm->iptV);
    } else {
        printf("check_input: error - incorrect number of variant template files (%d). expected 1. aborting...\n", varFile_supplied);
        abort();
    }
    
    //check that one or fewer barcodes file were supplied
    if (brcdFile_supplied == 1) {
        fprintf(out_fp, "barcodes: %s\n", nm->iptB);
    } else if (brcdFile_supplied > 1) {
        printf("check_input: error - incorrect number of barcode files (%d). expected 1. aborting...\n", brcdFile_supplied);
        abort();
    }
    
    //close input details file
    if (fclose(out_fp) == EOF) {
        printf("print_input: error - error occurred when closing variant output file. Aborting program...\n");
        abort();
    }
    
    return;
}

/* print_output: print variants (without barcode) to output file */
void print_output(names * nm, basemap * bmap, int vTmpCnt, int varCnt, char * out_dir)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
    int c = 0; //constant insertion index
    int d = 0; //constant deletion index
    
    int printed_vb = 0; //flag that first vb of a variant template was printed
    
    int filtered_tot = 0; //total number of variants after filtering
    
    for (i = 0; i < MAXREF; i++) {
        filtered_tot += bmap[i].cnt[FILTERED];
    }
    
    FILE * out_fp = NULL;                                       //output file pointer
    char out_nm[512] = {0};                                     //output file name
    sprintf(out_nm, "./%s/%s_variants.txt", out_dir, nm->vTmp); //construct variant output file name
    
    char checkname[MAX_LINE] = {0};
    
    char vb_nm[MAX_VBASE_NAME_FIELD+1]={0}; //array for assembling variant attributes during name construction
    
    //generate output file
    if ((out_fp = fopen(out_nm, "w")) == NULL) {
        printf("print_output: ERROR - could not generate variants output file. Aborting program...\n");
        abort();
    }
        
    fprintf(out_fp, "variants:%d\n", filtered_tot); //print the number of variants in the output file (excluding reference variants)
    
    fprintf(out_fp, "WT:%s\t%s\n", bmap[0].wt->nm, bmap[0].wt->sq); //print wt name and sequence
    
    //print variants to the output file
    for (i = 0, j = 0; i < v_indx; i++) {
        if (strstr(vrnts[i].nm, "REF") != NULL) {
            
            //printing a reference sequence, append the number of targets
            //that were generated by that reference to "REF"
            if (j < vTmpCnt) { //check that basemap array bounds are not exceeded
                
                sprintf(checkname, "REF:%s_%s", bmap[j].wt->nm, bmap[j].nm);
                if (strcmp(vrnts[i].nm, checkname)) {
                    printf("print_output: error - unexpected reference target name. aborting...\n");
                    abort();
                }
                
                fprintf(out_fp, "%s|TPR:%d", vrnts[i].nm, bmap[j].cnt[FILTERED]); //print transcripts per reference
                
                //print list of variable bases in current reference
                fprintf(out_fp, "|VBS:");
                for (k = 0, printed_vb = 0; bmap[j].nts[k] && k < MAXLEN; k++) {
                    if (bmap[j].nts[k] == '*') {  //if nucletide is variable
                        if (printed_vb) {         //if a variable base has already been printed
                            fprintf(out_fp, "_"); //print an underscore delimiter
                        } else {                  //otherwise
                            printed_vb = 1;       //set flag that the first variable base was printed
                        }
                        
                        mk_vbase_nm(&bmap[j], k, vb_nm, MAX_VBASE_NAME_FIELD+1, '\0'); //assemble vbase name
                        fprintf(out_fp, "%s", vb_nm); //print vbase name to file
                    }
                }
                
                if (bmap[j].ci_cnt || bmap[j].d_cnt) {
                    fprintf(out_fp, "|const:");
                    
                    //print constant insertions
                    for (c = 0; c < bmap[j].ci_cnt; c++) {
                        mk_vbase_nm(&bmap[j], bmap[j].ci_ix[c], vb_nm, MAX_VBASE_NAME_FIELD+1, '\0'); //assemble vbase name
                        fprintf(out_fp,"%s", vb_nm);                      //print vbase_name to file
                        if (c+1 < bmap[j].ci_cnt || bmap[j].d_cnt) { //check if there are more constants to print
                            fprintf(out_fp, "_");                    //if so, print '_' delimiter
                        }
                        
                    }
                    
                    //print constant deletions
                    for (d = 0; d < bmap[j].d_cnt; d++) {
                        mk_vbase_nm(&bmap[j], bmap[j].d_ix[d], vb_nm, MAX_VBASE_NAME_FIELD+1, '\0'); //assemble vbase name
                        fprintf(out_fp, "%s", vb_nm);    //print vbase_name to file
                        if (d+1 < bmap[j].d_cnt) { //check if there are more deletions to print
                            fprintf(out_fp, "_");  //if so, print '_' delimiter
                        }
                    }
                }
                
                fprintf(out_fp, "\t%s\n", vrnts[i].sq); //print reference sequence
                j++;                                    //increment basemap index
                
            } else {
                printf("print_output: error - >%d reference targets were generated for %d variant templates. aborting...\n", j, vTmpCnt);
                abort();
            }
        } else {
            //printing a variant sequence
            fprintf(out_fp, "%s\t%s\n", vrnts[i].nm, vrnts[i].sq);
        }
    }
    
    //close output file
    if (fclose(out_fp) == EOF) {
        printf("print_output: error - error occurred when closing variant output file. Aborting program...\n");
        abort();
    }
    
    return;
}



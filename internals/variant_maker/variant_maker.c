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
#include "read_bcFile.h"
#include "print_output.h"

/* check_input: check that correct input files were supplied for
 variant generation and print input file information to file*/
void check_input(names * nm, int varFile_supplied, int brcdFile_supplied, int append_barcode, char * out_dir);

//variant storage
fasta *vrnts = NULL;  //pointer to fasta structures used when generating variants
uint64_t v_indx = 0;  //index of current variant

FILE * prcs_ofp = NULL;           //output file pointer for processing messages
char prcs_out_nm[MAX_LINE] = {0}; //name of processing message output file
char out_msg[MAX_LINE] = {0};     //output message

//priming sites
//const char c3sc1[34] = "ATTCGGTGCTCTTCTCTTCGGCCTTCGGGCCAA";
//const char vra3[27]  = "GATCGTCGGACTGTAGAACTCTGAAC";
const char c3sc1[34] = "attcggtgctcttctcttcggccttcgggccaa";
const char vra3[27]  = "gatcgtcggactgtagaactctgaac";

/* end of global variables */

int main(int argc, char *argv[])
{
    extern fasta *vrnts;       //pointer to fasta structures used when generating variants
    extern uint64_t v_indx;    //index of current variant
    
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
    int append_priming = 0;    //flag to append C3-SC1 and VRA3 priming sites
    int make_fasta = 0;        //flag to make fasta file
    
    char cstm_lnkr[MAX_LINKER+1] = {0}; //custom linker sequence
    
    char usr_resp[4] = {0};         //storage for user response
    char discard[MAX_LINE+1] = {0}; //array for flushing stdin
    int resp_provided = 0;          //flag that valid reponse was provided
        
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int ret = 0; //variable for storing snprintf return value
    
    /* ******* parse options using getopt_long ******* */
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"mk-variants",     required_argument,  0,  'v'},  //variant template input file, sets MAKE_VARIANTS mode
            /* WARNING: mk-barcodes mode should only be run on systems with >= 64 GB of RAM */
            {"mk-barcodes",     required_argument,  0,  'b'},  //flag to make barcode file
            {"append_priming",  no_argument,        0,  'p'},  //append C3-SC1 and VRA3 priming sites
            {"append_barcode",  required_argument,  0,  'a'},  //flag to append barcodes to variants
            {"custom_linker",   required_argument,  0,  'l'},  //use custom linker/exclude linker
            {"make_fasta",      no_argument,        0,  'f'},  //make fasta file
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "v:b:pa:l:f", long_options, &option_index);
        
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
                if (mode) {
                    printf("variant_maker: ERROR - cannot run variant maker in more than one mode. aborting...\n");
                    abort();
                }
                
                printf("\nWARNING: barcode generation should only be performed on machines with >=64GB of RAM.\n\nproceed with barcode generation (yes/no)? \n");
                
                while (!resp_provided) {
                    scanf("%3s", usr_resp);
                    
                    if (!strcmp(usr_resp, "yes")) {
                        printf("proceeding with barcode generation\n\n");
                        resp_provided = 1;
                        
                    } else if (!strcmp(usr_resp, "no")) {
                        printf("aborting barcode generation\n\n");
                        abort();
                        
                    } else {
                        printf("invalid response. please enter \"yes\" or \"no\".\n");
                        get_line(discard, stdin);
                        
                    }
                }
                
                mode = MAKE_BARCODES;                       //set mode to MAKE_BARCODES
                brcds2mk = atoi(argv[optind-1]);            //get number of barcodes to make
                
                //check that brcds2mk does not exceed barcode generation cap
                if (brcds2mk > MAX_BRCDS_2_MK) {
                    printf("main: error - the maximum number of barcodes to make is currently %d. aborting...\n", MAX_BRCDS_2_MK);
                    abort();
                }
                
                mk_brcds(brcds2mk); //generate barcodes
                return 1;           //end program after barcodes are made
                break;              //leaving this in case I decide to remove the return later
                
            //append C3-SC1 and VRA3 priming sites to variants
            case 'p':
                append_priming = 1;
                break;
            
            //append pre-generated barcodes to variants
            case 'a':
                append_barcode = 1;                         //set flag to append barcode,
                get_file(&fp_brcd, argv[optind-1]);         //open barcode file
                strcpy(nm.iptB, argv[optind-1]);            //copy barcode file name to names struct
                get_sample_name(argv[optind-1], nm.brcd);   //get barcode file name
                brcdFile_supplied++;                        //increment number of barcode files supplied
                break;
            
            //set custom linker
            case 'l':
                if (strlen(argv[optind-1]) <= MAX_LINKER) { //check that custom linker is not too long
                    strcpy(cstm_lnkr, argv[optind-1]);      //store custom linker string
                } else {
                    printf("variant_maker: custom linker string exceeds maximum length (%d characters), aborting...\n", MAX_LINKER);
                    abort();
                }
                
                if (strcmp(cstm_lnkr, "exclude")) {  //if the custom linker argument is not 'exclude'
                    for (i = 0; cstm_lnkr[i]; i++) {    //check that the custom linker sequence only
                        if (!isDNAbase(cstm_lnkr[i])) { //contains DNA bases
                            printf("variant_maker: ERROR - custom linkers should only contain DNA bases. if the linker should be excluded, provide the argument 'exclude'\n");
                            abort();
                        }
                    }
                }
                break;
            
            //make fasta file
            case 'f':
                make_fasta = 1;
                break;
                
            default: printf("error: unrecognized option. Aborting program...\n"); abort();
        }
    }
    /* ******* end of options parsing ******* */
    
    /* ******* make output directory ******* */
    char out_dir[MAX_LINE] = {0}; //array to store output directory name
    
    //construct output directory name
    ret = snprintf(out_dir, MAX_LINE, "%s_out", nm.vTmp);
    if (ret >= MAX_LINE || ret < 0) {
        printf("variant_maker: error - error when constructing output directory name. aborting...\n");
        abort();
    }
    mk_out_dir(out_dir); //make output directory
    
    //generate input details file
    ret = snprintf(prcs_out_nm, MAX_LINE, "./%s/%s_processing.txt", out_dir, nm.vTmp);
    if (ret >= MAX_LINE || ret < 0) {
        printf("variant_maker: error - error when constructing processing messages file name. aborting...\n");
        abort();
    }
    
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
    
    print_output(&nm, &bmap[0], vTmpCnt, varCnt, out_dir, append_priming, append_barcode, fp_brcd, cstm_lnkr, make_fasta);
    sprintf(out_msg, "\n%llu variants were generated from %d variant template(s)\n", (long long unsigned int)(v_indx-vTmpCnt), vTmpCnt);
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
    FILE * out_fp = NULL;         //output file pointer
    char out_nm[MAX_LINE] = {0};  //output file name
    
    int ret = 0; //variable for storing snprintf return value
    
    //construct input details file name
    ret = snprintf(out_nm, MAX_LINE, "./%s/%s_input.txt", out_dir, nm->vTmp);
    if (ret >= MAX_LINE || ret < 0) {
        printf("check_input: error - error when constructing input details file name. aborting...\n");
        abort();
    }
    
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


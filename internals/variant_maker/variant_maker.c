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

#include "../utils/debug.h"
#include "../utils/io_management.h"
#include "../utils/gen_utils.h"
#include "../seq_utils/ispair.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/t2u.h" //may not need this
#include "../seq_utils/basemap.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"
#include "./constant_seqs.h"
#include "./vmt_suffix.h"
#include "./read_varFile.h"
#include "./count_variants.h"
#include "./expand_variant_template.h"
#include "./filter_variants.h"

#include "make_barcodes.h"
#include "read_bcFile.h"
#include "print_output.h"

/* check_input: check that correct input files were supplied for
 variant generation and print input file information to file*/
void check_input(names * nm, int varFile_supplied, int brcdFile_supplied, int lib_type, int lib_type_set, int append_priming, int append_barcode, int first_bc_2_use, int bcs_per_var, char * cstm_lnkr, int make_fasta, char * out_dir);

//variant storage
fasta *vrnts = NULL;  //pointer to fasta structures used when generating variants
uint64_t v_indx = 0;  //index of current variant

FILE * prcs_ofp = NULL;           //output file pointer for processing messages
char prcs_out_nm[MAX_LINE] = {0}; //name of processing message output file
char out_msg[MAX_LINE] = {0};     //output message

char * fwd2use = NULL;
char * rev2use = NULL;

/* end of global variables */

int main(int argc, char *argv[])
{
    extern int debug; //flag for running debug mode
    
    extern fasta *vrnts;       //pointer to fasta structures used when generating variants
    extern uint64_t v_indx;    //index of current variant
    
    extern FILE * prcs_ofp;            //output file pointer for processing messages
    extern char prcs_out_nm[MAX_LINE]; //name of processing message output file
    extern char out_msg[MAX_LINE];     //output message
    
    extern char pra1_sc1[41];
    extern char c3sc1[34];
    extern char vra3[27];
    extern char RLA29synch_3p11[34]; //v2
    
    extern char vmt_suffix[4]; //variant maker target file suffix
    
    FILE * fp_var = NULL;      //pointer to variant parameters file
    FILE * fp_brcd = NULL;     //pointer to barcode sequences file
    
    struct names nm = {{0}};   //input names storage
    
    int mode = 0;              //run mode
    int varFile_supplied = 0;  //flag that variant template file was provided
    int brcdFile_supplied = 0; //flag that barcode file was provided
    int brcds2mk = 0;          //number of barcodes to make
    int append_barcode = 0;    //flag to append barcodes to variants
    int first_bc_2_use = 1;    //first barcode id to use. default = 1
    int append_priming = 0;    //flag to append C3-SC1 and VRA3 priming sites
    int bcs_per_var = 1;       //number of barcodes to apply to each variant
    int make_fasta = 0;        //flag to make fasta file
    int lib_type = -1;         //flag to indicate library type
    int lib_type_set = 0;      //flag that library type was set
        
    char usr_resp[4] = {0};         //storage for user response
    char discard[MAX_LINE+1] = {0}; //array for flushing stdin
    int resp_provided = 0;          //flag that valid reponse was provided
        
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int ret = 0; //variable for storing snprintf return value
    
    char cstm_lnkr[MAX_LINKER+1] = {"exclude"}; //array for storing custom linker sequence
    char *lnkr = NULL; //pointer to linker that will be used during library generation
    
    /* ******* parse options using getopt_long ******* */
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"mk-variants",     required_argument,  0,  'v'},  //variant template input file, sets MAKE_VARIANTS mode
            /* WARNING: mk-barcodes mode should only be run on systems with >= 128 GB of RAM */
            {"mk-barcodes",     required_argument,  0,  'b'},  //flag to make barcode file
            {"library-type",    required_argument,  0,  't'},  //type of library to be generated
            {"append-priming",  no_argument,        0,  'p'},  //append priming sites
            {"append-barcode",  required_argument,  0,  'a'},  //flag to append barcodes to variants
            {"first-bc-2-use",  required_argument,  0,  'i'},  //first barcode id to use
            {"brcds-per-vrnt",  required_argument,  0,  'c'},  //number of barcodes to apply to each variant
            {"custom-linker",   required_argument,  0,  'l'},  //use custom linker
            {"make_fasta",      no_argument,        0,  'f'},  //make fasta file
            {"debug",           no_argument,        0,  'd'},  //run debug mode
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "v:b:t:pa:i:c:l:fd", long_options, &option_index);
        
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
                
                printf("\nWARNING: barcode generation should only be performed on machines with >=128GB of RAM.\n\nproceed with barcode generation (yes/no)? \n");
                
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
                    printf("variant_maker: error - the maximum number of barcodes to make is currently %d. aborting...\n", MAX_BRCDS_2_MK);
                    abort();
                }
                
                break;
                
            case 't':
                if (!strcmp(argv[optind-1], "TECdisplay")) {          //making TECdisplay library
                    fwd2use = &c3sc1[0];
                    rev2use = &vra3[0];
                    lnkr = &cstm_lnkr[0];
                    lib_type = TECDISPLAY_LIB;
                    lib_type_set = 1;
                    
                } else if (!strcmp(argv[optind-1], "TECprobe-MUX")) { //making TECprobe-MUX library
                    fwd2use = &pra1_sc1[0];
                    rev2use = &vra3[0];
                    lnkr = &RLA29synch_3p11[0];
                    lib_type = TECPROBE_MUX_LIB;
                    lib_type_set = 1;
                    
                } else {
                    printf("variant_maker: ERROR - unrecognized library type. aborting...\n");
                    abort();
                }
                break;
                
            //append priming sites to variants in fasta file
            case 'p':
                if (append_priming == 1) {
                    printf("variant_maker: error - append_priming option can only be supplied once. aborting...\n");
                    abort();
                } else {
                    append_priming = 1;
                }
                break;
            
            //append pre-generated barcodes to variants
            case 'a':
                append_barcode = 1;                         //set flag to append barcode,
                get_file(&fp_brcd, argv[optind-1]);         //open barcode file
                strcpy(nm.iptB, argv[optind-1]);            //copy barcode file name to names struct
                get_sample_name(argv[optind-1], nm.brcd);   //get barcode file name
                brcdFile_supplied++;                        //increment number of barcode files supplied
                break;
            
            //set first barcode id to use
            case 'i':
                for (i = 0; argv[optind-1][i]; i++) {
                    if (!isdigit(argv[optind-1][i])) {
                        printf("variant_maker: ERROR - argument for first barcode to use must be composed of digits. aborting...\n");
                        abort();
                    }
                }
                first_bc_2_use = atoi(argv[optind-1]);
                break;
                
            //set number of barcodes per variant
            case 'c':
                for (i = 0; argv[optind-1][i]; i++) {
                    if (!isdigit(argv[optind-1][i])) {
                        printf("variant_maker: ERROR - argument for number of barcodes per variant must be composed of digits. aborting...\n");
                        abort();
                    }
                }
                
                bcs_per_var = atoi(argv[optind-1]); //set barcodes per variant
                
                if (!bcs_per_var) { //check that bcs_per_var is not set to 0
                    printf("variant_maker: error - number of barcodes per variant cannot be set to 0. aborting...");
                    abort();
                    
                } else if (bcs_per_var > MAX_BRCDS_PER_VAR) { //check that bcs_per_var is not > max
                    printf("variant_maker: error - number of barcodes per variant (%d) exceeds maximum (%d). aborting...\n", bcs_per_var, MAX_BRCDS_PER_VAR);
                    abort();
                }
                
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
                
            case 'd':
                printf("turning debug mode on\n");
                debug = 1;
                break;
                
            default: printf("error: unrecognized option. Aborting program...\n"); abort();
        }
    }
    /* ******* end of options parsing ******* */
    
    if (mode == MAKE_BARCODES) {
        mk_brcds(brcds2mk); //generate barcodes
        return 1;           //end program after barcodes are made
    }
    
    /* ******* make output directory ******* */
    char out_dir[MAX_LINE] = {0}; //array to store output directory name
    
    //construct output directory name
    ret = snprintf(out_dir, MAX_LINE, "%s_out", nm.vTmp);
    if (ret >= MAX_LINE || ret < 0) {
        printf("variant_maker: error - error when constructing output directory name. aborting...\n");
        abort();
    }
    mk_out_dir(out_dir); //make output directory
    
    //check and print input files
    check_input(&nm, varFile_supplied, brcdFile_supplied, lib_type, lib_type_set, append_priming, append_barcode, first_bc_2_use, bcs_per_var, cstm_lnkr, make_fasta, out_dir);
    
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
        
    print_output(&nm, &bmap[0], vTmpCnt, varCnt, out_dir, append_priming, append_barcode, fp_brcd, first_bc_2_use, bcs_per_var, lnkr, make_fasta, lib_type);
    
    if (append_barcode) {
        sprintf(out_msg, "\n%llu variants (%llu variants, %d barcode(s) per variant) were generated from %d variant template(s)\n", (long long unsigned int)((v_indx-vTmpCnt) * bcs_per_var), (long long unsigned int)(v_indx-vTmpCnt), bcs_per_var, vTmpCnt);
    } else {
        sprintf(out_msg, "\n%llu variants were generated from %d variant template(s)\n", (long long unsigned int)(v_indx-vTmpCnt), vTmpCnt);
    }
    
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
void check_input(names * nm, int varFile_supplied, int brcdFile_supplied, int lib_type, int lib_type_set, int append_priming, int append_barcode, int first_bc_2_use, int bcs_per_var, char * cstm_lnkr, int make_fasta, char * out_dir)
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
    
    //check that library type was defined if either append_priming or append_barcode options were set
    if ((append_priming || append_barcode) && !lib_type_set) {
        //TODO: expand library types listed in error message when new methods are ready to be made public
        printf("check_input: error - option to append priming sites and/or barcode was provided, but library type was not defined. define library type as TECdisplay using the --library-type/-t option\n");
        abort();
    }
    
    //check that make_fasta option was set if priming sites are to be appended
    if (append_priming && !make_fasta) {
        printf("check_input: error - priming sites are only appended to sequences within the fasta output file, and the 'make fasta' option (--make-fasta/-a) was not turned on. aborting...\n");
        abort();
    }
    
    //check that one or fewer barcodes file were supplied
    if (brcdFile_supplied == 1) {
        fprintf(out_fp, "barcodes: %s\n", nm->iptB);
    } else if (brcdFile_supplied > 1) {
        printf("check_input: error - incorrect number of barcode files (%d). expected 1. aborting...\n", brcdFile_supplied);
        abort();
    }
    
    //check that barcode file was provided alongside first barcode id to use option and that
    if (first_bc_2_use > 1) {
        fprintf(out_fp, "first BC id: %d\n", first_bc_2_use);
        
        if (!brcdFile_supplied) {
            printf("check_input: error - first barcode to use was set but no barcode file was provided. aborting...\n");
            abort();
        }
    //check that first barcode id to use was not set to zero
    } else if (!first_bc_2_use) {
        printf("check_input: error - the first barcode id to use cannot be 0. aborting...\n");
        abort();
    }
    
    //check that bcs_per_var is set to 1 if append barcode option was not turned on
    if (!append_barcode && bcs_per_var != 1) {
        printf("check_input: error - barcodes per variant value was set but append_barcode option was not provided. aborting...\n");
        abort();
    }
    
    
    //if lib_type is TECprobe-MUX, check that custom linker was not provided
    if ((lib_type == TECPROBE_MUX_LIB) && strcmp(cstm_lnkr, "exclude")) {
        printf("check_input: error - custom linker sequences cannot be used in TECprobe-MUX libraries. aborting...\n");
        abort();
    }
    
    //close input details file
    if (fclose(out_fp) == EOF) {
        printf("check_input: error - error occurred when closing input details file. Aborting program...\n");
        abort();
    }
    
    return;
}


//
//  filter_by_characteristics.c
//  
//
//  Created by Eric Strobel on 11/11/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <getopt.h>
#include <unistd.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"
#include "../utils/io_management.h"
#include "../utils/debug.h"
#include "../seq_utils/sto_suffix.h"
#include "../variant_maker/vmt_suffix.h"

#include "./filter_by_characteristics_defs.h"
#include "./filter_by_characteristics_structs.h"

#include "./set_seq_attributes/parse_seq_ipt.h"
#include "./set_seq_attributes/store_TECdisplay_data.h"
#include "./set_seq_attributes/init_seq_attributes.h"
#include "./set_seq_attributes/parse_descriptor_file.h"
#include "./set_seq_attributes/set_attributes.h"
#include "./print_output/print_annotated_data.h"
#include "./compare_variants/parse_comparison_file.h"
#include "./compare_variants/cmpr_smlr_seqs.h"
#include "./compare_variants/cmpr_smlr_structs.h"

char path2RNAStructure[MAX_LINE+1] = {0}; //path to RNAStructure Fold executable
void check_input(int msa_path_provided, int des_path_provided, int tecd_path_provided, int path2RNAStructure_provided);

int main(int argc, char *argv[])
{
    extern char path2RNAStructure[MAX_LINE+1]; //path to RNAStructure directory
    extern char sto_suffix[4];                 //stockholm format multiple sequence alignment file suffix
    extern char vmt_suffix[4];                 //variant maker targets multiple sequence alignment file suffix
    
    int path2RNAStructure_provided = 0; //tracks RNAstructure Fold paths provided
            
    char msa_path[MAX_LINE+1] = {0}; //path to multiple sequence alignment
    char des_path[MAX_LINE+1] = {0}; //path to descriptor file
    
    char msa_smpl_nm[MAX_LINE+1] = {0}; //multiple sequence alignment sample name
    char des_smpl_nm[MAX_LINE+1] = {0}; //descriptor file sample name
    
    char * msa_sffx = NULL; // pointer to multiple sequence alignment file suffix
    
    char tecd_path[MAX_LINE+1] = {0}; //path to TECdisplay data
    char tecp_path[MAX_LINE+1] = {0}; //path to TECprobe data

    char cmpv_path[MAX_LINE+1] = {0}; //path to comparison values file
    
    int msa_path_provided = 0; //tracks msa path options provided
    int des_path_provided = 0; //tracks descriptor path options provided
    
    int tecd_path_provided = 0; //tracks whether TECdisplay data was provided
    int tecp_path_provided = 0; //tracks whether TECprobe data was provided
    
    int cmpv_path_provided = 0; //tracks whether comparison values file was provided
    
    int msa_typ = FILE_TYPE_INIT; //multiple sequence alignment file type
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    /* ******* parse options using getopt_long ******* */
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"input",             required_argument,  0,  'i'}, //input sequences file
            {"TECdisplay",        required_argument,  0,  'D'}, //TECdisplay data
            {"comparison",        required_argument,  0,  'C'}, //comparison file
            {"descriptor",        required_argument,  0,  'd'}, //descriptor file
            {"RNAstructure-path", required_argument,  0,  'R'}, //RNAstructure fold path
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "i:D:C:d:R:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            
            //multiple sequence alignment input file
            case 'i':
                strcpy(msa_path, argv[optind-1]);             //store multiple sequence aligment file path
                get_sample_name(argv[optind-1], msa_smpl_nm); //get multiple sequence alignment sample name
                
                msa_sffx = &msa_path[strlen(msa_path)]; //set pointer to end of multiple sequence aligment file path string
                while ((uint64_t)(msa_sffx) != (uint64_t)(msa_path) && msa_sffx[-1] != '.') { //find '.' char ahead of file suffix
                    msa_sffx = &msa_sffx[-1];
                }

                if (msa_sffx[-1] == '.') { //if '.' char ahead of file suffix was found, set file type
                    if (!strcmp(msa_sffx, sto_suffix)) {
                        msa_typ = STO_FILE;
                    } else if (!strcmp(msa_sffx, vmt_suffix)) {
                        msa_typ = VMT_FILE;
                    } else {
                        printf("filter_by_characterstics: error - unexpected multiple sequence alignment format. expected stockholm (.%s) or variant maker targets (.%s). aborting...\n", sto_suffix, vmt_suffix);
                        abort();
                    }
                    
                } else {
                    printf("filter_by_characteristics: error - did not detect suffix for multiple sequence alignment input file. aborting...\n");
                    abort();
                }
                
                msa_path_provided++; //increment msa_path_provided flag (a check that only one multiple sequence alignment file
                break;               //was provided is performed in the 'check_input' function below)
            
            //TECdisplay data input file
            case 'D':
                strcpy(tecd_path, argv[optind-1]); //store path to TECdisplay data
                tecd_path_provided++;              //increment tecd_path_provided flag (a check that only one TECdisplay data file
                break;                             //was provided is performed in the 'check_input' function below)
            
            case 'C':
                strcpy(cmpv_path, argv[optind-1]); //store path to comparison values file
                cmpv_path_provided++;              //increment cmpv_path_provided flag (a check that only one comparison file
                break;                             //was provided is performed in the 'check_input' function below)
                
            //descriptor input file
            case 'd':
                strcpy(des_path, argv[optind-1]);             //store descriptor file path
                get_sample_name(argv[optind-1], des_smpl_nm); //get descriptor sample name
                des_path_provided++;                          //increment des_path_provided flag (a check that only one descriptor
                break;                                        //file was provided is performed in the 'check_input' function below)
            
            //path to RNA structure
            case 'R':
                strcpy(path2RNAStructure, argv[optind-1]); //store path to RNAStructure directory
                path2RNAStructure_provided++;              //increment path2RNAStructure_provided flag (a check that only one path
                break;                                     //was provided is performed in the 'check_input' function below)
                
            default: printf("error: unrecognized option. Aborting program...\n"); abort();
        }
    }
    /* ******* end of options parsing ******* */
    
    //check input files
    check_input(msa_path_provided, des_path_provided, tecd_path_provided, path2RNAStructure_provided);
    
    //process input files
    int seq_cnt = 0;                   //number of sequences
    int des_cnt = 0;                   //number of descriptors
    char * tecd_data_hdr = NULL;       //TECdisplay data header string
    char ** prsd_tecd_data_hdr = NULL; //parsed TECdisplay data header strings
    
    sequence_attributes * sq_att = {0};       //pointer for allocating sequence attributes structures
    seq_cnt = parse_seq_ipt(msa_path, msa_typ, &sq_att); //alloc sq_att mem and store input sequences
    
    //parse comparison file
    comparison_values * cmp = NULL; //pointer for allocating comparison values
    int cmp_cnt = 0; //number of comparison values in comparison file
    
    if (cmpv_path_provided) {
        cmp_cnt = parse_comparison_file(cmpv_path, &cmp);
    }
    
    //parse TECdisplay data file
    if (tecd_path_provided) {                                              //if TECdisplay data path was provided
        store_TECdisplay_data(sq_att, &tecd_data_hdr, &prsd_tecd_data_hdr, seq_cnt, tecd_path, cmp); //store TECdisplay data
    }
    
    //parse descriptor file
    descriptor * des = NULL;                       //pointer for allocating descriptor memory
    int typ_cnt[DSCRPTR_TYPES] = {0};              //variable for storing descriptor type counts
    des_cnt = parse_descriptor_file(des_path, &des, &typ_cnt[0]); //store descriptor information
    
    //generate output directory
    int ret = 0;                    //snprintf return value
    char out_dir[MAX_LINE+1] = {0}; //output directory name
    
    ret = snprintf(out_dir, MAX_LINE, "%s_%s_out", msa_smpl_nm, des_smpl_nm); //generate out dir name
    if (ret >= MAX_LINE || ret < 0) {
        printf("filter_by_characteristics: error - failed to generate output directory name. aborting...\n");
        abort();
    }
    
    mk_out_dir(out_dir); //make output directory
    chdir(out_dir);      //change to output directory
    
    //set sequence attributes
    descriptor_bank des_bnk = {0}; //descriptor structure bank
    
    init_seq_attributes(&sq_att, des, &des_bnk, seq_cnt, des_cnt, &typ_cnt[0]); //initialize seq_attributes struct
    set_attributes(sq_att, des, seq_cnt, des_cnt, path2RNAStructure);           //set attributes for each input sequence
    print_annotated_data(sq_att, seq_cnt, des_cnt, tecd_data_hdr);              //print output file containing annotated data
    
    //TODO: make running function conditional
    //TODO: change second arg to be variable input through comparison file
    for (i = 0; i < cmp_cnt; i++) {
        cmpr_smlr_structs(sq_att, des, seq_cnt, des_cnt, "IH1_38", &cmp[i]);
    }
    
    //cmpr_smlr_seqs(msa_path, msa_typ, sq_att, seq_cnt);
    
    //print number of structures predicted for each input descriptor that involves structure predictions. the number of
    //predictions can exceed the number of input sequences if a PRX_DG or DST_DG descriptor specifies that all predicted structures
    //should be kept. only one descriptor within a descriptor file can be set to keep all predicted structures
    for (i = 0; i < des_cnt; i++) {
        if (des[i].typ == PRX_DG || des[i].typ == DST_DG) {
            printf("%s\t%s\t%d structures predicted\n", des[i].typ_s, des[i].nm, des[i].tot_sp);
        }
    }
    printf("\n");
}

/* check_input: check validity of inputs */
void check_input(int msa_path_provided, int des_path_provided, int tecd_path_provided, int path2RNAStructure_provided)
{
    //check that exactly 1 multiple sequence alighment file was provided
    if (msa_path_provided != 1) {
        printf("check_input: expected one multiple sequence alignment file to be provided, detected %d. aborting...\n", msa_path_provided);
        abort();
    }
    
    //check that exactly 1 descriptor file was provided
    if (des_path_provided != 1) {
        printf("check_input: expected one descriptor file to be provided, detected %d. aborting...\n", des_path_provided);
        abort();
    }
    
    //check that <=1 TECdisplay data path was provided
    if (tecd_path_provided > 1) {
        printf("check input: expected one or fewer TECdisplay data files, detected %d. aborting...", tecd_path_provided);
        abort();
    }
    
    //check that exactly 1 path to RNAStructure Fold was provided
    if (path2RNAStructure_provided != 1) {
        printf("check_input: expected one path to RNAStructure Fold to be provided, detected %d. aborting...\n", path2RNAStructure_provided);
        abort();
    }
    
    return;
}

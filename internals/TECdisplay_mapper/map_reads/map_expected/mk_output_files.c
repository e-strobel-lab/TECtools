//
//  mk_output_files.c
//  
//
//  Created by Eric Strobel on 5/3/23.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../TECdisplay_mapper_structs.h"
#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_output_column_headers.h"

#include "../../../utils/io_management.h"
#include "../../../utils/gen_utils.h"

#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"
#include "../../../seq_utils/mapping_metrics.h"

#include "mk_output_files.h"

/* print_output_header: print data output file header line */
void print_output_header(FILE * out_fp, char * out_nm)
{
    extern const char TECdsply_clmn_hdrs[4][32]; //column headers from TECdisplay_output_column_headers.c
    
    int i = 0; //general purpose index
    
    //print header line. headers are printed in the order:
    //1. variant_id
    //2. bound read count
    //3. unbound read count
    //4. fraction bound
    for (i = 0; i < TDSPLY_HDR_CNT; i++) {                                  //for every column header
        if (i == TDSPLY_VID_HDR && i == 0) {                                //if at variant_id column
            if (out_nm[0]) {                                                //if an output name was supplied
                fprintf(out_fp, "%s_%s", out_nm, TECdsply_clmn_hdrs[i]);    //print "<out_nm>_id"
            } else {
                fprintf(out_fp, "variant_%s", TECdsply_clmn_hdrs[i]);       //otherwise, "print variant_id"
            }
            
        } else if (i == TDSPLY_BND_HDR || i == TDSPLY_UNB_HDR || i == TDSPLY_FRC_HDR) { //if printing data column
            if (out_nm[0]) {                                                //if an output name was supplied
                fprintf(out_fp, "\t%s_%s", out_nm, TECdsply_clmn_hdrs[i]);  //print "<tab><out_nm>_<clmn_hdr>"
            } else {
                fprintf(out_fp, "\t%s", TECdsply_clmn_hdrs[i]);             //otherwise, print "<tab><clmn_hdr>"
            }
            
        } else { // this should not be reachable                                        //error exceed max column
            printf("print_output_header: error - header index exceeded the maximum. aborting...\n");
            abort();
        }
    }
    
    fprintf(out_fp, "\n"); //print newline at the end of header line
}

/* print_output: print data output file */
void print_output(void * trgts, target_params * trg_prms, TDSPLY_names * nm, int mode)
{
    int i = 0; //general purpose index
    
    FILE *out_fp = {NULL};            //array for output file pointer
    char out_nm[256];                 //array for output name
    
    int ret = 0; //variable for storing snprintf return value
    
    //construct output file name
    if (nm->out_nm[0]) {
        ret = snprintf(out_nm, 256, "%s_out.txt", nm->out_nm);
        if (ret >= 256 || ret < 0) {
            printf("print_output: error - error when constructing output file name. aborting...\n");
            abort();
        }
    } else {
        sprintf(out_nm, "out.txt");
    }
    
    //generate output file
    if ((out_fp = fopen(out_nm, "w")) == NULL) {
        printf("print_output: error - could not open output file. Aborting program...\n");
        abort();
    }
        
    //print header line. headers are printed in the order:
    //1. variant_id
    //2. bound read count
    //3. unbound read count
    //4. fraction bound
    print_output_header(out_fp, nm->out_nm);
    
    //print data. as above for headers, data is printed in the order:
    //1. variant_id
    //2. bound read count
    //3. unbound read count
    //4. fraction bound
    for (i = 0; i < trg_prms->t_cnt; i++) {
        if (mode == STD_TDSPLY) {
            print_target_line(&(((target *)trgts)[i]), out_fp);
            
        } else if (mode == BRCD_TDSPLY) {
            print_compact_target_line(&(((compact_target *)trgts)[i]), out_fp);
            
        } else {
            printf("print_output: error - unrecognized mode. aborting...\n");
            abort();
        }
    }
    
    //close output file
    if (fclose(out_fp) == EOF) {
        printf("print_output: error - error occurred when closing output file. Aborting program...\n");
        abort();
    }
        
    return;
}

/* print_target_line: print target data line to output file*/
void print_target_line(target * trgt, FILE * ofp)
{
    extern const char TECdsply_clmn_hdrs[4][32]; //column headers from TECdisplay_output_column_headers.c
    
    int i = 0; //general purpose index
    
    opt_mx_trg * crnt_trg_val = NULL; //pointer for dereferencing target values
    opt_ref * crnt_ref_val = NULL;    //pointer for dereferencing reference values
    
    if (!trgt->mul) {                                        //if target is non-redundant
        crnt_trg_val = (opt_mx_trg *)(trgt->opt);            //set pointer to target values
        crnt_ref_val = (opt_ref *)(crnt_trg_val->ref->opt);  //set pointer to reference values
        
        for (i = 0; i < TDSPLY_HDR_CNT; i++) { //for every column header
            switch (i) {
                case TDSPLY_VID_HDR: //printing variant id column
                    fprintf(ofp, "%s", trgt->id);                   //print variant id
                    if (crnt_ref_val->cnstnts != NULL) {            //if ref target contains constant seq
                        fprintf(ofp, "_%s", crnt_ref_val->cnstnts); //print constants string
                    }
                    break;
                    
                case TDSPLY_BND_HDR: //printing bound read count
                    fprintf(ofp, "\t%d", crnt_trg_val->bnd);
                    break;
                    
                case TDSPLY_UNB_HDR: //printing unbound read count
                    fprintf(ofp, "\t%d", crnt_trg_val->unb);
                    break;
                    
                case TDSPLY_FRC_HDR: //printing fraction bound
                    fprintf(ofp, "\t%.4f",
                            (double)(crnt_trg_val->bnd)/(double)(crnt_trg_val->bnd + crnt_trg_val->unb));
                    break;
                    
                default: //this should not be reachable
                    printf("print_output: error - header index exceeded the maximum. aborting...\n");
                    abort();
                    break;
            }
        }
        fprintf(ofp, "\n"); //print newline at the end of data line
    }
}

/* print_compact_target_line: print compact_target data line to output file */
void print_compact_target_line(compact_target * ctrg, FILE * ofp)
{
    extern const char TECdsply_clmn_hdrs[4][32]; //column headers from TECdisplay_output_column_headers.c
    
    int i = 0; //general purpose index
    
    opt_BC * BC_val = NULL; //pointer for dereferencing optional barcode target values
    
    if (!ctrg->mul) {                   //if target is non-redundant
        BC_val = (opt_BC *)(ctrg->opt); //set pointer to target values
        
        for (i = 0; i < TDSPLY_HDR_CNT; i++) { //for every column header
            switch (i) {
                case TDSPLY_VID_HDR: //printing variant id column
                    fprintf(ofp, "%s", ctrg->cid);
                    break;
                    
                case TDSPLY_BND_HDR: //printing bound read count
                    fprintf(ofp, "\t%d", BC_val->chnl[BND]);
                    break;
                    
                case TDSPLY_UNB_HDR: //printing unbound read count
                    fprintf(ofp, "\t%d", BC_val->chnl[UNB]);
                    break;
                    
                case TDSPLY_FRC_HDR: //printing fraction bound
                    fprintf(ofp, "\t%.4f",
                            (double)(BC_val->chnl[BND])/(double)(BC_val->chnl[BND] + BC_val->chnl[UNB]));
                    break;
                    
                default: //this should not be reachable
                    printf("print_output: error - header index exceeded the maximum. aborting...\n");
                    abort();
                    break;
            }
        }
        fprintf(ofp, "\n"); //print newline at the end of data line
    }
}

/* print_metrics: print read mapping metrics */
void print_metrics(target_params * trg_prms, mapping_metrics * met, TDSPLY_names * nm)
{
    //generate output file
    FILE * out_fp = NULL; //output file pointer
    if ((out_fp = fopen("mapping.txt", "w")) == NULL) {
        printf("print_metrics: error - could not generate metrics file. Aborting program...\n");
        abort();
    }
    
    char out_str[MAX_LINE] = {0}; //array to store output strings
    
    int ret = 0; //variable for storing return value
    
    printf("\n\n");

    //print input files
    sprintf(out_str, "input files\n");
    printf2_scrn_n_fl(out_fp, out_str);
    
    ret = snprintf(out_str, MAX_LINE, "read 1:  %s\n", nm->file[READ1]);
    if (ret >= MAX_LINE || ret < 0) {
        printf("print_metrics: error - error when printing read 1 name to mapping record. aborting...\n");
        abort();
    }
    printf2_scrn_n_fl(out_fp, out_str);
    
    ret = snprintf(out_str, MAX_LINE, "read 2:  %s\n", nm->file[READ2]);
    if (ret >= MAX_LINE || ret < 0) {
        printf("print_metrics: error - error when printing read 2 name to mapping record. aborting...\n");
        abort();
    }
    printf2_scrn_n_fl(out_fp, out_str);
    
    ret = snprintf(out_str, MAX_LINE, "targets: %s\n", nm->trgs);
    if (ret >= MAX_LINE || ret < 0) {
        printf("print_metrics: error - error when printing targets name to mapping record. aborting...\n");
        abort();
    }
    printf2_scrn_n_fl(out_fp, out_str);
    
    //print merged sample name
    ret = snprintf(out_str, MAX_LINE, "\nmerged sample name:\n%s\n", nm->mrg);
    if (ret >= MAX_LINE || ret < 0) {
        printf("print_metrics: error - error when printing merged name string to mapping record. aborting...\n");
        abort();
    }
    printf2_scrn_n_fl(out_fp, out_str);
    
    if (nm->out_nm[0]) {
        ret = snprintf(out_str, MAX_LINE, "\noutput name:\n%s\n", nm->out_nm);
        if (ret >= MAX_LINE || ret < 0) {
            printf("print_metrics: error - error when printing output name to mapping record. aborting...\n");
            abort();
        }
        printf2_scrn_n_fl(out_fp, out_str);
    }
    
    //print overall read mapping metrics
    sprintf(out_str, "\n%d reads were assessed\n", met->reads_processed);
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "%d reads (%7.4f%%) mapped to a target\n", met->matches, 100*((float)(met->matches)/(float)(met->reads_processed)));
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "reads mapped to %d of %d (%6.2f%%) non-redundant targets\n", trg_prms->mapped2, trg_prms->nr_cnt, 100*((float)(trg_prms->mapped2)/(float)(trg_prms->nr_cnt)));
    printf2_scrn_n_fl(out_fp, out_str);
    
    //print channel distribution metrics
    sprintf(out_str, "\nthe channel distribution for %d mapped reads was:\n", met->matches);
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "  bound:    %9d reads (%5.2f%%)\n", met->chan_count[BND], 100*((float)(met->chan_count[BND])/(float)(met->matches)));
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "  unbound:  %9d reads (%5.2f%%)\n", met->chan_count[UNB], 100*((float)(met->chan_count[UNB])/(float)(met->matches)));
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "  unmapped: %9d reads (%5.2f%%)\n", met->chan_count[ERR], 100*((float)(met->chan_count[ERR])/(float)(met->matches)));
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "    of %9d bound reads,\n", met->chan_count[BND]);
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "      %9d (%5.2f%%) were identified with a full barcode match and\n", met->full_match[BND],
           100*((float)(met->full_match[BND])/(float)(met->chan_count[BND])));
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "      %9d (%5.2f%%) with a partial (3/4) barcode match\n", met->part_match[BND],
           100*((float)(met->part_match[BND])/(float)(met->chan_count[BND])));
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "    of %9d unbound reads,\n", met->chan_count[UNB]);
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "      %9d (%5.2f%%) were identified with a full barcode match and\n", met->full_match[UNB],
           100*((float)(met->full_match[UNB])/(float)(met->chan_count[UNB])));
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "      %9d (%5.2f%%) with a partial (3/4) barcode match\n", met->part_match[UNB],
           100*((float)(met->part_match[UNB])/(float)(met->chan_count[UNB])));
    printf2_scrn_n_fl(out_fp, out_str);
    
    //close metrics output file
    if ((fclose(out_fp)) == EOF) {
        printf("map_expected_reads: error - error occurred when closing mapping file. Aborting program...\n");
        abort();
    }
}

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
    for (i = 0; i < TDSPLY_HDR_CNT; i++) {                                              //for every column header
        if (i == TDSPLY_VID_HDR && i == 0) {                                            //if at variant_id column
            if (out_nm[0]) {                                                            //if an output name was supplied
                fprintf(out_fp, "%s_%s", out_nm, TECdsply_clmn_hdrs[i]);                //print "<out_nm>_id"
            } else {                                                                    //otherwise, "print variant_id"
                fprintf(out_fp, "variant_%s", TECdsply_clmn_hdrs[i]);
            }
            
        } else if (i == TDSPLY_BND_HDR || i == TDSPLY_UNB_HDR || i == TDSPLY_FRC_HDR) { //if printing data column
            fprintf(out_fp, "\t%s", TECdsply_clmn_hdrs[i]);                             //print <tab><clmn_hdr>
            
            
        } else { // this should not be reachable                                        //error exceed max column
            printf("print_output_header: error - header index exceeded the maximum. aborting...\n");
            abort();
        }
    }
    
    fprintf(out_fp, "\n"); //print newline at the end of header line
}

/* print_output: print data output file */
void print_output(target * trgts, target_params * trg_prms, names * nm)
{
    extern const char TECdsply_clmn_hdrs[4][32]; //column headers from TECdisplay_output_column_headers.c
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    FILE *out_fp = {NULL};            //array for output file pointer
    char out_nm[256];                 //array for output name
    opt_mx_trg * crnt_trg_val = NULL; //pointer for dereferencing target values
    opt_ref * crnt_ref_val = NULL;    //pointer for dereferencing reference values
    
    //construct output file name
    if (nm->out_nm[0]) {
        sprintf(out_nm, "%s_out.txt", nm->out_nm);
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
    for (i = 0; i < trg_prms->t_cnt; i++) {                    //for every target
        if (!trgts[i].mul) {                                   //if target is non-redundant
            crnt_trg_val = (opt_mx_trg *)trgts[i].opt;         //set pointer to target values
            crnt_ref_val = (opt_ref *)crnt_trg_val->ref->opt;  //set pointer to reference values
            
            for (j = 0; j < TDSPLY_HDR_CNT; j++) {                         //for every column header
                switch (j) {
                    case TDSPLY_VID_HDR: //printing variant id column
                        fprintf(out_fp, "%s", trgts[i].id);                //print variant id
                        if (crnt_ref_val->cnstnts != NULL) {               //if ref target contains constant seq
                            fprintf(out_fp, "_%s", crnt_ref_val->cnstnts); //print constants string
                        }
                        break;
                        
                    case TDSPLY_BND_HDR: //printing bount read count
                        fprintf(out_fp, "\t%d", crnt_trg_val->bnd);
                        break;
                        
                    case TDSPLY_UNB_HDR: //printing unbound read count
                        fprintf(out_fp, "\t%d", crnt_trg_val->unb);
                        break;
                        
                    case TDSPLY_FRC_HDR: //printing fraction bound
                        fprintf(out_fp, "\t%.4f",
                                (double)(crnt_trg_val->bnd)/(double)(crnt_trg_val->bnd + crnt_trg_val->unb));
                        break;
                        
                    default: //this should not be reachable
                        printf("print_output: error - header index exceeded the maximum. aborting...\n");
                        abort();
                        break;
                }
            }
            fprintf(out_fp, "\n"); //print newline at the end of data line
        }
    }
    
    //close output file
    if (fclose(out_fp) == EOF) {
        printf("print_output: error - error occurred when closing output file. Aborting program...\n");
        abort();
    }
        
    return;
}

/* print_metrics: print read mapping metrics */
void print_metrics(target * trgts, target_params * trg_prms, metrics * met, names * nm)
{
    //generate output file
    FILE * out_fp = NULL; //output file pointer
    if ((out_fp = fopen("mapping.txt", "w")) == NULL) {
        printf("print_metrics: error - could not generate metrics file. Aborting program...\n");
        abort();
    }
    
    char out_str[MAX_LINE] = {0}; //array to store output strings
    
    printf("\n\n");

    //print input files
    sprintf(out_str, "input files\n");
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "read 1:  %s\n", nm->file[READ1]);
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "read 2:  %s\n", nm->file[READ2]);
    printf2_scrn_n_fl(out_fp, out_str);
    
    sprintf(out_str, "targets: %s\n", nm->trgs);
    printf2_scrn_n_fl(out_fp, out_str);
    
    //print merged sample name
    sprintf(out_str, "\nmerged sample name:\n%s\n", nm->mrg);
    printf2_scrn_n_fl(out_fp, out_str);
    
    if (nm->out_nm[0]) {
        sprintf(out_str, "\noutput name:\n%s\n", nm->out_nm);
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

/* print_navigator_template: make template files for TECdisplay_navigator */
void print_navigator_template(target *refs, fasta * wt, target_params * trg_prms)
{
    FILE * out_fp;    //output file pointer
    char out_nm[256]; //array for output file name
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
    int l = 0; //general purpose index
    
    char tmp[BASE_FIELD+1] = {0}; //array to store base constraints for printing
    int found_end = 0;            //flag that the end of the variable base string was reached
    
    opt_ref * p_ref_val = NULL;   //pointer to for dereferencing reference values
    
    char navtmp_dir[22] = {"./navigator_templates"};

    mk_out_dir(navtmp_dir); //make directory for output files
    
    //generate a navigator template for every reference target
    for (i = 0; i < trg_prms->r_cnt; i++) {
        
        //construct navigator template name
        sprintf(out_nm, "%s/%s_navtmp.txt", navtmp_dir, refs[i].id);
        
        //generate output file
        if ((out_fp = fopen(out_nm, "w")) == NULL) {
            printf("print_output: error - could not open navigator template file. Aborting program...\n");
            abort();
        }
        
        p_ref_val = (opt_ref *)refs[i].opt; //dereference refs[i].opt to simplify code below
        
        fprintf(out_fp, "/seq\t%s\n", wt->sq);            //print wt sequence
        fprintf(out_fp, "/vbs\t%s\n", p_ref_val->ipt_sq); //print variable base template
        if (p_ref_val->cnstnts != NULL) {                 //print constant indels
            fprintf(out_fp, "/constant_indels:%s\n", p_ref_val->cnstnts);
        } else {
            fprintf(out_fp, "/constant_indels:none\n");
        }
        
        fprintf(out_fp, "<FILL IN filter name>\n");       //print message to fill in filter name
        
        //read vbases string and print each variable base as a single line
        //in the navigator template file. each base line begins with "base"
        //followed by a tab and right-aligned variable base information.
        for (found_end = 0, j = 0, k = 0; !found_end; j++) {
            if (p_ref_val->vbases[j] == '_' || !p_ref_val->vbases[j]) {
                if (k <= BASE_FIELD) { //if the BASE_FIELD(+1) width has not been exceeded
                    tmp[k] = '\0';     //terminate the tmp string
                    
                    //print the base entry to the navigator template file
                    fprintf(out_fp, "base\t");
                    for (l = 0; l < BASE_FIELD - k; l++) {
                        fprintf(out_fp, " ");
                    }
                    fprintf(out_fp, "%s\n", tmp);
                    k = 0;
                    
                    //check if the end of the vbases string was reached
                    if (!p_ref_val->vbases[j]) {
                        found_end = 1;
                    }
                    
                } else {
                    printf("print_navigator_template: error unexpected long (>%d char) base id in vbases line %s. aborting...", BASE_FIELD, p_ref_val->vbases);
                    abort();
                }
                
            } else {
                if (k < BASE_FIELD) {                //check that next char is within BASE_FIELD bounds
                    tmp[k++] = p_ref_val->vbases[j]; //copy variable base entry to tmp array
                } else {
                    printf("print_navigator_template: error unexpected long (>%d char) base id in vbases line %s. aborting...", BASE_FIELD, p_ref_val->vbases);
                    abort();
                }
            }
        }
        
        fprintf(out_fp, "<FILL IN any pair constraints>\n"); //print fill in pair constraints message
        fprintf(out_fp, "#\n"); //print constraint terminator and newline

        
        //close navigator template output file
        if ((fclose(out_fp)) == EOF) {
            printf("print_navigator template: error - error occurred when closing the navigator template file. Aborting program...\n");
            abort();
        }
    }
    
    return;
}

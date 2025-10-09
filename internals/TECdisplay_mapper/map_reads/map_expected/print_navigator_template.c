//
//  print_navigator_template.c
//  
//
//  Created by Eric Strobel on 10/9/25.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../TECdisplay_mapper_structs.h"
#include "../../TECdisplay_mapper_defs.h"

#include "../../../utils/io_management.h"
#include "../../../utils/gen_utils.h"

#include "print_navigator_template.h"

/* print_navigator_template: make template files for TECdisplay_navigator */
void print_navigator_template(target *refs, TDSPLY_fasta * wt, target_params * trg_prms)
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

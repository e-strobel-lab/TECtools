//
//  print_annotated_data.c
//  
//
//  Created by Eric Strobel on 3/5/26.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"

#include "../set_seq_attributes/set_attributes.h"

#include "print_annotated_data.h"

/* print_annotated_data: print output file that includes all sequences, associated data, and structure predictions */
void print_annotated_data(sequence_attributes * sq_att, int seq_cnt, int des_cnt, char * tecd_data_hdr)
{
    extern char lnkr[6]; //linker sequence used when predicting distal base pairing interactions
    
    FILE * ofp = NULL; //output file pointer
    
    //open output file
    if ((ofp = fopen("out.txt", "w")) == NULL) {
        printf("print_annotated_data: error - failed to open output file. aborting...");
        abort();
    }
    
    //print standard output file headers
    fprintf(ofp, "variant_id\tsequence\tmultiple_structs\tprofile_number");
    if (tecd_data_hdr != NULL) {
        fprintf(ofp, "\t%s", tecd_data_hdr);
    }
    
    int i = 0;  //general purpose index
    int j = 0;  //general purpose index
    int d = 0;  //descriptor index
    int sp = 0; //structProps index
    
    nucleotide_identity * p_nuc_id = NULL; //pointer to nucleotide_identity struct
    proximal_deltaG * p_prx_dG = NULL;     //pointer to proximal_deltaG struct
    distal_deltaG * p_dst_dG = NULL;       //pointer to distal_deltaG struct
    subsequence_length * p_ss_len = NULL;  //pointer to subsequence_length struct
    
    //print column headers for descriptor values
    for (d = 0; d < des_cnt; d++) {
        switch (sq_att[0].des[d].typ) {
            case NUC_ID:
                fprintf(ofp, "\tnuc_id:%s_seq", sq_att[0].des[d].nm);
                break;
                
            case PRX_DG:
                fprintf(ofp, "\tprx_dG:%s_seq", sq_att[0].des[d].nm);
                fprintf(ofp, "\tprx_dG:%s_dot", sq_att[0].des[d].nm);
                fprintf(ofp, "\tprx_dG:%s_dot_ann", sq_att[0].des[d].nm);
                fprintf(ofp, "\tprx_dG:%s_dG",  sq_att[0].des[d].nm);
                break;
                
            case DST_DG:
                fprintf(ofp, "\tdst_dG:%s_seq", sq_att[0].des[d].nm);
                fprintf(ofp, "\tdst_dG:%s_dot", sq_att[0].des[d].nm);
                fprintf(ofp, "\tdst_dG:%s_dot_ann", sq_att[0].des[d].nm);
                fprintf(ofp, "\tdst_dG:%s_dG",  sq_att[0].des[d].nm);
                break;
                
            case SS_LEN:
                fprintf(ofp, "\tss_len:%s_seq", sq_att[0].des[d].nm);
                fprintf(ofp, "\tss_len:%s_len", sq_att[0].des[d].nm);
                break;
                
            default:
                printf("print_annotated_data: error - unrecognized sequence attribute type. aborting...\n");
                abort();
                break;
        }
    }
    fprintf(ofp, "\n");
    
    //up to one descriptor within the descriptor file is permitted to keep alternative predicted structures
    //that approximate the minimum free energy structure. these alternative structures are accommodated by
    //printing the same variant multiple times, each with one of the alternative structures. variants with
    //multiple entries are indicated by the number in the 'profile' column
    
    int max_sp_cnt = 1;           //number of structProps associated with a variant
    int tmp_sp_cnt = 0;           //temp variable for number of structProps associated with a descriptor of a variant
    structProps * crnt_sp = NULL; //pointer to the current structProps struct
    
    for (i = 0; i < seq_cnt; i++) { //for every sequence
        
        //determine maximum sp_cnt
        for (d = 0, max_sp_cnt = 1; d < des_cnt; d++) { //for every descriptor
            
            tmp_sp_cnt = 0; //set tmp_sp_cnt to 0
            
            switch (sq_att[i].des[d].typ) {
                    
                case PRX_DG:
                    p_prx_dG = (proximal_deltaG *)(sq_att[i].att[d]); //point proximal_deltaG structure pointer to attribute
                    tmp_sp_cnt = p_prx_dG->sp_cnt;                    //set tmp_sp_cnt
                    break;
                    
                case DST_DG:
                    p_dst_dG = (distal_deltaG *)(sq_att[i].att[d]);   //point distal_deltaG structure pointer to attribute
                    tmp_sp_cnt = p_dst_dG->sp_cnt;                    //set_tmp_sp_cnt
                    break;
                    
                default:
                    break;
            }
            
            if (tmp_sp_cnt > max_sp_cnt) {
                max_sp_cnt = tmp_sp_cnt;
            }
        }
                
        for (sp = 0; sp < max_sp_cnt; sp++) { //for every structProps
            
            //print name, sequence, and profile number information
            fprintf(ofp, "%s\t%s\t%s\t%d", sq_att[i].nm, sq_att[i].sq1, (max_sp_cnt > 1) ? "Y" : "N", sp+1);
            
            if (tecd_data_hdr != NULL) {              //if TECdisplay data was provided
                fprintf(ofp, "\t%s", sq_att[i].tecd); //print the TECdisplay data
            }
            
            for (d = 0; d < des_cnt; d++) {  //for each descriptor
                switch (sq_att[i].des[d].typ) {
                    case NUC_ID:
                        p_nuc_id = (nucleotide_identity *)(sq_att[i].att[d]); //point nucleotide_identity structure pointer to attribute
                        fprintf(ofp, "\t%s", p_nuc_id->sq);                   //print sequence
                        break;
                        
                    case PRX_DG:
                        p_prx_dG = (proximal_deltaG *)(sq_att[i].att[d]); //point proximal_deltaG structure pointer to attribute
                        crnt_sp = &p_prx_dG->sp;                          //set pointer to structProps
                        
                        if (sp > 0 && p_prx_dG->sp_cnt > 1) {  //if not first iteration and there is more than one structProps
                            for (j = 0; j < sp; j++) {         //iterate through the structProps linked list
                                crnt_sp = crnt_sp->nxt;        //and set crnt_sp to the next structProps
                            }
                        }
                        
                        fprintf(ofp, "\t%s",   p_prx_dG->sq);    //print sequence used for prediction
                        fprintf(ofp, "\t%s",   crnt_sp->db);     //print dot-bracket structure
                        fprintf(ofp, "\t%s",   crnt_sp->db_an);  //print dot-bracket annotated base-pair structure
                        fprintf(ofp, "\t%.1f", crnt_sp->dG);     //print deltaG
                        break;
                        
                    case DST_DG:
                        p_dst_dG = (distal_deltaG *)(sq_att[i].att[d]); //point distal_deltaG structure pointer to attribute
                        crnt_sp = &p_dst_dG->sp;                        //set pointer to structProps
                        
                        if (sp > 0 && p_dst_dG->sp_cnt > 1) {  //if not first iteration and there is more than one structProps
                            for (j = 0; j < sp; j++) {         //iterate through the structProps linked list
                                crnt_sp = crnt_sp->nxt;        //and set crnt_SP to the next structProps
                            }
                        }
                        
                        fprintf(ofp, "\t%s%s%s", p_dst_dG->up, lnkr, p_dst_dG->dn); //print sequence used for prediction
                        fprintf(ofp, "\t%s",     p_dst_dG->sp.db);                  //print dot-bracket structure
                        fprintf(ofp, "\t%s",     p_dst_dG->sp.db_an);               //print dot-bracket structure
                        fprintf(ofp, "\t%.1f",   p_dst_dG->sp.dG);                  //print deltaG
                        
                        break;
                        
                    case SS_LEN:
                        p_ss_len = (subsequence_length *)(sq_att[i].att[d]); //point subsequence_length structure pointer to attribute
                        fprintf(ofp, "\t%s", p_ss_len->sq);                  //print subsequence
                        fprintf(ofp, "\t%d", p_ss_len->len);                 //print length
                        break;
                        
                    default:
                        printf("print_annotated_data: error - unrecognized sequence attribute type. aborting...\n");
                        break;
                }
            }
            fprintf(ofp, "\n");
        }
    }
    
    //close file
    if (fclose(ofp) == EOF) {
        printf("print_annotated_data: error - failed to open output file. aborting...\n");
        abort();
    }
}

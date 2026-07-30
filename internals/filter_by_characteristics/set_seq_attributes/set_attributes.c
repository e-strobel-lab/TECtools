//
//  set_attributes.c
//  
//
//  Created by Eric Strobel on 12/1/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"
#include "../../seq_utils/mk_fasta.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"

#include "./parse_descriptor_file.h"
#include "./get_msa_subseq.h"
#include "./run_RNAStructure.h"
#include "./parse_ct_file.h"

#include "set_attributes.h"

char tmp_fasta_path[11] = "tmp_fasta/";
char ct_output_path[11] = "ct_output/";

char lnkr[6] = "aaaaa";        //linker for joining distal subsequences during RNA structure prediction

/* set_attributes: manages sequence attribute setting */
void set_attributes(sequence_attributes * sq_att, descriptor * des, int seq_cnt, int des_cnt, char * path2RNAStructure)
{
    //TODO: consider only generating if RNA structure will be called
    mk_out_dir(tmp_fasta_path); //make output directory for temporary fasta files
    mk_out_dir(ct_output_path); //make output directory for ct files
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    //set attribute values for each input sequence
    for (i = 0; i < seq_cnt; i++) {
        for (j = 0; j < des_cnt; j++) {
            switch (sq_att->des[j].typ) {
                case NUC_ID:
                    set_nuc_id((nucleotide_identity *)(sq_att[i].att[j]), &sq_att[i].des[j], &sq_att[i]);
                    break;
                    
                case PRX_DG:
                    set_prx_dG((proximal_deltaG *)(sq_att[i].att[j]), &sq_att[i].des[j], &sq_att[i], path2RNAStructure);
                    break;
                    
                case DST_DG:
                    set_dst_dG((distal_deltaG *)(sq_att[i].att[j]), &sq_att[i].des[j], &sq_att[i], path2RNAStructure);
                    break;
                    
                case SS_LEN:
                    set_ss_len((subsequence_length *)(sq_att[i].att[j]), &sq_att[i].des[j], &sq_att[i]);
                    break;
                    
                default:
                    printf("set_attributes: error - unrecognized sequence attribute type. aborting...\n");
                    abort();
                    break;
            }
        }
    }
}

/* set_nuc_id: set nucleotide_identity structure values */
void set_nuc_id(nucleotide_identity * nuc_id, descriptor * des, sequence_attributes * sq_att)
{
    if (des->typ != NUC_ID) { //throw error if descriptor type is not nucleotide identity
        printf("set_nuc_id: error - descriptor is not nucleotide identity type. aborting...\n");
        abort();
    }
    
    nuc_id->sq_att = sq_att; //set pointer to associated sequence attributes structure
    
    //store multiple sequence alignment subsequence in nucleotide_identity structure
    get_msa_subseq(&nuc_id->sq, sq_att->sq0, des->wndw[0].b1, des->wndw[0].b2, BOUND2_LENGTH);
    //printf("%s\n", nuc_id->sq);
    
    return;
}

/* set_prx_dG: set proximal_deltaG structure values */
void set_prx_dG(proximal_deltaG * prx_dG, descriptor * des, sequence_attributes * sq_att, char * path2RNAStructure)
{
    extern char path2fold[MAX_LINE+1];
    extern char tmp_fasta_path[11];
    extern char ct_output_path[11];
    extern char fa_sffx[4];
    
    char command[MAX_LINE+1] = {0}; //array for storing system commands
    
    if (des->typ != PRX_DG) { //throw error if descriptor type is not proximal deltaG
        printf("set_prx_dG: error - descriptor is not proximal deltaG type. aborting...\n");
        abort();
    }
    
    prx_dG->sq_att = sq_att; //set pointer to associated sequence attributes structure
        
    //store multiple sequence alignment subsequence in proximal deltaG structure
    get_msa_subseq(&prx_dG->sq, sq_att->sq0, des->wndw[0].b1, des->wndw[0].b2, BOUND2_INDEX);
        
    //generate fasta file that will be used for RNA structure prediction
    char fa_nm[MAX_LINE+1] = {0}; //array to store fasta name
    int ret = 0;                  //return value of snprintf
    ret = snprintf(fa_nm, MAX_LINE, "%s_%s", sq_att->nm, des->nm);
    if (ret >= MAX_LINE || ret < 0) {
        printf("set_prx_dG: error - failed to generate fasta name. aborting...\n");
        abort();
    }
    mk_fasta_file(fa_nm, prx_dG->sq, tmp_fasta_path);
        
    //run structure prediction using RNAstructure Fold algorithm
    char path2ct[MAX_LINE+1]; //path to connectivity table file
    run_RNAStructure_Fold(path2RNAStructure, fa_nm, tmp_fasta_path, ct_output_path, fa_sffx, path2ct, MAX_LINE, des->act.ptyp);
    
    con_table ct = {0}; //connectivity table structure
    
    prx_dG->sp_cnt = parse_ct_file(&ct, path2ct, fa_nm); //parse connectivity table and set sp_count for attribute
    des->tot_sp += prx_dG->sp_cnt;                       //increment total number of structProps for current descriptor
    
    set_structProps(sq_att, prx_dG, PRX_DG, &prx_dG->sp, &ct, prx_dG->sp_cnt);  //set structProps values
    free_con_table_mem(&ct, prx_dG->sp_cnt);                                    //free allocated con_table memory
    
    snprintf(command, MAX_LINE, "rm %s* %s*", tmp_fasta_path,  ct_output_path); //make command for removing temp files
    system(command);                                                            //remove temp files
    
    return;
}

/* set_dst_dG: set distal_deltaG structure values */
void set_dst_dG(distal_deltaG * dst_dG, descriptor * des, sequence_attributes * sq_att, char * path2RNAStructure)
{
    extern char path2fold[MAX_LINE+1];
    extern char tmp_fasta_path[11];
    extern char ct_output_path[11];
    extern char lnkr[6];
    extern char fa_sffx[4];
    
    char cat_sq[MAX_LINE+1] = {0};  //array to store concatenated subsequences
    
    char command[MAX_LINE+1] = {0}; //array for storing system commands
    
    if (des->typ != DST_DG) { //throw error if descriptor type is not distal deltaG
        printf("set_dst_dG: error - descriptor is not distal deltaG type. aborting...\n");
        abort();
    }
    
    dst_dG->sq_att = sq_att; //set pointer to associated sequence attributes structure
    
    //store multiple sequence alignment subsequences in distal deltaG structure
    get_msa_subseq(&dst_dG->up, sq_att->sq0, des->wndw[0].b1, des->wndw[0].b2, BOUND2_INDEX);
    get_msa_subseq(&dst_dG->dn, sq_att->sq0, des->wndw[1].b1, des->wndw[1].b2, BOUND2_INDEX);
    
    //printf("%s\n", dst_dG->up);
    //printf("%s\n", dst_dG->dn);
    
    //generate sequence in which distal subsequences are joined by an unstructured linker,
    //which will be used for RNA structure prediction
    int ret = 0; //return value of snprintf
    ret = snprintf(cat_sq, MAX_LINE, "%s%s%s", dst_dG->up, lnkr, dst_dG->dn);
    if (ret >= MAX_LINE || ret < 0) {
        printf("set_dst_dG: error - failed to construct sequence for secondary structure prediction. aborting...\n");
        abort();
    }
    
    //generate fasta file that will be used for RNA structure prediction
    char fa_nm[MAX_LINE+1] = {0};
    ret = snprintf(fa_nm, MAX_LINE, "%s_%s", sq_att->nm, des->nm);
    if (ret >= MAX_LINE || ret < 0) {
        printf("set_prx_dG: error - failed to generate fasta name. aborting...\n");
        abort();
    }
    mk_fasta_file(fa_nm, cat_sq, tmp_fasta_path);
        
    //run structure prediction using RNAstructure Fold algorithm
    char path2ct[MAX_LINE+1]; //path to connectivity table file
    run_RNAStructure_Fold(path2RNAStructure, fa_nm, tmp_fasta_path, ct_output_path, fa_sffx, path2ct, MAX_LINE, des->act.ptyp);
    
    con_table ct = {0}; //connectivity table structure
    
    dst_dG->sp_cnt = parse_ct_file(&ct, path2ct, fa_nm); //parse connectivity table and set sp_count for attribute
    des->tot_sp += dst_dG->sp_cnt;                       //increment total number of structProps for current descriptor
    
    set_structProps(sq_att, dst_dG, DST_DG, &dst_dG->sp, &ct, dst_dG->sp_cnt);  //set structProps values
    free_con_table_mem(&ct, dst_dG->sp_cnt);                                    //free allocated con_table memory
    
    snprintf(command, MAX_LINE, "rm %s* %s*", tmp_fasta_path,  ct_output_path); //make command for removing temp files
    system(command);                                                            //remove temp files
        
    return;
}

/* set_ss_len: set subsequence_length structure values */
void set_ss_len(subsequence_length * ss_len, descriptor * des, sequence_attributes * sq_att)
{
    if (des->typ != SS_LEN) { //throw error if descriptor type is not subsequence length
        printf("set_ss_len: error - descriptor is not subsequence length type. aborting...\n");
        abort();
    }
    
    ss_len->sq_att = sq_att; //set pointer to associated sequence attributes structure
    
    //store multiple sequence alignment subsequence and subsequence length in subsequence length structure
    get_msa_subseq(&ss_len->sq, sq_att->sq0, des->wndw[0].b1, des->wndw[0].b2, BOUND2_INDEX);
    ss_len->len = strlen(ss_len->sq);
    
    //printf("%s\n", ss_len->sq);
    
    return;
}

/* set_structProps: set structProps values */
void set_structProps(sequence_attributes * sq_att, void * att, int att_typ, structProps * sp, con_table * ct, int ct_cnt)
{
    structProps * crnt_sp = sp; //pointer to current structProps structure
    con_table * crnt_ct = ct;   //pointer to current connectivity table
    
    int i = 0; //general purpose index
    
    for (i = 0; i < ct_cnt; i++) { //for every ct table that was generated
        
        if (i > 0) { //if not processing first ct tabel
            
            //allocate new structProps structure in linked list
            if ((crnt_sp->nxt = calloc(1, sizeof(*crnt_sp->nxt))) == NULL) {
                printf("set_structProps: error - failed to allocate structure properties memory");
                abort();
            }
            crnt_sp = crnt_sp->nxt; //point crnt_sp to new structProps structure
            crnt_ct = crnt_ct->nxt; //point crnt_ct to next connectivity table
        }
        
        crnt_sp->sq_att = sq_att;   //set pointer to parent sequence attributes structure
        crnt_sp->att = att;         //set pointer to parent attributes structure
        crnt_sp->att_typ = att_typ; //set attributes type
        
        crnt_sp->dG = crnt_ct->dG;  //copy deltaG
            
        //allocate memory for and copy dot bracket string
        if ((crnt_sp->db = malloc((strlen(ct->db)+1) * sizeof(*(crnt_sp->db)))) == NULL) {
            printf("set_structProps: error - failed to allocate memory for dot-bracket notation structure. aborting...\n");
            abort();
        }
        strcpy(crnt_sp->db, crnt_ct->db);
        
        //allocate memory for and copy annotated dot bracket string
        if ((crnt_sp->db_an = malloc((strlen(ct->db_an)+1) * sizeof(*(crnt_sp->db_an)))) == NULL) {
            printf("set_structProps: error - failed to allocate memory for annotated dot-bracket notation structure. aborting...\n");
            abort();
        }
        strcpy(crnt_sp->db_an, crnt_ct->db_an);
    }
    
    return;
}

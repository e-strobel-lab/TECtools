//
//  mk_TDSPLY_test_data.c
//  
//
//  Created by Eric Strobel on 4/26/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../variant_maker/variant_maker_defs.h"

#include "../TECdisplay_mapper_structs.h"
#include "../TECdisplay_mapper_defs.h"

#include "../map_reads/map_expected/parse_vmt_trgts.h"

#include "../../utils/debug.h"

#include "mk_TDSPLY_test_data.h"

extern int debug;

/* mk_TDSPLY_test_data: coordinates test data generation */
int mk_TDSPLY_test_data(TDSPLY_names * nm, target *refs, target *trgts, target_params *trg_prms)
{
    
    srand(time(NULL)); //seed pseudorandom number generator
    
    char out_dir[MAX_LINE] = {0}; //array to store output directory
    
    int ret = 0; //variable for storing snprintf output name
    
    //construct output directory name
    ret = snprintf(out_dir, MAX_LINE, "%s_test_data", nm->trgs_prefix);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_TDSPLY_test_data: error - error when constructing output directory name. aborting...\n");
        abort();
    }
    mk_out_dir(out_dir); //make output directory
    
    FILE * out_rd1 = NULL;    //read 1 file pointer
    FILE * out_rd2 = NULL;    //read 2 file pointer

    //generate read 1 output file
    ret = snprintf(nm->file[READ1], MAX_LINE, "./%s/%s_testdata_R1.fq", out_dir, nm->trgs_prefix);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_TDSPLY_test_data: error - error when constructing read one output file name. aborting...\n");
        abort();
    }
    if ((out_rd1 = fopen(nm->file[READ1], "w")) == NULL) {
        printf("mk_TDSPLY_test_data: error - could not generate test data read one file. Aborting program...\n");
        abort();
    }
    
    //generate read 2 output file
    ret = snprintf(nm->file[READ2], MAX_LINE, "./%s/%s_testdata_R2.fq", out_dir, nm->trgs_prefix);
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_TDSPLY_test_data: error - error when constructing read two output file name. aborting...\n");
        abort();
    }
    if ((out_rd2 = fopen(nm->file[READ2], "w")) == NULL) {
        printf("mk_TDSPLY_test_data: error - could not generate test data read two file. Aborting program...\n");
        abort();
    }
    
    //channel barcode generation variables
    static const int mk_chnl[CHNL_CLASSES] = {  UNB,  BND,  UNB,  BND,  ERR};  //channel type
    static const int mk_mtch[CHNL_CLASSES] = { FULL, FULL, PART, PART, FULL};  //channel match type
    
    opt_ref * crnt_ref_val = NULL;    //pointer to target values for current reference sequence
    
    char chnl_bc[5] = {0};            //array for storing channel barcode, append to head of read 1
    char vb_map[MAX_LINE] = {0};      //array to store map of variable base locations
    char mut_insrt[MAX_LINE] = {0};   //array to store mutated insert sequence

    int ref_indx = 0; //reference index
    int vb_indx  = 0; //variable base array index
    int tpr_indx = 0; //targets per reference index
    int trg_indx = 0; //target index
    int chnl_cd  = 0; //channel code
    int mut_cd = 0;   //mutation code
    int i = 0;        //general purpose index
        
    //generate test data reads. test data is only made using non-redundant targets.
    //five mappable reads and 100 unmappable reads are made for each non-redundant
    //target. the five mappable reads contain the following randomized channel barcodes:
    
    //  1. bound, full match
    //  2. bound, partial match
    //  3. unbound, full match
    //  4. unbound, partial match
    //  5. unmappable
    
    //the 100 unmappable reads contain a random 1 nt substitution, insertion, or deletion
    //and a randomly selected, randomized, mappable channel barcode
     
    for (ref_indx = 0, trg_indx = 0; ref_indx < trg_prms->r_cnt; ref_indx++) { //for every reference sequence
        
        crnt_ref_val = (opt_ref *)refs[ref_indx].opt; //dereference reference target optional values
        
        //make variable base map. the variable base map marks the location of
        //variable bases in the current reference sequence. it is used below
        //when generating unmappable reads to avoid making changes at variable
        //bases that would map to another target
        for (i = 0, vb_indx = 0; refs[ref_indx].sq[i]; i++) {
            if (i == crnt_ref_val->vb_pos[vb_indx] && vb_indx <crnt_ref_val->vb_cnt) { //index is for variable base
                vb_map[i] = '*'; //set variable base locations to '*'
                vb_indx++;       //increment variable base index
            } else {             //index is for a non-variable base
                vb_map[i] = '.'; //set all other locations to '.'
            }
        }
        vb_map[i] = '\0'; //terminate variable base map arra
        
        //printf("%s\n%s\n", refs[ref_indx].sq, vb_map);
        
        //for every non-redundant target associated with the current reference target
        for (tpr_indx = 0; tpr_indx < trg_prms->t_per_r[ref_indx]; tpr_indx++, trg_indx++) {
            if (!trgts[trg_indx].mul) { //target is not redundant
                
                //make native sequence reads with variable channel codes
                for (i = 0; i < CHNL_CLASSES; i++) {
                    mk_rndmzd_TDSPLY_bc(chnl_bc, mk_chnl[i], mk_mtch[i]); //generate randomized channel barcode
                    
                    //make reads from native sequence
                    print_TDSPLY_fq(out_rd1, out_rd2, &trgts[trg_indx].id[0], &trgts[trg_indx].sq[0], chnl_bc, NAT);
                }
                
                //make mutant sequence reads
                for (i = 0; i < 100; i++) {                      //make 100 mutants for every target
                    chnl_cd = rand() % 2;                        //select random mappable channel barcode
                    mk_rndmzd_TDSPLY_bc(chnl_bc, chnl_cd, FULL); //generate randomized channel barcode
                    strcpy(mut_insrt, trgts[trg_indx].sq);       //copy target sequence to mut array
                    mut_cd = mutate_insrt(vb_map, mut_insrt);    //generate mutant target sequence
                    
                    //make reads from mutated sequence
                    print_TDSPLY_fq(out_rd1, out_rd2, &trgts[trg_indx].id[0], &mut_insrt[0], chnl_bc, mut_cd);
                }
            }
        }
    }
    
    //close read 1 output file
    if ((fclose(out_rd1)) == EOF) {
        printf("mk_TDSPLY_test_data: error - error occurred when closing test data read 1 file. Aborting program...\n");
        abort();
    }
    
    //close read 2 output file
    if ((fclose(out_rd2)) == EOF) {
        printf("mk_TDSPLY_test_data: error - error occurred when closing test data read 2 file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

/* mk_rndmzd_TDSPLY_bc: generate a randomized channel barcode with variable channel and match settings.
 chnl variable:
 BND = bound barcode (RYYY)
 UNB = unbound barcode  (YRRR)
 ERR = undetermined barcode (<3/4 position match to expected sequence)
 
 mtch variable:
 FULL = 4/4 position match to expected sequence
 PART = 3/4 position match to expected sequence
 */
void mk_rndmzd_TDSPLY_bc(char * bc, int chnl, int mtch)
{
    int i = 0;  //general purpose index
    
    static const char pur[3] = "AG";    //array for randomizing purine (R) bases
    static const char pyr[3] = "TC";    //array for randomizing pyrimidine (Y) bases
    static const char all[5] = "ATGC";  //array for general randomization (N bases)
    
    int bnd_cnt = 0;    //counter for matches to bound barcode
    int unb_cnt = 0;    //counter for matches to unbound barcode
    
    int bc_made = 0;    //flag to indicate that random barcode meets chnl and mtch specifications
    
    if (chnl == BND && mtch == FULL) {        //make bound full match barcode
        bc[0] = pur[rand() % 2];
        bc[1] = pyr[rand() % 2];              //simple randomization of each base such that
        bc[2] = pyr[rand() % 2];              //the resulting sequence matches RYYY
        bc[3] = pyr[rand() % 2];
        bc[4] = '\0';
    } else if (chnl == UNB && mtch == FULL) { //make unbound full match barcode
        bc[0] = pyr[rand() % 2];
        bc[1] = pur[rand() % 2];              //simple randomization of each base such that
        bc[2] = pur[rand() % 2];              //the resulting sequences matches YRRR
        bc[3] = pur[rand() % 2];
        bc[4] = '\0';
    } else { //procedure for making unmappable and partial match barcodes
        bc_made = 0;
        while (!bc_made) { //repeat until random barcode meets chnl and mtch specifications
            bnd_cnt = 0;   //initialize bound barcode base count to zero
            unb_cnt = 0;   //initialize unbound barcode base count to zero
            for (i = 0; i < 4; i++) {    //generate a random 4-mer
                bc[i] = all[rand() % 4]; //selecte a random base
                
                //track matches to untreated and modified barcodes
                //bound   = A or G at i < 1, T or C at i >= 1
                //unbound = T or C at i < 1, A or G at i >= 1
                switch (bc[i]) {
                    case 'A': (i < 1) ? bnd_cnt++ : unb_cnt++; break;
                    case 'T': (i < 1) ? unb_cnt++ : bnd_cnt++; break;
                    case 'G': (i < 1) ? bnd_cnt++ : unb_cnt++; break;
                    case 'C': (i < 1) ? unb_cnt++ : bnd_cnt++; break;
                    default:
                        printf("mk_rndmzd_bc: error - unexpected character in channel barcode. aborting\n");
                        abort();
                        break;
                }
            }
            bc[i] = '\0';
            
            //check if random barcode meets chnl and mtch specifications
            if (chnl == ERR) { //ERR barcode (<3/4 position match to expected sequence)
                if (bnd_cnt < TDSPLY_MIN_MATCH && unb_cnt < TDSPLY_MIN_MATCH) {
                    bc_made = 1;
                }
            } else if (chnl == BND && mtch == PART) { //partial barcodes require a 3/4 match
                if (bnd_cnt == TDSPLY_MIN_MATCH) {
                    bc_made = 1;
                }
            } else if (chnl == UNB && mtch == PART) { //partial barcodes require a 3/4 match
                if (unb_cnt == TDSPLY_MIN_MATCH) {
                    bc_made = 1;
                }
            } else {
                printf("mk_rndmzd_bc: error - error occurred when generating randomized channel barcode. aborting\n");
                abort();
            }
        }
    }
    
    return;
}

/* print_TDSPLY_fq: construct read sequences and print to fastq file */
void print_TDSPLY_fq(FILE * out_rd1, FILE * out_rd2, char * var_id, char * insrt2use, char * chnl_bc, int end_rnd_typ) {
    
    static int std_cnt = 0; //number of read pairs generated
    
    int i = 0; //general purpose index
    
    char rd1[MAX_LINE] = {0}; //array for generating read 1 sequence
    char rd2[MAX_LINE] = {0}; //array for generating read 2 sequence
    
    const char UMI[13] = "CATCATCATCAT"; //12 nt UMI placeholder, append after channel barcode
    const char prm[51] = "TTATCAAAAAGAGTATTGACTCTTTTACCTCTGGCGGTGATAATGGTTGC"; //PRA1 promoter
    const char ldr[34] = "ATTCGGTGCTCTTCTCTTCGGCCTTCGGGCCAA"; //C3SC1 leader sequence, append 5' end of insert
    
    sprintf(rd2, "%s%s%s%s%s", chnl_bc, UMI, prm, ldr, insrt2use); //construct read 2 sequence
    reverse_complement(rd1, rd2, REVCOMP);                         //revcomp read 2 to obtain read 1 sequence
    //printf("%s\n", rd2);
    
    char end_rnd_code[4][4] = {"SUB", "INS", "DEL", "NAT"}; //code for end randomization stype
    char * end_rnd_ptr = NULL;
    
    //set end randomization type
    switch (end_rnd_typ) {
        case SUB: end_rnd_ptr = &end_rnd_code[SUB][0]; break;
        case INS: end_rnd_ptr = &end_rnd_code[INS][0]; break;
        case DEL: end_rnd_ptr = &end_rnd_code[DEL][0]; break;
        case NAT: end_rnd_ptr = &end_rnd_code[NAT][0]; break;
        default:
            break;
    }
    
    //print fastq line 1 (read id)
    //read id contains:
    //1. prefix 'testdata_TDsply' to indicate read is for use as TECdisplay test data
    //2. id of the target used to generate the read
    //3. end randomization type (NAT, SUB, INS, DEL)
    //4. read id in hexadecimal
    //5. generalized Illumina read suffix
    fprintf(out_rd1, "@testdata_TDsply=%s_%s_0x%08x_R1 1:N:0:INDEX\n", var_id, end_rnd_ptr, std_cnt);
    fprintf(out_rd2, "@testdata_TDsply=%s_%s_0x%08x_R2 2:N:0:INDEX\n", var_id, end_rnd_ptr, std_cnt);
    
    //print fastq line 2 (read sequence)
    fprintf(out_rd1, "%s\n", rd1);
    fprintf(out_rd2, "%s\n", rd2);
    
    //print fastq line 3
    fprintf(out_rd1, "+\n");
    fprintf(out_rd2, "+\n");
    
    //print fastq line 4 (qscore, set as string of 'I' that matches sequence length)
    int len = strlen(rd1);
    for (i = 0; i < (len); i++) {
        fprintf(out_rd1, "I");
        fprintf(out_rd2, "I");
    }
    fprintf(out_rd1, "\n");
    fprintf(out_rd2, "\n");
    
    std_cnt++; //increment read output counter
}

/* mutate_insrt: randomly mutate test data sequencing read insert */
int mutate_insrt(char * vb_map, char * mut)
{
    int mut_pos = rand() % strlen(mut); //select random position for mutation
    int mut_cd = rand() % 3;            //select random type of mutation
    int mut_bs = 0;                     //index for random base
    int cmptbl_sub = 0;                 //flag that mutation positions is compatible with substitution
    int cmptbl_ins = 0;                 //flag that mutation position is compatible with insertion
    
    char * sub_ptr = NULL;              //pointer to substitution base array
    char * ins_ptr = NULL;              //pointer to insertion base array
    char subA[4] = "TGC";               //base array for A substitutions
    char subT[4] = "AGC";               //base array for T substitutions
    char subG[4] = "ATC";               //base array for G substitutions
    char subC[4] = "ATG";               //base array for C substitutions
    char all4[4] = "ATGC";              //base array for insertions
    
    char seg1[MAX_LINE] = {0};          //array for storing first segment of insert
    char tmp2[MAX_LINE] = {0};          //array for storing second segment of insert
    char * seg2 = NULL;                 //pointer to start of second insert segment
        
    if (mut_cd == SUB) {                            //making substitution
        cmptbl_sub = 0;                             //initialize cmptbl_sub to 0
        while (!cmptbl_sub) {                       //find a sub position that...
            if (vb_map[mut_pos] != '*') {           //does not match a variable base location because
                cmptbl_sub = 1;                     //a sub at a variable base may map to another target
            } else {
                mut_pos = rand() % strlen(mut);     //select new positions until a non-variable position is found
            }
        }
        
        switch (mut[mut_pos]) {                     //set pointer to substitution array
            case 'A': sub_ptr = &subA[0]; break;
            case 'T': sub_ptr = &subT[0]; break;
            case 'G': sub_ptr = &subG[0]; break;
            case 'C': sub_ptr = &subC[0]; break;
            default:
                printf("mutate_insrt: error - unrecognized base %c\n", mut[mut_pos]);
                break;
        }
        
        mut_bs = rand() % (strlen(sub_ptr));  //select random base
        mut[mut_pos] = sub_ptr[mut_bs];       //substitute mutant base at mutation position
        
    } else if (mut_cd == INS) {                     //making insertion. the insertion is made after mut_pos
        cmptbl_ins = 0;                             //intitialize cmptbl_ins to 0
        while (!cmptbl_ins) {                       //find an insertion position that...
            if (mut_pos < (strlen(mut)-1)) {        //is not the last position of array because
                cmptbl_ins = 1;                     //trailing insertions are mappable
            } else {
                mut_pos = rand() % strlen(mut);     //select new positions until non-last position is found
            }
        }
        
        strcpy(seg1, mut);         //copy mut to seg1 array
        strcpy(tmp2, mut);         //copy mut to seg2 array
        seg1[mut_pos + 1] = '\0';  //terminate seg1 one nt after the mutation position
        seg2 = &tmp2[mut_pos+1];   //set seg2 pointer to one nt after mutation position
    
        //disallow insertions that duplicate the last base, which would yield a mappable read
        if (mut_pos == (strlen(mut)-2)) { //if insertion is at the second to last position
            switch (mut[mut_pos+1]) {     //select a substitution array that avoids trailing base duplication
                case 'A': ins_ptr = &subA[0]; break;
                case 'T': ins_ptr = &subT[0]; break;
                case 'G': ins_ptr = &subG[0]; break;
                case 'C': ins_ptr = &subC[0]; break;
                default:
                    printf("mutate_insrt: error - unrecognized base %c\n", mut[mut_pos]);
                    break;
            }
        } else { //otherwise, any base is OK
            ins_ptr = &all4[0];
        }
        
        mut_bs = rand() % (strlen(ins_ptr));                  //select random base
        sprintf(mut, "%s%c%s", seg1, ins_ptr[mut_bs], seg2);  //construct insert

    } else if (mut_cd == DEL) {                 //making deletion
        strcpy(seg1, mut);                      //copy mut string to seg1 array
        if (!mut_pos) {                         //if deleting first base...
            sprintf(mut, "%s", &seg1[1]);       //print seg1 to mut array starting at index1
        } else {                                //if deleting non-first base...
            strcpy(tmp2, mut);                  //copy mut string to tmp2 array
            seg1[mut_pos] = '\0';               //terminate seg1 string at mutation position
            seg2 = &tmp2[mut_pos+1];            //start seg2 one nt after mutation position
            sprintf(mut, "%s%s", seg1, seg2);   //construct insert
        }
    }
    
    return mut_cd; //return mutation code
}

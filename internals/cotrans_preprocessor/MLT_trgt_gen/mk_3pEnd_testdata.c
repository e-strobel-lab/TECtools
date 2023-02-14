//
//  mk_3pEnd_testdata.c
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../seq_utils/revcomp.h"

#include "mk_3pEnd_testdata.h"

/* mk_3pEnd_testdata: generate individual fasta files for every intermediate
 length of the input sequence, starting at min_len+1 */
int mk_3pEnd_testdata(char name[MAX_LINE], char seq[MAX_LINE], target3p_genVals * ends, int randomize_end, int mltplir, int min_len)
{
    int i = 0;
    int j = 0;
    int k = 0;
    
    FILE * out_rd1 = NULL;	//read 1 file pointer
    FILE * out_rd2 = NULL;	//read 2 file pointer
    
    //open output files
    char outR1_name[MAX_LINE] = {0}; //read 1 name array
    char outR2_name[MAX_LINE] = {0}; //read 2 name array
    
    //test data without end randomization contains the string NATIVE in the name
    //test data with end randomization contains the string RANDOM in the name
    char nm_NATIVE[7] = "NATIVE";
    char nm_RANDOM[7] = "RANDOM";
    char * nm_typ_ptr = NULL;
    
    if (randomize_end) { //set name type pointer based on randomize_end option
        nm_typ_ptr = &nm_RANDOM[0];
    } else {
        nm_typ_ptr = &nm_NATIVE[0];
    }
    
    //generate read 1 output file
    sprintf(outR1_name, "./%s_targets/test_data/%s_ends_%02dntLong_testdata%s_R1.fq", name, name, ends->len, nm_typ_ptr);
    if ((out_rd1 = fopen(outR1_name, "w")) == NULL) {
        printf("mk_3pEnd_testdata:error: could not generate read one file. Aborting program...\n");
        abort();
    }
    
    //generate read 2 output file
    sprintf(outR2_name, "./%s_targets/test_data/%s_ends_%02dntLong_testdata%s_R2.fq",name, name, ends->len, nm_typ_ptr);
    if ((out_rd2 = fopen(outR2_name, "w")) == NULL) {
        printf("mk_3pEnd_testdata:error: could not generate read two file. Aborting program...\n");
        abort();
    }
    
    /*
     reads are generated for every intermediate transcript length from min_len+1 to the full
     length input sequence. read 2 is assembled from the insert sequence as provided and the
     SC1 leader, UMI, and chnl_bc sequences as they are stored in their respective arrays;
     read 1 is the reverse complement of read 2. this configuration matches what is present
     in real data.
     
     the order of each elemeent in read 2 is:
     [channel barcode (chnl_bc)] - [SC1 leader (ldr)] - [insert (insrt or insrt_rnd)] - [UMI (UMI)]
     
     channel barcode randomization (always performed):
     untreated = YYYRR
     modified  = RRRYY
     
     the exact channel barcode sequence in each read is always randomized, but the resulting reads
     have a predictable composition:
     
     45% of reads contain an untreated barcode
     45% of reads contain a modified barcode
     10% of reads contain an unidentifiable barcode
     90% of reads that contain an untreated or modified barcode match the expected sequence fully (5/5)
     10% of reads that contain an untreated or modified barcode match the expected sequence partially (4/5)
     
     
     3' end seed sequence randomization (optional):
     the region of each read that will be used for 3' end mapping can be optionally modified to contain
     1 nt substitutions, insertions, and deletions. this is useful for assessing how frequently reads
     that do not match the native sequence will map to other transcript lengths and, when this happens,
     how distant the mapped transcript length is from the expected transcript length
     
     when 3' end randomization is performed, only the 3'-most ends.len (sequence length of ends targets)
     will be randomized because these are the nucleotides that will be used for 3' end mapping. constraints
     for each randomization type are described in the rndmz_end function.
     
     */
    
    static const int mk_nmbr[CHNL_CLASSES] = {  405,  405,   45,   45,  100};
    static const int mk_chnl[CHNL_CLASSES] = {  UNT,  MOD,  UNT,  MOD,  ERR};
    static const int mk_mtch[CHNL_CLASSES] = { FULL, FULL, PART, PART, FULL};
    
    char chnl_bc[6] = {0};			//array for storing channel barcode, append to head of read 1
    
    char insrt[MAX_LINE] = {0};		//read 2 insert sequence
    char insrt_rnd[MAX_LINE] = {0};	//read 2 insert sequence with 3' end subs/indels
    char * insrt2use = NULL;		//insert to use when generating read 2 sequence
    
    int ipt_len = strlen(seq);		//length of input sequence
    int end_rnd_typ = NAT;			//code for end type (NAT, SUB, INS, DEL)
    
    
    for (i = min_len+1; i <= ipt_len; i++) { //perform loop for every intermediate transcript
        
        //get sequence for current intermediate transcript
        for (j = 0; j < i; j++) {
            insrt[j] = toupper(seq[j]);
        }
        insrt[j] = '\0';
        
        //generate reads for current intermediate transcript
        for (j = 0; j < CHNL_CLASSES; j++) {
            for (k = 0; k < (mk_nmbr[j] * mltplir); k++) {
                if (randomize_end) { //if randomize end option is on
                    end_rnd_typ = rndmz_end(insrt_rnd, insrt, seq, ends, min_len); //add a sub or indel in insrt_rnd array
                    insrt2use = &insrt_rnd[0];	//use insrt_rnd to make read 2
                } else {						//otherwise
                    insrt2use = &insrt[0];		//use insrt (no randomization) to make read 2
                }
                
                mk_rndmzd_bc(chnl_bc, mk_chnl[j], mk_mtch[j]); //generate channel barcode
                print_fq(out_rd1, out_rd2, insrt2use, chnl_bc, i, end_rnd_typ); //print read to output fastq files
            }
        }
    }
    
    if ((fclose(out_rd1)) == EOF) {
        printf("mk_3pEnd_testdata: error - error occurred when closing test data read 1 file. Aborting program...\n");
        abort();
    }
    
    if ((fclose(out_rd2)) == EOF) {
        printf("mk_3pEnd_testdata: error - error occurred when closing test data read 2 file. Aborting program...\n");
        abort();
    }
    
    char command[MAX_LINE] = {0};
    
    sprintf(command, "gzip ./%s_targets/test_data/%s_ends_%02dntLong_testdata%s_R*.fq",
            name, name, ends->len, nm_typ_ptr);
    system(command);
    
    return 1;
}



/* mk_rndmzd_bc: generate a randomized channel barcode with variable channel and match settings.
 chnl variable:
 UNT = untreated barcode (YYYRR)
 MOD = modified barcode  (RRRYY)
 ERR = undetermined barcode (<4/5 position match to expected sequence)
 
 mtch variable:
 FULL = 5/5 position match to expected sequence
 PART = 4/5 position match to expected sequence
 */
void mk_rndmzd_bc(char * bc, int chnl, int mtch)
{
    int i = 0;
    
    static const char pur[3] = "AG";	//array for randomizing purine (R) bases
    static const char pyr[3] = "TC";	//array for randomizing pyrimidine (Y) bases
    static const char all[5] = "ATGC";	//array for general randomization (N bases)
    
    int unt_cnt = 0;	//counter for matches to untreated barcode
    int mod_cnt = 0;	//counter for matches to modified barcode
    
    int bc_made = 0;	//flag to indicate that random barcode meets chnl and mtch specifications
    
    if (chnl == UNT && mtch == FULL) {	//make untreated full match barcode
        bc[0] = pyr[rand() % 2];		//simple randomization of each base such that
        bc[1] = pyr[rand() % 2];		//the resulting sequence matches YYYRR
        bc[2] = pyr[rand() % 2];
        bc[3] = pur[rand() % 2];
        bc[4] = pur[rand() % 2];
        bc[5] = '\0';
    } else if (chnl == MOD && mtch == FULL) {
        bc[0] = pur[rand() % 2];		//make modified full match barcode
        bc[1] = pur[rand() % 2];		//simple randomization of each base such that
        bc[2] = pur[rand() % 2];		//the resulting sequences matches RRRYY
        bc[3] = pyr[rand() % 2];
        bc[4] = pyr[rand() % 2];
        bc[5] = '\0';
    } else { //procedure for making unmappable and partial match barcodes
        bc_made = 0;
        while (!bc_made) { //repeat until random barcode meets chnl and mtch specifications
            unt_cnt = 0;
            mod_cnt = 0;
            for (i = 0; i < 5; i++) {	//generate a random 5-mer
                bc[i] = all[rand() % 4];
                
                //track matches to untreated and modified barcodes
                //untreated = T or C at i <= 2, A or G at i == 3 and i == 4
                //modified  = A or G at i <= 2, T or C at i == 3 and i == 4
                switch (bc[i]) {
                    case 'A': (i <= 2) ? mod_cnt++ : unt_cnt++; break;
                    case 'T': (i <= 2) ? unt_cnt++ : mod_cnt++; break;
                    case 'G': (i <= 2) ? mod_cnt++ : unt_cnt++; break;
                    case 'C': (i <= 2) ? unt_cnt++ : mod_cnt++; break;
                    default:
                        printf("mk_rndmzd_bc: error - unexpected character in channel barcode. aborting\n");
                        abort();
                        break;
                }
            }
            bc[i] = '\0';
            
            //check if random barcode meets chnl and mtch specifications
            if (chnl == ERR) { //ERR barcode (<4/5 position match to expected sequence)
                if (unt_cnt < 4 && mod_cnt < 4) {
                    bc_made = 1;
                }
            } else if (chnl == UNT && mtch == PART) { //partial barcodes require a 4/5 match
                if (unt_cnt == 4) {
                    bc_made = 1;
                }
            } else if (chnl == MOD && mtch == PART) { //partial barcodes require a 4/5 match
                if (mod_cnt == 4) {
                    bc_made = 1;
                }
            } else {
                printf("mk_rndmzd_bc: error - error occurred when generating randomized channel barcode. aborting\n");
                abort();
            }
        }
    }
}



/* rndmz_end: randomize 3' end region of testdata insert sequence
 NOTE: whereas 3' end target generation aims to avoid generating
 indels that match native sequences, randomized test data generation
 does not.*/
int rndmz_end(char * out, char * ipt, char * seq, target3p_genVals * ends, int min_len) {
    
    int len = strlen(ipt);	//length of input sequence
    int last_indx = len-1;	//last index of input sequence
    
    int end_rnd_typ = (rand() % 3);	//variant type
    int rnd_pos = 0; //randomization position, the index that will be randomized is last_indx-rnd_pos
    
    //rnd_pos index must be < end target length (ends->len)
    //therefore, rnd_pos is set as <random number> % <end target length>
    //
    //randomized position cannot be 0 for deletion variants
    //of the shortest transcript length because this will
    //yield unmappable testdata reads
    if (len == (min_len+1) && end_rnd_typ == DEL) {	  //making deletion in shortest transcript
        while (!(rnd_pos = (rand() % ends->len))) {;} //generate random numbers until rnd_pos is non-zero
    } else {
        rnd_pos = rand() % ends->len; //otherwise, any position is fine
    }

    //variables for handling strings when
    //generating 3' end variants
    char left[MAX_LINE] = {0};
    char tmp[MAX_LINE] = {0};
    char * r_ptr = NULL;
    
    //variable arrays for making random subs and indels
    static const char all[5] = "ATGC";
    static const char notA[4] = "TGC";
    static const char notT[4] = "AGC";
    static const char notG[4] = "ATC";
    static const char notC[4] = "ATG";
    
    if (end_rnd_typ == SUB) { 				//making substitution variant
        strcpy(out, ipt);					//copy input insert to output
        switch (out[last_indx-rnd_pos]) {	//make a substitution in output string
            case 'A': out[last_indx-rnd_pos] = notA[rand() % 3]; break;
            case 'T': out[last_indx-rnd_pos] = notT[rand() % 3]; break;
            case 'G': out[last_indx-rnd_pos] = notG[rand() % 3]; break;
            case 'C': out[last_indx-rnd_pos] = notC[rand() % 3]; break;
            default:
                printf("rndmz_end: error - unexpected non-DNA character. aborting\n");
                abort();
                break;
        }
        
    } else if (end_rnd_typ == INS) {		   //making insertion variant
        if (rnd_pos) {						   //insertion is not at the transcript 3' end
            strcpy(left, ipt);				   //copy input sequence to left string
            left[last_indx-rnd_pos+1] = '\0';  //truncate left string 1 nt to the left of the insertion
            
            strcpy(tmp, ipt);				   //copy input sequence to tmp string
            r_ptr = &tmp[last_indx-rnd_pos+1]; //set pointer to index 1 nt to the right of the insertion
            
            //construct insertion variant. the random nucleotide
            //insertion is made between the left and r_ptr string.
            sprintf(out, "%s%c%s", left, all[rand() % 4], r_ptr);
            
        } else { //adding insertion after 3' end nucleotide
            strcpy(out, ipt);					//copy inputsequence to output array
            out[last_indx+1] = all[rand() % 4]; //append randomized nt to 3' end of output seq
            out[last_indx+2] = '\0';			//append null character to terminate output string
        }
        
    } else if (end_rnd_typ == DEL) {
        strcpy(tmp, ipt); 					  //copy input to tmp array
        
        r_ptr = &tmp[last_indx-rnd_pos+1];	  //set pointer to the right of the deleted base
        tmp[last_indx-rnd_pos] = '\0';		  //set base to delete to zero
        if (rnd_pos) {						  //deleting internal base
            sprintf(out, "%s%s", tmp, r_ptr);
        } else {							  //deleting terminal 3' end base
            strcpy(out, tmp);				  //copy ipt to tmp array
        }
        
    } else {
        printf("rndmz_end: error - unexpected variant types. aborting\n");
        abort();
    }
    
    //print_var_map(ipt, out, rnd_pos, end_rnd_typ);
    
    return end_rnd_typ;
}



/* print_fq: construct read sequences and print to fastq file */
void print_fq(FILE * out_rd1, FILE * out_rd2, char * insrt2use, char * chnl_bc, int end3p, int end_rnd_typ) {
    
    static int cnt = 0; //number of read pairs generated
    
    int i = 0;
    
    char rd1[MAX_LINE] = {0};	//array for generating read 1 sequence
    char rd2[MAX_LINE] = {0};	//array for generating read 2 sequence
    
    static const char ldr[21] = "ATGGCCTTCGGGCCAA"; //SC1 leader sequence, append 5' end of insert
    static const char UMI[10] = "CATCATCAT";		//9 nt UMI placeholder, append to head of read 1
    
    sprintf(rd2, "%s%s%s%s", chnl_bc, ldr, insrt2use, UMI); //construct read 2 sequence
    reverse_complement(rd1, rd2, REVCOMP); 					//revcomp read 2 to obtain read 1 sequence
    
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
    //1. prefix 'testdata' to indicate read is for use as test data
    //2. length of transcript intermediate used to generate the read, format is '3pEND=XXX' where 'XXX' is
    //   the length of the transcript intermediate
    //3. end randomization type (NAT, SUB, INS, DEL)
    //4. read id in hexadecimal
    //5. generalize Illumina read suffix
    fprintf(out_rd1, "@testdata_3pEnd=%03d_%s_0x%07x_R1 1:N:0:INDEX\n", end3p, end_rnd_ptr, cnt);
    fprintf(out_rd2, "@testdata_3pEnd=%03d_%s_0x%07x_R2 2:N:0:INDEX\n", end3p, end_rnd_ptr, cnt);
    
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
    
    cnt++;
}



/* print_var_map: print map that illustrates 3' end randomization
 this is for debugging purposes and is not typically used. it is
 not accessible unless the code is edited directly.
 */
void print_var_map(char * ipt, char * out, int rnd_pos, int mode)
{
    char pipes[MAX_LINE] = {0};	//array for storing pipe character line (used to number nts from 3' end)
    char nums[MAX_LINE] = {0};	//array for storing nt number line
    char star[MAX_LINE] = {0};	//array for storing star charactar line (used to mark variant position)
    
    int len = strlen(ipt);		//length of input sequence
    int last_indx = len-1;		//last index of input sequence
    
    int i = 0;
    int j = 0;
    
    for (i = 0; i < len; i++) { //set pipes, nums, and star arrays to spaces
        pipes[i] = ' ';
        nums[i] = ' ';
        star[i] = ' ';
    }
    pipes[i] = '\0';
    nums[i] = '\0';
    star[i] = '\0';
    
    star[last_indx-rnd_pos] = '*'; //set the nt that was randomized to * in star array
    
    //populate pipes and nums array
    //pipes are placed every 10 nts starting at the 3' end
    //nums are used to mark how many nts inset from the 3' end a pipe is
    char tens[5] = "0123"; //tens characters for use when populating nums array
    for (i = last_indx, j = 0; j <= 3; i-=10, j++) {
        pipes[i] = '|';
        if (i == last_indx) {
            nums[i] = '0';
        } else if (i > 0) {
            nums[i] = tens[j];
            nums[i+1] = '0';
        }
    }
    
    switch (mode) {
        case SUB: printf("\n\nSUB at 3' end index -%02d\n",rnd_pos); break;
        case INS: printf("\n\nINS at 3' end index -%02d\n",rnd_pos); break;
        case DEL: printf("\n\nDEL at 3' end index -%02d\n",rnd_pos); break;
        default:
            break;
    }
    printf("%s\n%s\n%s\n%s\n%s\n",nums, pipes, ipt, star, out);
    
}

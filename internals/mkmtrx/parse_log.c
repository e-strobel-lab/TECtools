//
//  parse_log.c
//  
//
//  Created by Eric Strobel on 10/11/22.
//

#include "parse_log.h"

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "./cotrans_mtrx.h"
#include "./mkmtrx_defs.h"
#include "./mkmtrx_structs.h"
#include "./find_nxt_dir_entry.h"




/* parse_log: open and parse log file to get bowtie2 alignment statistics
 for modified and untreated channels */
int parse_log(char * file_loc, alignment_stats * algn)
{
    //printf("in parse log\n");
    DIR * sub_dirp = NULL;           //pointer to transcript length analysis directory
    
    char log_loc[MAX_LINE+2] = {0};    //location of shapemapper log file
    FILE * log_p = NULL;             //pointer to shapemapper log file
    
    char line[MAX_LINE] = {0};       //array to store lines from log file
    
    char ntry_name[MAX_LINE] = {0};  //full name of next directory entry
    char ntry_prfx[MAX_LINE] = {0};  //prefix of next directory entry, currently used to get transcript length string
    
    
    //open transcript length analysis directory
    if ((sub_dirp = opendir(file_loc)) == NULL) {
        printf("parse_log: error - failed to open analysis directory. aborting...\n");
        abort();
    }
    
    int fnd_log = 0; //flag that log file was found
    
    int modified = 0;         //flag that bowtie2 output for modified channel was found
    int untreated = 0;        //flag that bowtie2 output for untreated channel was found
    int chnl = 0;             //variable to store channel code for use as an index
    int chnls_processed = 0;  //number of channels that have been processd
    int proceed = 0;          //flag to proceed with reading bowtie2 output lines
    
    int reading_paired = 0;   //reading bowtie2 output for paired reads
    int reading_unpaired = 0; //reading bowtie2 output for unpaired reads
    
    char digits[64] = {0};    //array to store read count value strings
    
    
    if ((fnd_log = find_nxt_dir_entry(sub_dirp, "_shapemapper_log.txt", ntry_name, ntry_prfx))) { //find shapemapper log file
        
        //test if log location string will exceed array bounds
        if ((strlen(file_loc) + strlen(ntry_name)) > MAX_LINE) {
            printf("parse_log: error - unexpected long shapemapper2 log path. aborting...");
            abort();
        } else {
            sprintf(log_loc,"%s/%s", file_loc, ntry_name); //construct a location string for the shapemapper log file
        }
        
        
        //open shapemapper log file
        if ((log_p = fopen(log_loc, "r")) == NULL) {
            printf("parse_log: error - failed to open profile file. aborting\n");
            abort();
            
        } else {
            while (fgets(line, MAX_LINE, log_p) != NULL) { //until all lines of the shapemapper log file have been read
                
                if (strstr(line, "BowtieAligner (sample: Modified) output message:")) { //found modified read bowtie output
                    modified = 1;   //set flag that modified read bowtie output was found
                    chnl = MOD;     //set chnl index variable to modified
                } else if (strstr(line, "BowtieAligner (sample: Untreated) output message:")) { //found untreated read bowtie output
                    untreated = 1;  //set flag that untreated read bowtie output was found
                    chnl = UNT;     //set chnl index variable to untreated
                }
                
                if (modified + untreated > chnls_processed) { //if a new bowtie output was found
                    
                    proceed = 1;          //initialize proceed flag to true
                    reading_paired = 0;   //initialize reading_paired flag to false
                    reading_unpaired = 0; //initialize reading_unpaird flag to false
                    
                    //proceed is set to false the bowtie output line that states
                    //the overall alignment rate is found
                    
                    while (proceed) { //while there are still bowtie output lines to parse
                        if (fgets(line, MAX_LINE, log_p) == NULL) {
                            printf("parse_log: error - unanticipated NULL line in bowtie output message. file is truncated. aborting...\n");
                            abort();
                        }
                        get_leading_digits(line, digits); //store read count for current line in digits string

                        /* the code below affirmatively identifies each relevant line of the bowtie
                         output message. lines that contain identical string identifiers are disambiguated
                         by setting the reading_paired and reading_unpaired flags to true or false
                         depending on what section of the output message is being read. the read count
                         for each relevant line is converted to an int or double and stored in the
                         alignment_stats structure */
                        
                        if (strstr(line, "reads; of these:")) { //total read count
                            algn->tot_P[chnl] = atoi(digits);
                            
                        } else if (strstr(line, "were paired; of these:")) { //paired read count
                            algn->prd[chnl] = atoi(digits);
                            reading_paired = 1;   //currently reading output for paired reads, set reading_paired to true..
                            reading_unpaired = 0; //and reading_unpaired to false
                            
                        } else if (strstr(line, "aligned concordantly 0 times") && reading_paired) {
                            algn->prd_c0[chnl] = atoi(digits);
                            
                        } else if (strstr(line, "aligned concordantly exactly 1 time") && reading_paired) {
                            algn->prd_c1[chnl] = atoi(digits);
                            
                        } else if (strstr(line, "aligned concordantly >1 times")  && reading_paired) {
                            algn->prd_cM[chnl] = atoi(digits);
                            
                        } /*else if (strstr(line, "pairs aligned concordantly 0 times; of these:")) {
                            algn->[chnl] = atoi(digits);
                            
                        }*/ else if (strstr(line, "aligned discordantly 1 time")  && reading_paired) {
                            algn->prd_d1[chnl] = atoi(digits);
                            
                        } else if (strstr(line, "pairs aligned 0 times concordantly or discordantly; of these:")  && reading_paired) {
                            algn->prd_0[chnl] = atoi(digits);
                            
                        } else if (strstr(line, "mates make up the pairs; of these:")  && reading_paired) {
                            algn->prd_0m[chnl] = atoi(digits);
                            
                        } else if (strstr(line, "aligned 0 times") && reading_paired) {
                            algn->prd_0m0[chnl] = atoi(digits);
                            
                        } else if (strstr(line, "aligned exactly 1 time") && reading_paired) {
                            algn->prd_0m1[chnl] = atoi(digits);
                            
                        } else if (strstr(line, "aligned >1 times") && reading_paired) {
                            algn->prd_0mM[chnl] = atoi(digits);
                            
                        } else if (strstr(line, "were unpaired; of these:")) { //unpaired read count
                            algn->unp[chnl]  = atoi(digits);
                            reading_paired = 0;   //currently reading output for unpaired reads, set reading_unpaired to true..
                            reading_unpaired = 1; //and reading_paired to false
                            
                        } else if (strstr(line, "aligned 0 times") && reading_unpaired) {
                            algn->unp_0[chnl] = atoi(digits);
                            
                        } else if (strstr(line, "aligned exactly 1 time") && reading_unpaired) {
                            algn->unp_1[chnl] = atoi(digits);
                            
                        } else if (strstr(line, "aligned >1 times") && reading_unpaired) {
                            algn->unp_M[chnl] = atoi(digits);
                            
                        } else if (strstr(line, "overall alignment rate")) {
                            algn->rate[chnl] = atof(digits);
                            proceed = 0;
                        }
                    }
                    chnls_processed++; //increment counter to indicate a channel has been processed
                }
            }
        }
    } else {
        printf("parse_log: error - could not find shapemapper log file in directory %s. aborting...\n", file_loc);
        abort();
    }
    
    return 1;
}

/* get_leading_digits: store the first encountered
 contiguous string of digits in the digits array */
int get_leading_digits(char * line, char * digits) {
    
    int i = 0;  //general purpose index
    int j = 0;  //general purpose index
    
    int fnd_digit = 0;
    
    //iterate to the first digit or negative sign in the line
    for (i = 0, fnd_digit = 0; !fnd_digit && line[i]; i++) {
        if (isdigit(line[i]) || (line[i] == '-' && isdigit(line[i+1]))) {
            fnd_digit = 1;
        }
    }
    i--; //decrement i to index of first digit
    
    if (fnd_digit) {            //found start of digit string
        
        //copy line to digits until the first non-digit-associated character is reached
        for (j = 0; (isdigit(line[i]) || line[i] == '-' || line[i] == '.') && line[i]; i++, j++) {
            digits[j] = line[i];
        }
        digits[j] = '\0';
        
        return 1;
        
    } else { //the line did not contain a string of digits
        digits[0] = '\0';
        return 0;
    }
    
}


/* calc_reads_algned: calculate the total number of reads that were
 assessed and the total number of reads that were aligned */
int calc_reads_algnd(alignment_stats * algn, cotrans_matrix * mtrx)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
        
    int crnt_reads[2] = {0}; //input reads for current transcript length
    int crnt_algnd[2] = {0}; //aligned reads for current transcript length
    
    for (i = mtrx->tl[MIN]; i <= mtrx->tl[MAX]; i++) { //for each transcript length
        
        //initialize total and aligned read counters to 0
        crnt_reads[0] = crnt_reads[1] = 0;
        crnt_algnd[0] = crnt_algnd[1] = 0;
        
        for (j = 0; j < 2; j++) { //for each channel
            
            //calculated total reads
            crnt_reads[j] += (algn[i].prd[j] * 2); //paired reads contain two reads, multiply by two
            crnt_reads[j] += (algn[i].unp[j] * 2); //unpaired reads contain one merged read, multiply by two
                        
            //calculate aligned reads
            crnt_algnd[j] += (algn[i].prd_c1[j] * 2); //paired reads, multiply by two
            crnt_algnd[j] += (algn[i].prd_cM[j] * 2); //paired reads, multiply by two
            crnt_algnd[j] += (algn[i].prd_d1[j] * 2); //paired reads, multiply by two
            crnt_algnd[j] += algn[i].prd_0m1[j];      //singletons
            crnt_algnd[j] += algn[i].prd_0mM[j];      //singletons
            crnt_algnd[j] += (algn[i].unp_1[j] * 2);  //merged reads, multiply by two
            crnt_algnd[j] += (algn[i].unp_M[j] * 2);  //merged reads, multiply by two
            
            algn[i].tot_C[j] = crnt_reads[j]; //set calculated total reads variable in alignment_stats struct
            algn[i].algnd[j] = crnt_algnd[j]; //set aligned reads variable in alignment_stats struct
            algn[i].calc[j] = ((double)(crnt_algnd[j])/(double)(crnt_reads[j]))*100; //calculate per read alignment rate
        }
    }
    
    return 1;
}

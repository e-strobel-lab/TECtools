//
//  bypass_fastp.c
//  
//
//  Created by Eric Strobel on 7/18/25.
//

#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

#include "../../../variant_maker/make_barcodes.h"

#include "../../../utils/io_management.h"

#include "bypass_fastp.h"

/* bypass_fastp: bypass fastp and perform simple read processing during testdata analysis. useful for systems in which fastp is not easily installed. */
void bypass_fastp(char * fq1, char * fq2, FILE ** ifp)
{
    printf("bypassing fastp\n");
    
    FILE *fqin[READ_MAX] = {NULL};  //pointers for input fastq files
    FILE *fqout[READ_MAX] = {NULL}; //pointers for output fastq files
    
    //open fastq input files
    if ((fqin[READ1] = fopen(fq1, "r")) == NULL) {
        printf("bypass_fastp: error - could not open input read one file. Aborting program...\n");
        abort();
    }
    
    if ((fqin[READ2] = fopen(fq2, "r")) == NULL) {
        printf("bypass_fastp: error - could not open input read two file. Aborting program...\n");
        abort();
    }
    
    //open processed fastq output files
    if ((fqout[READ1] = fopen("./split/R1out.fq", "w")) == NULL) {
        printf("bypass_fastp: error - could not open input read one file. Aborting program...\n");
        abort();
    }
    
    if ((fqout[READ2] = fopen("./split/R2out.fq", "w")) == NULL) {
        printf("bypass_fastp: error - could not open input read two file. Aborting program...\n");
        abort();
    }
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
    
    
    char * pRtrm[READ_MAX] = {NULL}; //pointer to UMI-trimmed read sequences
    char * pQtrm[READ_MAX] = {NULL}; //pointer to UMI-trimmed read quality scores
    
    char umi[READ_MAX][MAX_BARCODE_LEN+1] = {{0}}; //array to store UMIs
    
    char * idUP[READ_MAX] = {NULL};  //pointer to upstream segment of read ID
    char * idDN[READ_MAX] = {NULL};  //pointer to downstream segment of read ID
    
    //general variables
    int got_line[READ_MAX] = {1};    //flag to indicate success of the get_line function
    int proceed = 1;                 //flag to indicate that read processing should proceed
    
    //read variables
    char read1[FQ_LINES][MAX_LINE] = {{0}};  //array for storing all four lines of one read 1 fastq entry
    char read2[FQ_LINES][MAX_LINE] = {{0}};  //array for storing all four lines of one read 2 fastq entry
    
    for (i = 0; proceed; i++) {
        for (j = 0; j < FQ_LINES; j++) {    //zero first index of fastq line storage
            read1[j][0] = read2[j][0] = '\0';
        }
        
        //copy fastq file lines for each read to read1 and read2 arrays
        for (j = 0; j < FQ_LINES && proceed; j++) {
            
            got_line[READ1] = get_line(&read1[j][0], fqin[READ1]); //get line for read 1
            got_line[READ2] = get_line(&read2[j][0], fqin[READ2]); //get line for read 2
            
            //test success of get_line for read 1 and read 2
            if (!(got_line[READ1] && got_line[READ2])) {
                //test failed. this should only occur at end of file and get_line will
                //throw an error if an unexpected input line is encountered
                
                //test that the end of read 1 and 2 fastq files were reached at the same time
                if (got_line[READ1] == got_line[READ2]) { //reached EOF for both files
                    proceed = 0; //exit read processing loop
                } else {
                    printf("bypass_fastp: error fastq input files do not contain the same number of lines. aborting...\n");
                    abort();
                }
            }
        }
        
        //store read 1 and read 2 UMIs
        for (j = 0; j < MAX_BARCODE_LEN &&
                    read1[LINE2][j] && read2[LINE2][j] &&
                    read1[LINE4][j] && read2[LINE4][j]; j++) {
            umi[READ1][j] = read1[LINE2][j];
            umi[READ2][j] = read2[LINE2][j];
        }
        umi[READ1][j] = '\0';
        umi[READ2][j] = '\0';
        
        //set pointer to char after last UMI char in both reads and quality scores
        pRtrm[READ1] = &read1[LINE2][j];
        pRtrm[READ2] = &read2[LINE2][j];
        
        pQtrm[READ1] = &read1[LINE4][j];
        pQtrm[READ2] = &read2[LINE4][j];
        
        //find space that delimits upstream and downstream segments of read ids
        for (j = 0; read1[LINE1][j] && read1[LINE1][j] != ' '; j++) {;}
        for (k = 0; read2[LINE1][k] && read2[LINE1][k] != ' '; k++) {;}
        
        //split upstream and downstream segments of read ids
        read1[LINE1][j] = '\0';
        read2[LINE1][k] = '\0';
        
        //set pointer to start of upstream read id
        idUP[READ1] = &read1[LINE1][0];
        idUP[READ2] = &read2[LINE1][0];
        
        //set pointer to start of downstream read id
        idDN[READ1] = &read1[LINE1][j+1];
        idDN[READ2] = &read2[LINE1][k+1];

        //print fastq files
        fprintf(fqout[READ1], "%s:%s_%s %s\n", idUP[READ1], umi[READ1], umi[READ2], idDN[READ1]);
        fprintf(fqout[READ2], "%s:%s_%s %s\n", idUP[READ2], umi[READ1], umi[READ2], idDN[READ2]);
        
        fprintf(fqout[READ1], "%s\n", pRtrm[READ1]);
        fprintf(fqout[READ2], "%s\n", pRtrm[READ2]);
        
        fprintf(fqout[READ1], "+\n");
        fprintf(fqout[READ2], "+\n");
        
        fprintf(fqout[READ1], "%s\n", pQtrm[READ1]);
        fprintf(fqout[READ2], "%s\n", pQtrm[READ2]);
        
    }
    
    //close output fastq files
    if ((fclose(fqout[READ1])) == EOF || (fclose(fqout[READ2])) == EOF) {
            printf("bypass_fastp: error - error occurred when closing split fastq file. Aborting program...\n");
            abort();
    }
    
    //open newly generated fastq files as input files for demultiplexing
    if ((ifp[READ1] = fopen("./split/R1out.fq", "r")) == NULL) {
        printf("bypass_fastp: error - could not open fastp_bypass_R1.fq as read one file. Aborting program...\n");
        abort();
    }
    
    if ((ifp[READ2] = fopen("./split/R2out.fq", "r")) == NULL) {
        printf("bypass_fastp: error - could not open fastp_bypass_R2.fq as read two file.  Aborting program...\n");
        abort();
    }
}

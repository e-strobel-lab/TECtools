//
//  global_defs.h
//  
//
//  Created by Eric Strobel on 6/22/22.
//

#ifndef global_defs_h
#define global_defs_h

#include <stdio.h>

#define MAX_LINE 8192	//maximum line length

/* fastq definitions */
#define READ_MAX 2		//number of reads in a paired end seq run
#define READ1 0 		//specifies read 1
#define READ2 1			//specifies read 2
#define FQ_LINES 4		//number of lines in a fastq file
#define LINE1 0			//specifies fastq line 1
#define LINE2 1			//specifies fastq line 2
#define LINE3 2			//specifies fastq line 3
#define LINE4 3			//specifies fastq line 4

/* target type definitions */
#define FILE_TYPE_INIT 0    //file type initializer
#define VMT_FILE 1          //targets file is a variant maker targets file
#define FASTA_FILE 2        //targets file is a fasta file

#endif /* global_defs_h */

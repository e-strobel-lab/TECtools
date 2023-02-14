//
//  cotrans_io_management.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef io_management_h
#define io_management_h

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

#include "../global/global_defs.h"


/* get_file: open file input and assign to file pointer
 
 ***arguments***
 FILE **ifp: pointer to the file pointer that will be used to open file
 char *ipt: name of file that will be opened
 */
int get_file(FILE **ifp, char *ipt);



/* get_line: get line from file, place into array, remove trailing newline, and return
 line length if successful
 
 ***arguments***
 char *line: pointer to array that will be used to store line
 FILE *ifp: input file pointer
 */
int get_line(char *line, FILE *ifp);



/* get_sample_name: extract sample name from file name
 the sample name is expected to be flanked by /  on
 the left and . on the right. e.g. ./directory/samplename.suffix
 
 ***arguments***
 char * file_name: file name from which sample name will be extracted
 char * sample_name: pointer to array that will contain the sample name
 */
int get_sample_name(char * file_name, char * sample_name);



/* mk_out_dir: generate output directory
 
 ***arguments***
 char * name: name of output directory that will be generated
 */
int mk_out_dir(char * name);

#endif /* io_management_h */

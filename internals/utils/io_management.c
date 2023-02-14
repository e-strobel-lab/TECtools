//
//  io_management.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

#include "../global/global_defs.h"

#include "io_management.h"


/* get_file: open file input and assign to file pointer */
int get_file(FILE **ifp, char *ipt)
{
    if ((*ifp = fopen(ipt, "r")) == NULL) {
        printf("get_file: error - could not open %s. Aborting program...\n", ipt);
        abort();
    }
    return 1;
}


/* get_line: get line from file, place into array, remove trailing newline, and return
 line length if successful */
int get_line(char *line, FILE *ifp)
{
    /* function gets line and replaces terminal newline with null character.
     there is no need for buffering in this case because lines that exceed
     MAX_LINE should not exist and if they do, they are an error.
     the only acceptable mode of failure is to reach the end of the file
     without getting any preceeding characters */
    
    int i = 0;
    char c = 0;
    
    for (i = 0; (c = fgetc(ifp)) != '\n' && c != EOF && c &&  i < MAX_LINE; i++) {
        line[i] = c;
    }
    if (c == '\n' && i != MAX_LINE) {
        line[i] = '\0';			//remove trailing newline
        return i;				//success
    } else if (c == EOF) {		//reached end of file
        if (i == 0) {			//EOF is expected at the start of a line
            return 0;
        } else {				//unexpected EOF
            printf("get_line: error -last line in file not terminated by newline\n");
            abort();
        }
    }else if (!c) {				//unexpected null character
        printf("get_line: error - unanticipated null character\n");
        abort();
    } else if (i == MAX_LINE) {	//unexpected long line
        printf("get_line: error - unanticipated long line\n");
        abort();
    } else {
        printf("get_line: error - reached unreachable code. figure out why. aborting...\n");
        abort();
    }
}


/* get_sample_name: extract sample name from file name
 the sample name is expected to be flanked by /  on
 the left and . on the right. e.g. ./directory/samplename.suffix
 */
int get_sample_name(char * file_name, char * sample_name)
{
    int i = 0;
    int j = 0;
    
    int len = strlen(file_name);	//length of file name
    int nearest_dot = -1;			//index of left-most . character
    int start = -1;					//index of sample name start
    
    //identify index of left-most . character
    for (i = len; file_name[i] != '/' && i >= 0; i--) {
        if (file_name[i] == '.') {
            nearest_dot = i;
        }
    }
    
    //if / character precedes sample name, ignore
    start = (file_name[i] == '/') ? i+1 : 0;
    
    //copy sample name from file name using the bounds established above
    if (nearest_dot == -1) { //failure because no detectable suffix
        printf("get_sample_name: error - input file %s does not contain detectable suffix. aborting...\n", file_name);
        abort();
    } else {
        for (i = start, j = 0; i < nearest_dot && j < (MAX_LINE-1); i++, j++) {
            sample_name[j] = file_name[i];
        }
    }
    sample_name[j] = '\0';
    
    //throw error if sample name is too long for array
    if (j == MAX_LINE-1) {
        printf("get_sample_name: error - input filename %s is too long. aborting...\n", file_name);
        abort();
    }
    
    return 1;
}


/* mk_out_dir: generate output directory*/
int mk_out_dir(char *name)
{
    struct stat st = {0};
    
    if (stat(name, &st) == -1) {
        mkdir(name, 0777);
    } else {
        printf("mk_out_dir: error - failed to make output directory %s. aborting\n", name);
        abort();
    }
    
    return 1;
}

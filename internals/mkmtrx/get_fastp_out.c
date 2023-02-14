//
//  get_fastp_out.c
//  
//
//  Created by Eric Strobel on 10/12/22.
//

#include "get_fastp_out.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"

int get_fastp_out(char * prnt_dir, int * ipt_reads, int * pass_reads)
{
    FILE * fastp_fp = NULL;
    char fastp_fn[MAX_LINE] = {0};
    char line[MAX_LINE] = {0};
    int * val2get = NULL;
    char tmp_val[MAX_LINE] = {0};
    int i = 0;
    int j = 0;
    
    sprintf(fastp_fn, "%s/../fastp.json", prnt_dir);
    
    if ((fastp_fp = fopen(fastp_fn, "r")) == NULL) {
        printf("get_fastp_out: error - failed to find fastp.json file. aborting...\n");
        abort();
    }
    
    int val_cnt = 0;
    while (fgets(line, MAX_LINE, fastp_fp) != NULL && val_cnt < 2) {
        if (strstr(line, "before_filtering") != NULL) {
            val2get = ipt_reads;
        } else if (strstr(line, "after_filtering") != NULL) {
            val2get = pass_reads;
        }
        
        if (val2get != NULL) {
            if (fgets(line, MAX_LINE, fastp_fp) != NULL) {
                if (strstr(line, "total_reads") != NULL) {
                    for (i = 0; line[i] && line[i] != ':'; i++) { ;}
                    for (i++, j = 0; line[i] && isdigit(line[i]); i++, j++) {
                        tmp_val[j] = line[i];
                    }
                    tmp_val[j] = '\0';
                    *val2get = atoi(tmp_val);
                    
                } else {
                    printf("get_fastp_out: error - unanticipated format for fastp.json file. aborting...\n");
                    abort();
                }
            } else {
                printf("get_fastp_out: error - unanticipated file end. aborting\n");
                abort();
            }
            
            tmp_val[0] = '\0';
            val2get = NULL;
            val_cnt++;
        }
    }
    
    if (fclose(fastp_fp) == EOF) { //close shapemapper output file
        fprintf(fastp_fp, "get_fastp_out: error - failed to close fastp.json file. Aborting program...\n");
        abort();
    }
    
    printf("input reads:\t%d\n", *ipt_reads);
    printf("passed reads:\t%d\n", *pass_reads);
    
    return 1;
}

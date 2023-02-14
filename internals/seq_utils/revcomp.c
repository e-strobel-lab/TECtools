//
//  revcomp.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "revcomp.h"

/* reverse_complement: reverse complement or reverse input DNA sequence */
int reverse_complement(char *out, char *ipt, int mode) {
    int i = 0;
    int j = 0;
    
    int len = 0;		//length of sequence
    char rcBase = 0;	//stores complement of input sequence base
    int vbase = 0;		//flag indicating that sequence contains an n base
    
    len = strlen(ipt);	//determine length of input string
    
    for (i = 0, j = len-1; j >= 0; i++, j--) {
        if (mode == REVCOMP) {	//reverse complement mode
            switch (ipt[i]) {
                case 'A': rcBase = 'T'; break;
                case 'T': rcBase = 'A'; break;
                case 'G': rcBase = 'C'; break;
                case 'C': rcBase = 'G'; break;
                case 'N':
                    rcBase = 'N';
                    vbase = 1;
                    break;
                case 'a': rcBase = 't'; break;
                case 't': rcBase = 'a'; break;
                case 'g': rcBase = 'c'; break;
                case 'c': rcBase = 'g'; break;
                case 'n':
                    rcBase = 'n';
                    vbase = 1;
                    break;
                default:
                    //error
                    printf("reverse_complement: error - unexpected base. aborting...");
                    abort();
                    break;
            }
            out[j] = rcBase;
        } else if (mode == REVERSE) {	//reverse mode
            out[j] = ipt[i];
        } else {
            printf("error: revcomp - unexpected mode. aborting");
            abort();
        }
    }
    out[len] = '\0';
    
    if (ipt[i]) { //did not reach null character for input string i when j == 0
        printf("reverse_complement: error - output array is not large enough to reverse complement the entire input string. aborting...");
        abort();
    }
    
    return vbase;
}

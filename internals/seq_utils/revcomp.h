//
//  revcomp.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef revcomp_h
#define revcomp_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* revcomp mode definitions */
#define REVCOMP 2	//specifies reverse complement mode for rc function
#define REVERSE 1	//specifies reverse mode for rc function

/* reverse_complement: reverse complement or complement a DNA sequence
 
 ***arguments***
 char *out: pointer to output sequence char array
 char *ipt:	pointer to input sequence char array
 int mode:  specifies REVCOMP or REVERSE mode
 */

int reverse_complement(char *out, char *ipt, int mode);

#endif /* revcomp_h */

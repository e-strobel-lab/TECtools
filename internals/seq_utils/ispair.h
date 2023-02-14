//
//  ispair.h
//  
//
//  Created by Eric Strobel on 6/21/22.
//

#ifndef ispair_h
#define ispair_h

#include <stdio.h>

#define NO_PAIR 0
#define GU_PAIR 1
#define AT_PAIR 2
#define GC_PAIR 3

/* ispair: test whether two nucleotides constitute a G:C, A:U/T, or G:U/T pair.
   returns a distinct value for each type of pair and 0 for no pair. */
int ispair(char nt1, char nt2);

#endif /* ispair_h */

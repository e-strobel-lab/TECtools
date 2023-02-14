//
//  ispair.c
//  
//
//  Created by Eric Strobel on 6/21/22.
//

#include <stdio.h>

#include "ispair.h"

/* ispair: test whether two nucleotides constitute a G:C, A:U/T, or G:U/T pair.
   returns a distinct value for each type of pair and 0 for no pair. */
int ispair(char nt1, char nt2)
{
    if ((nt1 == 'G' && nt2 == 'C') ||  (nt1 == 'C' && nt2 == 'G')){
        return 3;
    } else if ((nt1 == 'A' && nt2 == 'T') ||  (nt1 == 'T' && nt2 == 'A')){
        return 2;
    } else if ((nt1 == 'G' && nt2 == 'T') ||  (nt1 == 'T' && nt2 == 'G')){
        return 1;
    } else if ((nt1 == 'A' && nt2 == 'U') ||  (nt1 == 'U' && nt2 == 'A')){
        return 2;
    } else if ((nt1 == 'G' && nt2 == 'U') ||  (nt1 == 'U' && nt2 == 'G')){
        return 1;
    } else {
        return 0;
    }
}

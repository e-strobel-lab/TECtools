//
//  t2u.c
//  
//
//  Created by Eric Strobel on 6/21/22.
//

#include <stdio.h>

#include "t2u.h"

/* t2u: converts t bases to u bases */
int t2u(char * sq)
{
    int i = 0;
    
    for (i = 0; sq[i]; i++) {
        switch (sq[i]) {
            case 'T':
                sq[i] = 'U';
                break;
            case 't':
                sq[i] = 'u';
                break;
            default:
                break;
        }
    }
    
    return 1;
}

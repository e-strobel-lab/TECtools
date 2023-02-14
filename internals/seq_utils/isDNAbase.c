//
//  isDNAbase.c
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#include <stdio.h>

#include "isDNAbase.h"

/* isbase: test whether character c is a specified DNA base (ATGC/atgc) */
int isDNAbase(char c)
{
    switch (c) {
        case 'A': return 1; break;
        case 'T': return 1; break;
        case 'G': return 1; break;
        case 'C': return 1; break;
        case 'a': return 1; break;
        case 't': return 1; break;
        case 'g': return 1; break;
        case 'c': return 1; break;
        default: break;
    }
    return 0;
}

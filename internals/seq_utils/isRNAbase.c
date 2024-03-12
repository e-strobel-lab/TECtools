//
//  isRNAbase.c
//  
//
//  Created by Eric Strobel on 2/13/24.
//

#include <stdio.h>

#include "isRNAbase.h"

int isRNAbase(char c)
{
    switch (c) {
        case 'A': return 1; break;
        case 'U': return 1; break;
        case 'G': return 1; break;
        case 'C': return 1; break;
        case 'a': return 1; break;
        case 'u': return 1; break;
        case 'g': return 1; break;
        case 'c': return 1; break;
        default: break;
    }
    return 0;
}

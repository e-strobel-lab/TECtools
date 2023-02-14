//
//  isIUPACbase.c
//  
//
//  Created by Eric Strobel on 6/21/22.
//

#include <stdio.h>

#include "isIUPACbase.h"

/* isIUPACbase: test that character matches IUPAC notation */
int isIUPACbase(char c) {
    switch (c) {
        case 'A': return 1; break;
        case 'T': return 1; break;
        case 'G': return 1; break;
        case 'C': return 1; break;
        case 'U': return 1; break;
        case 'N': return 1; break;
        case 'R': return 1; break;
        case 'Y': return 1; break;
        case 'K': return 1; break;
        case 'M': return 1; break;
        case 'S': return 1; break;
        case 'W': return 1; break;
        case 'B': return 1; break;
        case 'D': return 1; break;
        case 'H': return 1; break;
        case 'V': return 1; break;
        case 'a': return 1; break;
        case 't': return 1; break;
        case 'g': return 1; break;
        case 'c': return 1; break;
        case 'u': return 1; break;
        case 'n': return 1; break;
        case 'r': return 1; break;
        case 'y': return 1; break;
        case 'k': return 1; break;
        case 'm': return 1; break;
        case 's': return 1; break;
        case 'w': return 1; break;
        case 'b': return 1; break;
        case 'd': return 1; break;
        case 'h': return 1; break;
        case 'v': return 1; break;
        default: break;
    }
    return 0;
}

//
//  is_dgnrt_mtch.c
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#include <stdio.h>
#include <stdlib.h>

#include "./isIUPACbase.h"
#include "is_dgnrt_mtch.h"

/* is_dgnrt_mtch: test if a base is a match to a degenerate base */
int is_dgnrt_mtch(char bs, char dgnrt_bs)
{
    int match = 0;
    
    if (!isIUPACbase(bs)) {
        printf("is_dgnrt_mtch: error - unrecognized base %c. aborting...\n", bs);
        abort();
    }
    
    if (bs == dgnrt_bs) {    //self-match
        match = 1;
    } else {
        switch (dgnrt_bs) {
            case 'R': if (bs == 'A' || bs == 'G') {match = 1;} break;
            case 'Y': if (bs == 'T' || bs == 'C') {match = 1;} break;
            case 'M': if (bs == 'A' || bs == 'C') {match = 1;} break;
            case 'K': if (bs == 'G' || bs == 'T') {match = 1;} break;
            case 'S': if (bs == 'G' || bs == 'C') {match = 1;} break;
            case 'W': if (bs == 'A' || bs == 'T') {match = 1;} break;
                
            case 'B': if (bs == 'T' || bs == 'G' || bs == 'C' ||
                          bs == 'Y' || bs == 'K' || bs == 'S') {match = 1;} break;
                
            case 'D': if (bs == 'A' || bs == 'T' || bs == 'G' ||
                          bs == 'R' || bs == 'K' || bs == 'W') {match = 1;} break;
                
            case 'H': if (bs == 'A' || bs == 'T' || bs == 'C' ||
                          bs == 'Y' || bs == 'M' || bs == 'W') {match = 1;} break;
                
            case 'V': if (bs == 'A' || bs == 'G' || bs == 'C' ||
                          bs == 'R' || bs == 'M' || bs == 'S') {match = 1;} break;
                
            case 'N': match = 1; break;
                
            default:
                printf("is_dgnrt_mtch: error - unrecognized degenerate base %c. aborting...\n", dgnrt_bs);
                abort();
                break;
        }
    }
    return match;
}

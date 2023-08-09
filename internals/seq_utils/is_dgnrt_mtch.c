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

//TODO: consider having different success codes for matches against a degenerate base vs a non-degenerate base

/* is_dgnrt_mtch: test if a base is a match to a degenerate base */
int is_dgnrt_mtch(char bs1, char bs2)
{
    int match = 0;
    
    if (!isIUPACbase(bs1)) {
        printf("is_dgnrt_mtch: error - unrecognized base %c. aborting...\n", bs1);
        abort();
    }
    
    if (!isIUPACbase(bs2)) {
        printf("is_dgnrt_mtch: error - unrecognized base %c. aborting...\n", bs2);
        abort();
    }
    
    if (bs1 == bs2) { //self-match, this catches non degnerate bases too
        match = 1;
    } else {
        switch (bs2) {
            case 'R': if (bs1 == 'A' || bs1 == 'G') {match = 1;} break;
            case 'Y': if (bs1 == 'T' || bs1 == 'C') {match = 1;} break;
            case 'M': if (bs1 == 'A' || bs1 == 'C') {match = 1;} break;
            case 'K': if (bs1 == 'G' || bs1 == 'T') {match = 1;} break;
            case 'S': if (bs1 == 'G' || bs1 == 'C') {match = 1;} break;
            case 'W': if (bs1 == 'A' || bs1 == 'T') {match = 1;} break;
                
            case 'B': if (bs1 == 'T' || bs1 == 'G' || bs1 == 'C' ||
                          bs1 == 'Y' || bs1 == 'K' || bs1 == 'S') {match = 1;} break;
                
            case 'D': if (bs1 == 'A' || bs1 == 'T' || bs1 == 'G' ||
                          bs1 == 'R' || bs1 == 'K' || bs1 == 'W') {match = 1;} break;
                
            case 'H': if (bs1 == 'A' || bs1 == 'T' || bs1 == 'C' ||
                          bs1 == 'Y' || bs1 == 'M' || bs1 == 'W') {match = 1;} break;
                
            case 'V': if (bs1 == 'A' || bs1 == 'G' || bs1 == 'C' ||
                          bs1 == 'R' || bs1 == 'M' || bs1 == 'S') {match = 1;} break;
                
            case 'N': match = 1; break;
                
            case 'A': if (bs1 == 'N' || bs1 == 'R' || bs1 == 'M' || bs1 == 'W' ||
                          bs1 == 'D' || bs1 == 'H' || bs1 == 'V') {match = 1;} break;
            
            case 'T': if (bs1 == 'N' || bs1 == 'Y' || bs1 == 'K' || bs1 == 'W' ||
                          bs1 == 'B' || bs1 == 'D' || bs1 == 'H') {match = 1;} break;
                
            case 'G': if (bs1 == 'N' || bs1 == 'R' || bs1 == 'K' || bs1 == 'S' ||
                          bs1 == 'B' || bs1 == 'D' || bs1 == 'V') {match = 1;} break;
                
            case 'C': if (bs1 == 'N' || bs1 == 'Y' || bs1 == 'M' || bs1 == 'S' ||
                          bs1 == 'B' || bs1 == 'H' || bs1 == 'V') {match = 1;} break;
                
            default:
                printf("is_dgnrt_mtch: error - unrecognized degenerate base %c. aborting...\n", bs2);
                abort();
                break;
        }
    }
    
    return match;
}

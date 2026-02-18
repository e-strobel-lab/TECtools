//
//  mk_fasta.h
//  
//
//  Created by Eric Strobel on 1/17/23.
//

#ifndef mk_fasta_h
#define mk_fasta_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "./isIUPACbase.h"

int mk_fasta_file(char * nm, char * sq, char * path);

#endif /* mk_fasta_h */

//
//  variant_maker_defs.h
//  
//
//  Created by Eric Strobel on 8/3/22.
//

#ifndef variant_maker_defs_h
#define variant_maker_defs_h

#include <stdio.h>

#include "../global/global_defs.h"

//mode values
#define MAKE_VARIANTS 1     //make variants mode
#define MAKE_BARCODES 2     //make barcodes mode

#define MAXREF 8            //maximum number of references sequences
#define MAX_COMPLEX 26      //maximum number of complex pairs
#define MAX_BRCDS_2_MK 5000 //maximum number of barcodes to make

#define CALCULATED 0        //index for calculated variant count in basemap struct
#define EXPANDED   1        //index for expanded variant count in basemap struct
#define FILTERED   2        //index for filtered variant count in basemap struct  

//basemap pairs array values
#define CMPLX_PAIR_ROOT -1    //indicates base is the first mate in a complex pair
#define NO_PAIR_CONSTRAINT -2 //indicates no pairing constraint

#endif /* variant_maker_defs_h */

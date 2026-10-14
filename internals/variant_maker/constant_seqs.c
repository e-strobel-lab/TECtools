//
//  constant_seqs.c
//  
//
//  Created by Eric Strobel on 6/25/25.
//

#include "constant_seqs.h"

//priming sites

char pra1[51]            = "ttatcaaaaagagtattgactcttttacctctggcggtgataatggttgc"; //PRA1 promoter
char pra1_sc1[41]        = "acctctggcggtgataatggttgcatggccttcgggccaa";
char pra1_m25_to_m50[27] = "ttatcaaaaagagtattgactctttt";
char c3sc1[34]           = "attcggtgctcttctcttcggccttcgggccaa";
char vra3[27]            = "gatcgtcggactgtagaactctgaac";
char vra5[22]            = "ccttggcacccgagaattcca";
char RLA29synch_3p11[34] = "aaacaccaccaAcatcaccatcatcctgacgcg"; //v2
char RLA29synch_inFP[12] = "aaacaccacca"; //sequence in RNAP footprint when RNAP is at synch site

char sc1[17] = "atggccttcgggccaa"; //SC1 sequence

//array of barcoded C3SC sequences
char c3sc_BC[C3SC_BRCD_CNT][34] = {"ATTCGGTGCTCTTCTCTTCGGGGTTCGCCCCAA",
                                   "ATTCGGTGCTCTTCTCTTCGGGCTTCGGCCCAA",
                                   "ATTCGGTGCTCTTCTCTTCGGCGTTCGCGCCAA",
                                   "ATTCGGTGCTCTTCTCTTCGGCCTTCGGGCCAA",
                                   "ATTCGGTGCTCTTCTCTTCGCGGTTCGCCGCAA",
                                   "ATTCGGTGCTCTTCTCTTCGCGCTTCGGCGCAA",
                                   "ATTCGGTGCTCTTCTCTTCGCCGTTCGCGGCAA",
                                   "ATTCGGTGCTCTTCTCTTCGCCCTTCGGGGCAA",
                                   "ATTCGGTGCTCTTCTCTTCCGGGTTCGCCCGAA",
                                   "ATTCGGTGCTCTTCTCTTCCGGCTTCGGCCGAA",
                                   "ATTCGGTGCTCTTCTCTTCCGCGTTCGCGCGAA",
                                   "ATTCGGTGCTCTTCTCTTCCGCCTTCGGGCGAA",
                                   "ATTCGGTGCTCTTCTCTTCCCGGTTCGCCGGAA",
                                   "ATTCGGTGCTCTTCTCTTCCCGCTTCGGCGGAA",
                                   "ATTCGGTGCTCTTCTCTTCCCCGTTCGCGGGAA",
                                   "ATTCGGTGCTCTTCTCTTCCCCCTTCGGGGGAA"};

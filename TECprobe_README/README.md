# TECtools - TECprobe tools

This sub-directories within this directory contain README files for the following programs:

`cotrans_preprocessor` - Performs sequencing read preprocessing to prepare TECprobe-VL and TECprobe-LM data for analysis by ShapeMapper2 (https://github.com/Weeks-UNC/shapemapper2)

`process_TECprobeVL_profiles` - performs whole-dataset normalization for TECprobe-VL and TECprobe-LM data, merges replicate data, and assembles ShapeMapper2 output files into a single csv file that is compatible with TECprobe visualization tools.

`mkmtrx` - Assembles TECprobe-VL and TECprobe-LM data into matrix and other useful formats, reports alignment rates, and can be used to generate rdat files

`mtrx2cols` - Extracts reactivity trajectories for specific nucleotides from a reactivity matrix

`assemble_TECprobeLM_data` - assembles target transcript length reactivity profiles or transcript length distributions for each sample of a TECprobe-LM experiment into a single file.

`draw_intermediates` - performs reactivity-constrained minimum free energy RNA structure prediction for intermediate transcripts using the RNAstructure (Reuter, J. S. & Mathews, D. H., BMC Bioinformatics) Fold algorithm.

# TECtools

TECtools is a suite of tools for processing data from Transcription Elongation Complex Display (TECdisplay) and Transcription Elongation Complex RNA structure probing (TECprobe) experiments. README files for TECdisplay and TECprobe tools are available in the TECdisplay_README and TECprobe_README directories, respectively.

### TECdisplay tools:

`variant_maker` - Generates a targets file for user-specified sequence variants.

`TECdisplay_mapper` - Processes and maps sequencing reads to targets.

`merge_TECdisplay_replicates` - Merges data from replicate TECdisplay experiments.

`calc_FracBound_difference` - Calculates the difference in fraction bound for all variants in two TECdisplay datasets.

`TECdisplay_navigator` - Filters data for variants that match user-specified constraints.

`TECdisplay_Hnav` - Hierarchically filters data for variants that match a series of user-specified constraints.

`id2variant` - Reconstructs complete variant sequences from variant ids.


### TECprobe tools:

`cotrans_preprocessor` - Performs sequencing read preprocessing to prepare TECprobe-VL and TECprobe-LM data for analysis by ShapeMapper2 (https://github.com/Weeks-UNC/shapemapper2)

`process_TECprobeVL_profiles` - performs whole-dataset normalization for TECprobe-VL and TECprobe-LM data, merges replicate data, and assembles ShapeMapper2 output files into a single csv file that is compatible with TECprobe visualization tools.

`mkmtrx` - Assembles TECprobe-VL and TECprobe-LM data into matrix and other useful formats, reports alignment rates, and can be used to generate rdat files

`mtrx2cols` - Extracts reactivity trajectories for specific nucleotides from a reactivity matrix

`assemble_TECprobeLM_data` - assembles target transcript length reactivity profiles or transcript length distributions for each sample of a TECprobe-LM experiment into a single file.

`draw_intermediates` - performs reactivity-constrained minimum free energy RNA structure prediction for intermediate transcripts using the RNAstructure (Reuter, J. S. & Mathews, D. H., BMC Bioinformatics) Fold algorithm.


## Compiling TECtools

**TECtools will only run on Linux and MacOS systems. Windows is not currently supported**

### TECdisplay:

TECtools executables for TECdisplay data processing and analysis can be compiled using the following commands:

```
cd ./TECtools-1.1.0/build_TECdisplay

sh ./build_TECdisplay.sh
```

This will generate the executables `variant_maker`, `TECdisplay_mapper`, `merge_TECdisplay_replicates`, `calc_FracBound_difference`, `TECdisplay_navigator`, `TECdisplay_Hnav`, and `id2variant` in the build_TECdisplay directory.



### TECprobe:

TECtools executables for TECprobe data processing and analysis can be compiled using the following commands:

```
cd ./TECtools-1.1.0/build_TECprobe

sh ./build_TECprobe.sh
```

This will generate the executables `cotrans_preprocessor`, `process_TECprobeVL_profiles`, `mkmtrx`, `mtrx2cols`, `assemble_TECprobeLM_data`, and `draw_intermediates` in the build_TECprobe directory.

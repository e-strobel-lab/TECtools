# TECtools

TECtools is a suite of tools for processing data from Transcription Elongation Complex Display (TECdisplay) and Transcription Elongation Complex RNA structure probing (TECprobe) experiments. README files for TECdisplay and TECprobe tools are available in the TECdisplay_README and TECprobe_README directories, respectively.

### TECdisplay tools:

**variant_maker** - Generates a targets file for user-specified sequence variants.

**TECdisplay_mapper** - Processes and maps sequencing reads to targets.

**merge_TECdisplay_replicates** - Merges data from replicate experiments.

**TECdisplay_navigator** - Filters data for variants that match user-specified constraints.

**TECdisplay_Hnav** - Hierarchically filters data for variants that match a series of user-specified constraints.

**id2variant** - Reconstructs complete variant sequence from variant id.



### TECprobe tools:

**cotrans_preprocessor** - Performs sequencing read preprocessing to prepare data for analysis by ShapeMapper2 (https://github.com/Weeks-UNC/shapemapper2)

**mkmtrx** - Assembles data into matrix and other useful formats, reports alignment rates, and can be used to generate rdat files

**mtrx2cols** - Extracts reactivity trajectories for specific nucleotides from a reactivity matrix



## Compiling TECtools scripts

**TECtools will only run on Linux and MacOS systems. Windows is not currently supported**

### TECdisplay:

The `variant_maker`, `TECdisplay_mapper`, `merge_TECdisplay_replicates`, `TECdisplay_navigator`, `TECdisplay_Hnav`, and `id2variant` executables can be compiled using the following commands:

```
cd ./TECtools-1.0.0/build_TECdisplay

sh ./build_TECdisplay.sh
```

This will generate the executables `variant_maker`, `TECdisplay_mapper`, `merge_TECdisplay_replicates`, `TECdisplay_navigator`, `TECdisplay_Hnav`, and `id2variant` in the build_TECdisplay directory.



### TECprobe:

The `cotrans_preprocessor`, `mkmtrx`, and `mtrx2cols` executables can be compiled using the following commands:

```
cd ./TECtools-1.0.0/build_TECprobe

sh ./build_TECprobe.sh
```

This will generate the executables `cotrans_preprocessor`, `mkmtrx`, and `mtrx2cols` in the build_TECprobe directory.

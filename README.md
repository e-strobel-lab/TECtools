# TECtools_public

TECtools is a suite of tools for processing data from Transcription Elongation Complex RNA structure probing (TECprobe) data including:

**cotrans_preprocessor** - Performs sequencing read preprocessing to prepare data for analysis by ShapeMapper 2

**mkmtrx** - Assembles data into matrix and other useful formats, reports alignment rates, and can be used to generate rdat files

**mtrx2cols** - Extracts reactivity trajectories for specific nucleotides from a reactivity matrix

**TECtools is only intended for use on Linux and MacOS systems**

## Compiling TECtools scripts

The cotrans_preprocessor, mkmtrx, and mtrx2cols scripts can be compiled using the following commands:

```
cd ./TECtools_public/build_TECprobe

sh ./build_TECprobe.sh
```

This will generate the executables cotrans_preprocessor, mkmtrx, and mtrx2cols in the build_TECprobe directory.

## Processing TECprobe data using cotrans_preprocessor

cotrans_preprocessor performs sequencing read preprocessing to prepare data for analysis by ShapeMapper 2

### Set cotrans_preprocessor run mode (required)

```
-m/--mode <run_mode_specifier>

run_mode_specifiers:
  MAKE_FASTA - Generate FASTA file to be used for target generation
  MAKE_3pEND_TARGETS - Generate 3' end targets to be used for demultiplexing FASTQ files by transript length, and intermediate transcript targets to be used for sequencing read alignment.
  PROCESS_MULTI - Perform preprocessing for TECprobe-ML experiments
  PROCESS_SINGLE - Perform preprocessing for TECprobe-SL experiments
  MAKE_RUN_SCRIPT - Generate shell script for running ShapeMapper2
```

### MAKE_FASTA mode inputs and options

Required:

```
-n/--seq-name <name_of_sequence>
-s/--sequence <sequence>
```

### MAKE_3pEnd_TARGETS mode inputs and options

Required:
```
-A/--fasta  <fasta_file_containing_target_sequence>
```

Optional:
```
-E/--end-length 
-S/--min-length
-T/--test-data
-U/--multiplier
-R/--end-randomization
```

### PROCESS_MULTI and PROCESS_SINGLE mode inputs and options

Required for PROCESS_MULTI and PROCESS_SINGLE modes:

```
-i/--read1
-I/--read2
```

Required for PROCESS_MULTI mode only:

```
-e/--3pEnd
```

Required for PROCESS_SINGLE mode only:

```
-a/--fasta-ref
```

Optional:

```
-p/--fastp-path
-l/--limit
-t/--testdata
```

### MAKE_RUN_SCRIPT mode input

Required:
```
-c/--config
```

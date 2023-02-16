# TECtools

TECtools is a suite of tools for processing data from Transcription Elongation Complex RNA structure probing (TECprobe) data including:

**cotrans_preprocessor** - Performs sequencing read preprocessing to prepare data for analysis by ShapeMapper2 (https://github.com/Weeks-UNC/shapemapper2)

**mkmtrx** - Assembles data into matrix and other useful formats, reports alignment rates, and can be used to generate rdat files

**mtrx2cols** - Extracts reactivity trajectories for specific nucleotides from a reactivity matrix



## Compiling TECtools scripts

**TECtools will only run on Linux and MacOS systems. Windows is not currently supported**

The `cotrans_preprocessor`, `mkmtrx`, and `mtrx2cols` scripts can be compiled using the following commands:

```
cd ./TECtools_public/build_TECprobe

sh ./build_TECprobe.sh
```

This will generate the executables `cotrans_preprocessor`, `mkmtrx`, and `mtrx2cols` in the build_TECprobe directory.



## Processing TECprobe data using cotrans_preprocessor

`cotrans_preprocessor` performs sequencing read preprocessing to prepare data for analysis by `ShapeMapper2`

**`cotrans_preprocessor` requires that `fastp` (https://github.com/OpenGene/fastp) is installed**

**analysis of TECprobe using `cotrans_preprocessor` requires that `ShapeMapper2` (https://github.com/Weeks-UNC/shapemapper2) is installed** 

### Set cotrans_preprocessor run mode (required for all run modes)

```
-m/--mode <run_mode_specifier>

run_mode_specifiers:

  MAKE_FASTA          Generate FASTA file to be used for target generation
  
  MAKE_3pEND_TARGETS  Generate 3' end targets to be used for demultiplexing FASTQ files by 
                      transript length, and intermediate transcript targets to be used for 
                      sequencing read alignment. 3' end targets contain all native and 1 nt
                      substitution, insertion, and deletion variants of the last 14 nt 
                      (or more, if specified) of every possible intermediate transcript.
  
  PROCESS_MULTI       Perform preprocessing for TECprobe-ML experiments
  
  PROCESS_SINGLE      Perform preprocessing for TECprobe-SL experiments
  
  MAKE_RUN_SCRIPT     Generate shell script for running ShapeMapper2
```


### MAKE_FASTA mode inputs and options

Required:

```
-n/--seq-name <name_of_sequence>  Target name.

-s/--sequence <sequence>          Target sequence.
```


### MAKE_3pEnd_TARGETS mode inputs and options

Required:
```
-A/--fasta  <fasta_file>  Target sequence in FASTA format. Can be generated by running
                          cotrans_preprocessor in MAKE_FASTA mode.
```

Optional:
```
-E/--end-length <n>     Length of 3' end targets. Default=14, 3' end target length should ideally
                        be 14-16 nt, because the E. coli RNA polymerase footprint on RNA is 14-15 nt.
                        
-S/--min-length <n>     Minimum intermediate transcript target length. Default=20.

-T/--test-data          Generate test data to be used for validating transcript 3' end mapping.
                        By default, generates 1000 test_data reads per intermediate transcript
                        using native sequence with randomized channel barcodes.
                        
-U/--multiplier <n>     Multiply the number of test_data reads generated by <n>.

-R/--end-randomization  Make random 1 nt substitutions, insertions, and deletions in the last
                        <end_length> nts of test_data reads. End-randomization should be used
                        to validate correct 3' end mapping. Complete coverage of all 1 nt
                        variations typically requires that --multiplier be set to at least 30.
```


### PROCESS_MULTI and PROCESS_SINGLE mode inputs and options

Required for PROCESS_MULTI and PROCESS_SINGLE modes:

```
-i/--read1 <read1_fastq>      Read 1 FASTQ file input.

-I/--read2 <read2_fastq>      Read 2 FASTQ file input.
```

Required for PROCESS_MULTI mode only:

```
-e/--3pEnd <3'_end_targets>   3' end targets file input.
```

Required for PROCESS_SINGLE mode only:

```
-a/--fasta-ref <fasta_file>   FASTA file containing single target.
```

Optional:

```
-p/--fastp-path <fastp_path>  Path to fastp executable

-t/--testdata                 Enable test data analysis. This option should only be used
                              when analyzing test_data sequencing reads that were 
                              generated in MAKE_3pEND_TARGETS mode using the -T option.
```


### MAKE_RUN_SCRIPT mode input

Required:
```
-c/--config <config_file>     Run script configuration file input.
```

### Basic usage of cotrans_preprocessor

1.  Generate a target sequence FASTA file by running `cotrans_preprocessor` in `MAKE_FASTA` mode:
    
    `cotrans_preprocessor -m MAKE_FASTA -n <target_RNA_name> -s <target_RNA_sequence>`
    
    This will generate a FASTA file named `<target_RNA_name>.fa` with the sequence `<target_RNA_sequence>`.
    

2.  If analyzing TECprobe-ML data, generate 3' end and intermediate transcript targets by
    running `cotrans_preprocessor` in `MAKE_3pEND_TARGET` mode:
    
    without test_data:
    
    `cotrans_preprocessor -m MAKE_3pEND_TARGETS -A <target_RNA_fasta_file>`
    
    with native 3' end test_data generation:
    
    `cotrans_preprocessor -m MAKE_3pEND_TARGETS -A <target_RNA_fasta_file> -T`
    
    with randomized 3' end test_data generation:
    
    `cotrans_preprocessor -m MAKE_3pEND_TARGETS -A <target_RNA_fasta_file> -T -R -U 30`
    
    This will generate a directory `<target_RNA_name>_targets` that contains:
      * a 3' end targets file `<target_RNA_name>_ends_<n>ntLong.txt`, where <n> is the
        length of the 3' end targets
      * a sub-directory called `intermediate_transcripts` that contains a FASTA file for
        every intermediate transcript sequence
      * a sub-directory called `metrics` that contains target generation QC metrics files
      * if test_data generation was enabled, a sub-directory called test_data that contains
        test_data fastq files
    
3.  Process sequencing reads:

    If analyzing TECprobe-ML data, run `cotrans_preprocessor` in `PROCESS_MULTI` mode:
  
    `cotrans_preprocessor -m PROCESS_MULTI -i <read1_fastq_file> -I <read2_fastq_file> -e <3p_ends_target file>`
  
    If analyzing TECprobe-SL data, run `cotrans_preprocessor` in `PROCESS_SINGLE` mode:
  
    `cotrans_preprocessor -m PROCESS_SINGLE -i <read1_fastq_file> -I <read2_fastq_file> -a <target_RNA_FASTA_file>`
    
    This will generate
      * `fastp` output files `fastp.html` and `fastp.json`
      * `processing.txt`, which contains sequencing read processing metrics
      * `length_distribution.txt` which contains the number of reads that mapped to each transcript length
      * a directory called `split` that contains 
        * fastq files split by untreated/modified channel and transcript length (PROCESS_MULTI only)
        * `smooth_transition.sh` which can be used to concatenate fastq files for neighboring transcripts when performing neighboring transcript smoothing
      * `config.txt` which is used to configure ShapeMapper2 run script generation 
    
4. If performing neighboring transcript smoothing, run the command:

  `sh ./split/smooth_transition.sh`
  
  This will generate the directory `split_smooth` that contains fastq files in which sequencing reads for neighboring transcripts have been concatenated.
  
5. Make a directory for running ShapeMapper2, copy config.txt to that directory, and change to that directory
    
6. Generate a ShapeMapper2 run script:

   After providing all required information in the `config.txt` file, run `cotrans_preprocessor` in `MAKE_RUN_SCRIPT` mode:
   
   `cotrans_preprocessor -m MAKE_RUNSCRIPT -c config.txt`
   
   This will generate 
   * parsed_config.txt, which reports the settings that were parsed from config.txt
   * a shell script containing commands to run ShapeMapper2 analysis for every intermediate transcript
     
7. Run the shell script that was generated in Step 6 to start ShapeMapper2 processing.
  
    

## Assembling TECprobe-ML data using mkmtrx
 
`mkmtrx` assembles data into matrix and other useful formats, reports alignment rates, and can be used to generate rdat files.
  
###mkmtrx inputs and options
  
Required:
```
-m/--mode <run_mode_specifier>  Set run mode. Valid run_mode_specifier values are `MULTI` for
                                TECprobe-ML data and `SINGLE` for TECprobe-SL data
-i/--input <data_directory>     The directory in which the ShapeMapper2 run script was executed
```
 
Optional:
```
-w/--include-up-to <n>          Whitelist transcripts up to length <n> for inclusion in columnized reactivities 
-x/--exclude-term <n>           Exclude transcripts after length <n> from columnized reactivities
```

### Basic usage of mkmtrx

Run the command: 
`mkmtrx -m <run_mode_specifier> -i <data_directory>`
  
This will generate the following files in <data_directory>
  * '<data_directory>_reactivity.csv' which contains a matrix of raw reactivity values
  * '<data_directory>_untreated.csv'  which contains a matrix of untreated effective read depth values
  * '<data_directory>_modified.csv'   which contains a matrix of modified effective read depth values
  * '<data_directory>_alignment_totals.txt' which contains the number of reads that passed processing
    by 'cotrans_preprocessor' and the number of reads that aligned
  * '<data_directory>_alignment_rates.txt' which contains alignment rates for each transcript length
  * '<data_directory>_columns.txt' which contains a column of all reactivity values transcripts that
    were enriched by biotin-streptavidin roadblocks (which can be modified using the `--include-up-to` and
    '--exclude-term` options. The columnized reactivities are useful for plotting replicate correlation.
  * '<data_directory>_linebars.txt' which contains reactivites formatted for making overlapping bar plots

  
## Extract reactivity trajectories using mtrx2cols
    
`mtrx2cols` extracts reactivity trajectories for specific nucleotides from a reactivity matrix.
  
### mtrx2cols inputs and options
  
Required:
```
```
 
Optional:
```
```
  
### Basic usage of mtrx2cols

# TECtools - TECdisplay tools

[**variant_maker**](#generating-variant-targets-using-variant_maker) - Generates a targets file for user-specified sequence variants.

[**TECdisplay_mapper**](#processing-tecdisplay-data-using-tecdisplay_mapper) - Processes and maps sequencing reads to targets.

[**merge_TECdisplay_replicates**](#combining-tecdisplay-replicate-data-using-merge_tecdisplay_replicates) - Merges data from replicate experiments.

[**TECdisplay_navigator**](#filtering-tecdisplay-data-using-tecdisplay_navigator) - Filters data for variants that match user-specified constraints.

[**TECdisplay_Hnav**](#hierarchically-filtering-tecdisplay-data-using-tecdisplay_hnav) - Hierarchically filters data for variants that match a series of user-specified constraints.

[**id2variant**](#reconstructing-a-sequence-from-a-variant-id-using-id2variant) - Reconstructs complete variant sequence from variant id.


## Generating variant targets using variant_maker

`variant_maker` generates a file containing user-defined sequence variants that is used as a targets file by TECdisplay_mapper.

### variant_maker inputs and options

```
-v/--mk-variants <variant_tempate_file>   File containing template for variant sequence construction. 
```

### Basic usage of variant_maker

1. Write a variant template text file that contains the following lines:

```
/name<tab><name_of_sequence>      Name of sequence
/wtsq<tab><wild-type_sequence>    Wild-type sequence
/sxxx<tab><variant_template_seq>  Variant template sequence, 'xxx' is a numerical id
/pxxx<tab><pairing_constraint>    Base pair constraint, 'xxx' must match that of /s line; can supply multiple
#                                 End of variant template indicater
```

The `name` and `wtsq` lines are provided once for the entire file. Each variant template comprises:

- A variant template sequence line specified by `/sxxx`, where `xxx` is a numerical id. The variant template sequence       may use any IUPAC DNA base to specify positions that should be randomized.
     
- One or more pair constraint lines specified by `/pxxx`, where `xxx` is the same numerical id that was provided for the variant template sequence. Valid characters are `.` (no pair), `(` (1st pair partner), and `)` (2nd pair partner). If multiple pairs are supplied in one line, the pairs will close from the inside out. A variant template sequence must match all of its associated pairing constraints to be included in the output file.
    
- A `#` character that indicates the end of the variant template.

Multiple variant templates can be provided in a single file. The output file will contain sequences for every variant template that was specified. Redundant variant sequences that were specified by more than one variant template are not filtered at this stage and are excluded later when the variant template is used by TECdisplay_mapper.

2. Run the command:
`variant_maker -v <variant_templates_file>`

This generates a directory that contains the files:

- <variant_template_name>_input.txt, which contains the name of the variant templates input file

- <variant_template_name>_processing.txt, which contains a record of all processing messages

- <variant_template_name>_variants.txt, which contains a list of all variant sequences in the following format:

  header:
  ```
  variants:<number_of_variants>          number of variant sequences in file
  WT:<sequence_name><tab><wt_sequence>   wild-type name and sequence
  ```

  for each variant template:
  ```
  REF:<variant_template_name>|TPR:<number_of_variants>|VBS:<list_of_variable_bases>|const:<list_of_constant_indels>
  <variant_id><tab><variant_sequence>   line containing variant id and sequence for each variant of the current template
  ```

Variant ids comprise a list of underscore-delimited variable bases using with the following format specifications:
```
  - non-insertion variable base: <position><base>

  - insertion variable base:     <position_of_preceding_non-insertion_nt>i<consecutive_insertion_number><base>
```

Variant ids exclude constant insertions and deletions, which are recorded in the reference line of each variant template unders the specifier `const:`. Constant indels uses the following format:
```
  - constant deletion:  d<position><base>

  - constant insertion: c<position_of_preceding_non-insertion_nt>i<consecutive_insertion_number><base>
```

## Processing TECdisplay data using TECdisplay_mapper

text here



## Combining TECdisplay replicate data using merge_TECdisplay_replicates

text here



## Filtering TECdisplay data using TECdisplay_navigator

text here



## Hierarchically filtering TECdisplay data using TECdisplay_Hnav

text here



## Reconstructing a sequence from a variant id using id2variant

text here



## Processing TECprobe data using cotrans_preprocessor

`cotrans_preprocessor` performs sequencing read preprocessing to prepare data for analysis by `ShapeMapper2`. The expected run time is variable depending on the size of the  data set. A TECprobe-VL data set with a typical sequencing depth (60-100M paired end reads) will typically take 30-60 minutes to process.

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
  * `<data_directory>_reactivity.csv` which contains a matrix of raw reactivity values
  * `<data_directory>_untreated.csv`  which contains a matrix of untreated effective read depth values
  * `<data_directory>_modified.csv`   which contains a matrix of modified effective read depth values
  * `<data_directory>_alignment_totals.txt` which contains the number of reads that passed processing
    by 'cotrans_preprocessor` and the number of reads that aligned
  * `<data_directory>_alignment_rates.txt` which contains alignment rates for each transcript length
  * `<data_directory>_columns.txt` which contains a column of all reactivity values transcripts that
    were enriched by biotin-streptavidin roadblocks (which can be modified using the `--include-up-to` and
    `--exclude-term` options. The columnized reactivities are useful for plotting replicate correlation.
  * `<data_directory>_linebars.txt` which contains reactivites formatted for making overlapping bar plots

  
## Extract reactivity trajectories using mtrx2cols
    
`mtrx2cols` extracts reactivity trajectories for specific nucleotides from one or more reactivity matrices.
  
### mtrx2cols inputs and options
  
Required:
```
-m/--matrix <matrix_csv>  Reactivity matrix in csv format.
-f/--filter <filter_file> Filter file that specifies what nucleotides should be extractd and what 
                          transcript length windows should be output for each nucleotide. 
```
  
The format for the filter file is:

```
line1: 'nts=' followed by a list of comma-separated values specifying which nucleotide trajectores to extract.
  
line2: min=<n>,max=<m>, where <n> and <m> are the start and end of a transcript length window. Multiple transcript
       length window lines can be included.
```
  

For example, the filter
  
```
nts=68,69,70
min=20,max=137
min=169,max=172
  
```

will generate a file containing reactivity trajectories for nucleotides 67, 69, and 70 for windows from transcripts
20 to 137 and 169 to 172.
  
  
Optional:
```
-a/--alias <alias_file>   File containing aliases for input data. The --alias
                          option must be supplied **after** all reactivity matrix
                          csv files
-o/--output <output_name> Output file name
```
  
The format for the alias file is:
  
`<input_matrix_name><tab><alias>`
  
Aliases for multiple reactivity matrices can be supplied in the same alias file if more than one reactivity matrix csv is being processed.
  
### Basic usage of mtrx2cols
  
`mtrx2cols -m <matrix_csv1> -m <matrix_csv2> -a <aliases> -f <filter> -o <output_file_name>`
  
will generate a file `<output_file_name>.txt` that contains reactivity trajectories for the nucleotide specified by `<filter>`, with column names that use the aliases supplied by `<aliases>`.

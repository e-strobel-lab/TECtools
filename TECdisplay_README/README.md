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
-v/--mk-variants <variant_template_file>   File containing template for variant sequence construction. 
```

### Basic usage of variant_maker

1. Write a variant template text file that contains the following lines:

```
/name<tab><name_of_sequence>      Name of sequence
/wtsq<tab><wild-type_sequence>    Wild-type sequence
/sxxx<tab><variant_template_seq>  Variant template sequence, 'xxx' is a numerical id
/pxxx<tab><pairing_constraint>    Base pair constraint, 'xxx' must match that of /s line; can supply multiple pairing constraints
#                                 End of variant template indicater
```

The `name` and `wtsq` lines are provided once for the entire file. Each variant template comprises:

- A variant template sequence line specified by `/sxxx`, where `xxx` is a numerical id. The variant template sequence may use any IUPAC DNA base to specify positions that should be randomized and `-` to specify constant deletions. **The variant template sequence must align with with the wild-type sequence specified by `/wtsq`.** Insertions in the variant template sequence can be accommodated by including corresponding `.` characters in the wild-type sequence.
     
- One or more pairing constraint lines specified by `/pxxx`, where `xxx` is the same numerical id that was provided for the associated variant template sequence. Valid characters are `.` (no pair), `(` (1st pair partner), and `)` (2nd pair partner). If multiple pairs are supplied in one line, the pairs will close from the inside out. A variant template sequence must match all of its associated pairing constraints to be included in the output file.
    
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
  - variable non-insertion:        <position><base>

  - variable insertion (pos >= 1): <position_of_preceding_non-insertion_nt>i<consecutive_insertion_number><base>

  - variable insertion (pos < 1):  <position><base>
```

Variant ids exclude constant insertions and deletions, which are recorded in the reference line of each variant template unders the specifier `const:`. Constant indels uses the following format:
```
  - constant insertion (pos >= 1): c<position_of_preceding_non-insertion_nt>i<consecutive_insertion_number><base>

  - constant insertion (pos < 1):  c<position><base>

  - constant deletion:             d<position><base>
```



## Processing TECdisplay data using TECdisplay_mapper

`TECdisplay_mapper` coordinates sequencing read processing by `fastp` and maps sequencing reads to targets that were generated by `variant_maker`.  The expected run time is variable depending on the size of the  data set. A TECprobe-VL data set with a typical sequencing depth (200-300M paired end reads) will typically take 30-60 minutes to process.

**`TECdisplay_mapper` requires that `fastp` (https://github.com/OpenGene/fastp) is installed**

### TECdisplay_mapper inputs and options

```
-m/--mode <run_mode_specifier>     set run mode using one of following run mode specifiers:
                                   MAP_SEQ_READS  Process and map sequencing reads from a TECdisplay experiment.
                                   MAP_TEST_DATA  Generate and map test data for an input targets file.

-i/--read1 <read_1_fastq_file>     Read 1 FASTQ file input.
                                   Required for MAP_SEQ_READS mode.
                                   Must not be provided in MAP_TEST_DATA mode.

-i/--read2 <read_2_fastq_file>     Read 2 FASTQ file input.
                                   Required for MAP_SEQ_READS mode.
                                   Must not be provided in MAP_TEST_DATA mode.

-t/--targets <targets_file>        Targets file input (required). Targets file must be generated by variant_maker.

-o/--out_name <output_file_name>   Output file name (optional)

-p/--fastp_path <path_to_fastp>    Path to fastp executable (optional). Only needed if fastp was not added to PATH.

-v/--qual-variable <value 0 to 41> Minimum qscore for variable bases (optional). Default is 20.

-c/--qual-constant <value 0 to 41> Minimum qscore for constant bases (optional). Default is 0.
```

### Basic usage of TECdisplay_mapper

If analyzing TECdisplay data, run `TECdisplay_mapper` in `MAP_SEQ_READS` mode:

`TECdisplay_mapper -m MAP_SEQ_READS -i <read1_fastq> -I <read2_fastq> -t <targets>`

If generating and mapping test data, run `TECdisplay_mapper` in `MAP_TEST_DATA` mode:

`TECdisplay_mapper -m MAP_TEST_DATA -t <targets>`

In both `MAP_SEQ_READS` and `MAP_TEST_DATA` mode, the following files/directories will be generated:

- out.txt (or <out_name>_out.txt if -o option was used), which contains the number of reads that mapped to the bound and unbound channels for each variant and the fraction of reads tha mapped to the bound channel.
  
- mapping.txt, which contains a list of the input files that were provided and report of read mapping metrics.

- processed, which contains the output fastq files from fastp processing

- navigator_templates, which contains a template for use with `TECdisplay_navigator` and `TECdisplay_Hnav`.
  
- fastp_command.txt, which records the fastp command that was used

- fastp.html/fastp.json, fastp processing metrics


In `MAP_TEST_DATA` mode, the following additional files will be generated:

- test_data_analysis.tx, test data analysis metrics

- <targets_file_name>_test_data, which contains the test_data fastq files that were used for the test data analysis. 


## Combining TECdisplay replicate data using merge_TECdisplay_replicates

text here



## Filtering TECdisplay data using TECdisplay_navigator

`TECdisplay_navigator` filters TECdisplay data by user-supplied constraints on nucleotide identity and base pairing. 

### TECdisplay_navigator inputs and options

```
-v/--values <values_file_input>            Input values file (required). Multiple input values files can be supplied by
                                           providing the -v option more than once. Values files are typically the output
                                           of TECdisplay_mapper, but can be any tab-delimited file in which the first line
                                           contains column headears and the first column contains variant ids in the format                                                      used by variant_maker as long as the -n/--non-standard option is used.

-c/--constraints <constraints_file_input>  Constraints file (required), containing one or more constraints (format described below)
                                           that will be used to filter the TECdisplay data. When filtering for matches to the
                                           constraints (default), each constraint is handled separately and one output file, which
                                           contains matches, is generated for each constraint.

-x/--exclude                               Flag to exclude constraint matches (optional). Setting this flag filters for variants that
                                           do not match any of the constraints in the constraints file.

-o/--out-name <output_directory_name>      Output directory name (required).

-f/--out-prefix <output_file_prefix>       Output file prefix (required). Prefix to append to all output files.

-n/--non-standard                          Flag that input files are a non-standard format. Non-standard values files must be a
                                           tab-delimited file in which the first line contains column headears and the first column                                              contains variant ids in the format used by variant_maker.
```

`TECdisplay_mapper` generates a template constraint file that can be found in the directory `navigator_templates` in the `TECdisplay_mapper` output directory. 

Constraints files specify variable and constant bases in the format described above for variant_maker, which is provided here for convenience:

```
  - variable non-insertion:        <position><base>
  - constant non-insertion:        c<position><base>
  - variable insertion (pos >= 1): <position_of_preceding_non-insertion_nt>i<consecutive_insertion_number><base>
  - variable insertion (pos < 1):  <position><base>
  - constant insertion (pos >= 1): c<position_of_preceding_non-insertion_nt>i<consecutive_insertion_number><base>
  - constant insertion (pos < 1):  c<position><base>
  - constant deletion:             d<position><base>
```

All constraints files contain a header that links them to the targets that were used by TECdisplay_mapper to map the reads:

```
/seq<tab><wild_type_sequence>
/vbs<tab><variant_template_sequence>
/constant_indels:<list_of_constant_indels>
```

Following this header, a series of filters can be supplied using the format:

```
<filter_name>
base<tab><variable_base>     //Every variable base in the sequence is listed using the specifier
base<tab><variable_base>     //'base' and a variable base in the format defined above. By default,
base<tab><variable_base>     //each 'base' line contains the variable base that was specified in the
base<tab><variable_base>     //variant template sequence. The nucleotide identity of each base can be
base<tab><variable_base>     //changed to filter for variants that match the specified base.
base<tab><variable_base>     
pair<tab><base1>,<base2><tab><pair_type> //Pair constraints are optional. Each base in the pair must be
pair<tab><base1>,<base2><tab><pair_type> //provided in the format defined above. If the base is variable,
pair<tab><base1>,<base2><tab><pair_type> //it must match the corresponding base constraint in the lines
pair<tab><base1>,<base2><tab><pair_type> //above. All valid pair_type specifiers are defined below.
```

The following `pair_type` specifiers are valid:

```
ANY_PAIR      GC, AU, GU
WC_PAIR       GC, AU  
STRONG        GC
WEAK          AU, GU
WEAK_AU       AU
WEAK_GU       GU    
MISMATCH      Must be mismatch     
NO_CONSTRAINT Pair is not constrained (used as placeholder)
```

### Basic usage of TECdisplay_navigator

If filtering for matches to the supplied constraints, run `TECdisplay_navigator` using the command:

`TECdisplay_navigator -v <input_values_file> -c <constraints_file> -o <output_directory_name> -f <output_file_prefix>`

If filtering to variants that do not match the supplied constraints, run `TECdisplay_navigator` using the command:

`TECdisplay_navigator -v <input_values_file> -x <constraints_file> -o <output_directory_name> -f <output_file_prefix>`


## Hierarchically filtering TECdisplay data using TECdisplay_Hnav

`TECdisplay_Hnav` hierarchichally filters TECdisplay data by user-supplied constraints on nucleotide identity and base pairing by running TECdisplay_navigator sequentially on multiple user-supplied constraints. 

### TECdisplay_Hnav inputs and options

```
-v/--values <values_file_input>            Input values file (required). Multiple input values files can be supplied by
                                           providing the -v option more than once. Values files are typically the output
                                           of TECdisplay_mapper, but can be any tab-delimited file in which the first line
                                           contains column headears and the first column contains variant ids in the format                                                               used by variant_maker as long as the -n/--non-standard option is used.

-c/--constraints <constraints_file_input>  Inclusion constraints file input, containing one or more constraints (format                                                                   described above for TECdisplay_navigator). Constraints files supplied using the -c 
                                           option will be used for filter to matches for each constraint. As for TECdisplay_navigator
                                           each constraint is handled separately and one output file, which contains matches,
                                           is generated for each constraint. Multiple constraints files can be supplied. Data
                                           will be filtered by the constraints files in the order that they are supplied.

-x/--exclude <constraints_file_input>      Exclusion constraints file input, containing one or more constraints. Constraints
                                           files supplied using the -x option will be used to filter for variants that do not match                                                       any of the constraints in the exclusion constraitls file.
                                           

-o/--out-name <output_directory_name>      Output directory name (required).

-f/--out-prefix <output_file_prefix>       Output file prefix (required). Prefix to append to all output files.

-p/--path <TECdisplay_navigator_path>      Path to TECdisplay_navigator executable (optional). Only needed if TECdisplay_navigator
                                           has not been added to PATH.

-a/--aggregate                             Generate file in which the fraction bound columns of filtered variant output files are
                                           aggregated into a single file. This option can only be used when standard TECdisplay_mapper 
                                           output files are used as the input values files.
  
```

### Basic usage of TECdisplay_navigator

To hierarchically filter TECdisplay data using `TECdisplay_Hnav`, run the command:

`TECdisplay_Hnav -v <input_values_file> -c <constraints_file_1> -c <constraints_file_2> -c <constraints_file_3> -o <output_directory_name> -f <output_file_prefix>`

The command above will filter the input values file by each of the three constraints sequentially in the order that they are supplied. In each layer of filtering, every output file from the layer is filtered by the current constraint.  The output of each filtering layer is stored in a separate directory. **This can potentially lead to the generation of large numbers of files**. To safeguard against unintentionally generating large numbers of files, a prompt will indicate the number of files that will be created and ask if this is acceptable. Answering 'yes' will allow the analysis to proceed.

## Reconstructing a sequence from a variant id using id2variant

text here




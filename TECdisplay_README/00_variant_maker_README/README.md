## Generating variant targets using variant_maker

`variant_maker` generates a file containing user-defined sequence variants that is used as a targets file by `TECdisplay_mapper`.

### variant_maker inputs and options

```
-v/--mk-variants <variant_template_file>   File containing template for variant sequence construction (required). 
```

The format of variant template files is described below, and an example variant template file is provided in the `example_files` directory.

### Basic usage of variant_maker

1. Write a variant template text file that contains the following lines:

```
/name<tab><name_of_sequence>      Name of sequence
/wtsq<tab><wild-type_sequence>    Wild-type sequence
/sxxx<tab><variant_template_seq>  Variant template sequence, 'xxx' is a numerical id
/pxxx<tab><pairing_constraint>    Base pair constraint, 'xxx' must match that of the /s line; can supply multiple pairing constraints
#                                 End of variant template indicator
```

The `name` and `wtsq` lines are provided once for the entire file. Each variant template comprises:

- A variant template sequence line specified by `/sxxx`, where `xxx` is a numerical id. The variant template sequence may use any IUPAC DNA base to specify positions that should be randomized and `-` to specify constant deletions. **The variant template sequence must align with the wild-type sequence specified by `/wtsq`.** Insertions in the variant template sequence can be accommodated by including corresponding `.` characters in the wild-type sequence.
     
- One or more pairing constraint lines specified by `/pxxx`, where `xxx` is the same numerical id that was provided for the associated variant template sequence. Valid characters are `.` (no pair), `(` (1st pair partner), and `)` (2nd pair partner). If multiple pairs are supplied in one line, the pairs will close from the inside out. A variant template sequence must match all of its associated pairing constraints to be included in the output file.
    
- A `#` character that indicates the end of the variant template.

Multiple variant templates can be provided in a single file. The output file will contain sequences for every variant template that was specified. Redundant variant sequences that were specified by more than one variant template are not filtered at this stage and are excluded later when the variant template is used by `TECdisplay_mapper`.

2. Run the command:
`variant_maker -v <variant_template_file>`

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

Variant ids comprise a list of underscore-delimited variable bases with the following format specifications:
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

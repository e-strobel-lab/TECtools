## Filtering TECdisplay data using TECdisplay_navigator

`TECdisplay_navigator` filters TECdisplay data by user-supplied constraints on nucleotide identity and base pairing. 

### TECdisplay_navigator inputs and options

```
-v/--values <values_file_input>            Input values file (required). Multiple input values files can be supplied by
                                           providing the -v option more than once. Values files are typically the output
                                           of TECdisplay_mapper, but can be any tab-delimited file in which the first line
                                           contains column names and the first column contains variant ids in the format
                                           used by variant_maker as long as the -n/--non-standard option is used.

-c/--constraints <constraints_file_input>  Constraints file containing one or more constraints (format described below)
                                           that will be used to filter the TECdisplay data (required). When filtering for matches to the
                                           constraints (default), each constraint is handled separately and one output file, which
                                           contains matches, is generated for each constraint.

-x/--exclude                               Flag to exclude constraint matches (optional). Setting this flag filters for variants that
                                           do not match any of the constraints in the constraints file.

-o/--out-name <output_directory_name>      Output directory name (optional).

-f/--out-prefix <output_file_prefix>       Output file prefix (optional). Prefix to append to all output files.

-n/--non-standard                          Flag that input files are a non-standard format. Non-standard values files must be a
                                           tab-delimited file in which the first line contains column names and the first column
                                           contains variant ids in the format used by variant_maker.
```

`TECdisplay_mapper` generates a template constraint file that can be found in the directory `navigator_templates` in the `TECdisplay_mapper` output directory. The format of constraints files is described below, and an example constraint file is provided in the `example_files` directory.

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
base<tab><variable_base>     //changed to filter for variants that match the specified base. All variable
base<tab><variable_base>     //bases must be listed as base constraints even if they are not being constrained.
base<tab><variable_base>
pair<tab><base1>,<base2><tab><pair_type> //Pair constraints are optional. Each base in the pair must be
pair<tab><base1>,<base2><tab><pair_type> //provided in the format defined above. If the base is variable,
pair<tab><base1>,<base2><tab><pair_type> //it must match the corresponding base constraint in the lines
pair<tab><base1>,<base2><tab><pair_type> //above. All valid pair_type specifiers are defined below.
```

The following `pair_type` specifiers are valid:

```
Specifier     Allowed base pairs
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

If filtering for variants that do not match the supplied constraints, run `TECdisplay_navigator` using the command:

`TECdisplay_navigator -v <input_values_file> -c <constraints_file> -x -o <output_directory_name> -f <output_file_prefix>`

## Hierarchically filtering TECdisplay data using TECdisplay_Hnav

`TECdisplay_Hnav` hierarchichally filters TECdisplay data by user-supplied constraints on nucleotide identity and base pairing by running TECdisplay_navigator sequentially on multiple user-supplied constraints. 

### TECdisplay_Hnav inputs and options

```
-v/--values <values_file_input>            Input values file (required). Multiple input values files can be supplied by
                                           providing the -v option more than once. Values files are typically the output
                                           of TECdisplay_mapper, but can be any tab-delimited file in which the first line
                                           contains column names and the first column contains variant ids in the format
                                           used by variant_maker as long as the -n/--non-standard option is used.

-c/--constraints <constraints_file_input>  Relative path to an inclusion constraints file, containing one or more constraints
                                           (format described above for TECdisplay_navigator). Constraints files supplied using
                                           the -c option will be used to filter for matches to each constraint. As for
                                           TECdisplay_navigator, each constraint is handled separately and one output file,
                                           which contains matches, is generated for each constraint. Multiple constraints
                                           files can be supplied. Data will be filtered by the constraints files in the order
                                           that they are supplied.

-x/--exclude <constraints_file_input>      Relative path to an exclusion constraints file containing one or more constraints.
                                           Constraints files supplied using the -x option will be used to filter for variants
                                           that do not match any of the constraints in the exclusion constraints file. Multiple
                                           constraints files can be supplied. Data will be filtered by the constraints files in
                                           the order that they are supplied.

-o/--out-name <output_directory_name>      Output directory name (optional).

-f/--out-prefix <output_file_prefix>       Output file prefix (required). Prefix to append to all output files.

-p/--path <TECdisplay_navigator_path>      Absolute path to TECdisplay_navigator executable (optional). This must be the absolute
                                           path; a relative path is not valid. Only needed if TECdisplay_navigator has not been
                                           added to PATH. 

-a/--aggregate                             Generate file in which the fraction bound columns of filtered variant output files are
                                           aggregated into a single file (optional). This option can only be used when standard
                                           TECdisplay_mapper output files are used as the input values files.
  
```

### Basic usage of TECdisplay_Hnav

To hierarchically filter TECdisplay data using `TECdisplay_Hnav`, run the command:

`TECdisplay_Hnav -p <path_to_TECdisplay_navigator> -v <input_values_file> -c <constraints_file_1> -c <constraints_file_2> -c <constraints_file_3> -o <output_directory_name> -f <output_file_prefix>`

The command above will filter the input values file by each of the three constraints sequentially in the order that they are supplied. In each layer of filtering, every output file from the previous layer is filtered by the current constraint.  The output of each filtering layer is stored in a separate directory. **This can potentially lead to the generation of large numbers of files**. To safeguard against unintentionally generating large numbers of files, a prompt will indicate the number of files that will be created and ask if this is acceptable. Answering 'yes' will allow the analysis to proceed.

### Running TECdisplay_Hnav using the example data and constraints files
'TECdisplay_Hnav' can be run using the provided example data set and constraints file using the command: 

`TECdisplay_Hnav -p <path_to_TECdisplay_navigator> -v <path_to_01_pflZTP_U17_doseCurve_example_data.txt> -c <path_to_02_TECdisplay_nav_constraint_example_U17.txt> -o <example_Hnav_out> -f <pflZTP_U17>`

This command will generate the output directory 'example_Hnav_out'. The subdirectory 'layer_1' contains dose curves for the wild type, flip, strong, and weak pseudoknot variants of the pfl ZTP riboswitch Ultramer 17 TECdisplay library.


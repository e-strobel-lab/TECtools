## Assembling TECprobe-VL data using mkmtrx
 
`mkmtrx` assembles data into matrix and other useful formats, reports alignment rates, and can be used to generate rdat files.
  
### mkmtrx inputs and options
  
Required:
```
-m/--mode <run_mode_specifier>  Set run mode. Valid run_mode_specifier values are `MULTI` for
                                TECprobe-VL data and `SINGLE` for TECprobe-SL data
-i/--input <data_directory>     The directory in which the ShapeMapper2 run script was executed
-c/--reactivity-col <input>     Set reactivity value to output. Valid inputs are `REACTIVITY_PROFILE`
                                for unfilitered raw reactivity, `HQ_PROFILE` for high quality nucleotide
                                raw reactivity, and `NORM_PROFILE` for normalized reactivity.
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

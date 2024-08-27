## Combining TECdisplay replicate data using merge_TECdisplay_replicates

`merge_TECdisplay_replicates` combines the reads from replicate TECdisplay data sets and calculates a new fraction bound value.

### merge_TECdisplay_replicates inputs and options

```
-v/--values <values_file_input>   Input values file (required). Must be output file from TECdisplay_mapper.
                                  The -v option must be supplied for each input values file.

-o/--out-name <output_file_name>  Output file name.
```

### Basic usage of merge_TECdisplay_replicates

To merged replicate TECdisplay data using `merge_TECdisplay_replicates`, run the command:

`merge_TECdisplay_replicates -v <replicate_1_data> -v <replicate_2_data> -o <output_file_name>`

This will generate the files:

- <output_file_name>_merged.txt, which contains the merged TECdisplay data.
  
- <output_file_name>_merged_replicate_merge_record.txt, which records the names of the merged input files and the output file.

## Normalizing TECprobe-VL data using process_TECprobeVL_profiles  

`process_TECprobeVL_profiles` performs whole-dataset normalization for TECprobe-VL data (which is normalized on a per-transcript length basis during analysis) and assembles ShapeMapper2 (https://github.com/Weeks-UNC/shapemapper2) output file data into a single csv file that is compatible with TECprobe visualization tools. Normalization is performed exactly as in ShapeMapper2:

```
"Over the set of all RNAs and nucleotide positions without masked (lowercase) sequence, high background, or
low read depth, reactivities are normalized by dividing by the mean reactivity of the top 10% of reactivities
after reactivities above a threshold are excluded (see section 3.2.1 in Low and Weeks, 2010). That threshold
is selected from the largest value out of [1.5 × interquartile range, 90th percentile (if total seq length > 100)
or 95th percentile (if total seq length < 100)]." - ShapeMapper2 README
```

If more than one input directory is provided, `process_TECprobeVL_profiles` will merge the data into a single dataset. Therefore, only replicate data sets should be provided together. `process_TECprobeVL_profiles` checks the attributes of input data names and the data itself to confirm that input datasets are compatible.

### process_TECprobeVL_profiles inputs and options
  
Required:
```
-i/--input <data_directory>     Input TECprobe-VL data directory. This should be the directory in which the
                                ShapeMapper2 run script was executed. If more than one compatible input directory
                                is provided, the data will be merged.
```
 
Optional:
```
-o/--out_dir_name <input>       Output directory name. Default = 'dataset_norm_out'
-n/--sample-name <input>        Sample name. If no sample name is provided, a sample name will be automatically
                                generated from the input data names.
-e/--min-depth <n>              Minimum read depth for high-quality nucleotides. Default = 5000.
-b/--max-background <n>         Maximum background mutation rate for high-quality nucleotides. Default = 0.05.
```

### Basic usage of process_TECprobeVL_profiles

To process a single TECprobe-VL dataset, run the command: 
`process_TECprobeVL_profiles -i <data_directory>`

To merge three TECprobe-VL datasets, run the command: 
`process_TECprobeVL_profiles -i <data_directory> -i <data_directory> -i <data_directory>`

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

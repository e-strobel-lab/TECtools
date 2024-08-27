## Calculating the difference in fraction bound values using calc_FracBound_Difference

`calc_FracBound_difference` calculates the difference in fraction bound values for two TECdisplay datasets. Checks are performed to ensure that the two TECdisplay input datasets are compatible.

### calc_FracBound_difference inputs and options
  
Required:
```
-m/--minuend <TECdisplay data file>     minuend file input
-s/--subtrahend <TECdisplay data file>  subtrahend file input
-o/--out-name <output file name>        output file name
```

### Basic usage of calc_FracBound_difference

To calculate the difference in fraction bound values by subtracting `<input2>` from `<input1>`, run the command:
`calc_FracBound_difference -m <input1> -s <Input2> -o <output file name>`

This will generate the file `<output file name>`, which contains the TECdisplay data from both input files and the additional column "diff", which contains the difference values that were calculated.

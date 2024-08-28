## Reconstructing a sequence from a variant id using id2variant

`id2variant` reconstructs complete variant sequences from variant ids.

### id2variant inputs and options

```
-v/--variants <variants_input_file>    Variants input file (required). Format defined below.

-o/--out-name <output_directory_name>  Output directory name (optional).

-m/--min_nucleotide <value>            Minimum nucleotide to print (optional).

-x/--max_nucleotide <value>            Maximum nucleotide to print (optional).
```

An example variants input file is provided in the directory `example_files`. The variants input file must contain the following header, which can be copied from the `TECdisplay_navigator` constraint template that is generated by `TECdisplay_mapper`:

```
/seq<tab><wild_type_sequence>
/vbs<tab><variant_template_sequence>
/constant_indels:<list_of_constant_indels>
```

Following this header, the variant ids to be reconstructed are supplied on separate lines. If desired, it is possible to include an alias for each variant to simplify their identification by using the format:

`<variant_id><tab><alias>`

### Basic usage of id2variant

To reconstruct variant sequences from variant ids using `id2variant`, run the command:

`id2variant -v <variants_input_file> -o <output_directory_name>`

This will generate a directory that contains both individual FASTA files for each variant and an aggregate FASTA file that contains every variant.

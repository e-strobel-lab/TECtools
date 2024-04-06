#!/bin/sh
#
#  build_TECprobe.sh
#  
#
#  Created by Eric Strobel on 2/7/23.
#  
#
#build TECdisplay_mapper
gcc ../TECtools/internals/TECdisplay_mapper/TECdisplay_mapper_main.c ../TECtools/internals/utils/debug.c ../TECtools/internals/utils/gen_utils.c ../TECtools/internals/utils/io_management.c ../TECtools/internals/seq_utils/isDNAbase.c ../TECtools/internals/seq_utils/isIUPACbase.c ../TECtools/internals/seq_utils/revcomp.c ../TECtools/internals/seq_utils/seq2bin_hash.c ../TECtools/internals/TECdisplay_mapper/map_reads/UNV/call_fastp.c ../TECtools/internals/TECdisplay_mapper/TECdisplay_output_column_headers.c ../TECtools/internals/TECdisplay_mapper/map_reads/map_reads.c ../TECtools/internals/TECdisplay_mapper/map_reads/map_expected/parse_mx_trgts.c ../TECtools/internals/TECdisplay_mapper/map_reads/map_expected/mk_output_files.c ../TECtools/internals/TECdisplay_mapper/map_reads/UNV/prcs_chnl.c ../TECtools/internals/TECdisplay_mapper/testdata_analysis/mk_test_data.c ../TECtools/internals/TECdisplay_mapper/testdata_analysis/assess_test_data.c -lm -o TECdisplay_mapper
#
#
#
#
#build variant_maker
gcc ../TECtools/internals/variant_maker/variant_maker.c ../TECtools/internals/utils/gen_utils.c ../TECtools/internals/utils/io_management.c ../TECtools/internals/seq_utils/ispair.c ../TECtools/internals/seq_utils/isIUPACbase.c ../TECtools/internals/seq_utils/isDNAbase.c ../TECtools/internals/seq_utils/seq2bin_hash.c ../TECtools/internals/seq_utils/seq2bin_long.c ../TECtools/internals/seq_utils/test_possible_pairs.c ../TECtools/internals/seq_utils/basemap.c  ../TECtools/internals/variant_maker/count_variants.c ../TECtools/internals/variant_maker/read_varFile.c  ../TECtools/internals/variant_maker/expand_variant_template.c ../TECtools/internals/variant_maker/filter_variants.c ../TECtools/internals/variant_maker/mk_variant_nm.c  ../TECtools/internals/variant_maker/make_barcodes.c  -lm -o variant_maker
#
#
#
#
#build TECDisplay_navigator
gcc ../TECtools/internals/TECdisplay_navigator/TECdisplay_navigator.c ../TECtools/internals/utils/io_management.c ../TECtools/internals/utils/debug.c ../TECtools/internals/seq_utils/isDNAbase.c ../TECtools/internals/seq_utils/isIUPACbase.c ../TECtools/internals/seq_utils/ispair.c ../TECtools/internals/seq_utils/is_dgnrt_mtch.c ../TECtools/internals/seq_utils/test_possible_pairs.c ../TECtools/internals/seq_utils/basemap.c ../TECtools/internals/TECdisplay_mapper/TECdisplay_output_column_headers.c ../TECtools/internals/TECdisplay_navigator/merge_values_files.c ../TECtools/internals/TECdisplay_navigator/parse_TECdisplay_out_line.c ../TECtools/internals/TECdisplay_navigator/parse_reference.c ../TECtools/internals/TECdisplay_navigator/parse_constraints.c ../TECtools/internals/TECdisplay_navigator/read_vbase.c ../TECtools/internals/TECdisplay_navigator/search_4_vbase_match.c ../TECtools/internals/TECdisplay_navigator/filter_values.c  -lm -o TECdisplay_navigator
#
#
#
#
#build merge_TECdisplay_replicates
gcc ../TECtools/internals/merge_TECdisplay_replicates/merge_TECdisplay_replicates.c ../TECtools/internals/utils/io_management.c ../TECtools/internals/utils/debug.c ../TECtools/internals/utils/gen_utils.c ../TECtools/internals/TECdisplay_mapper/TECdisplay_output_column_headers.c ../TECtools/internals/TECdisplay_navigator/parse_TECdisplay_out_line.c ../TECtools/internals/TECdisplay_mapper/map_reads/map_expected/mk_output_files.c -lm -o merge_TECdisplay_replicates
#
#
#
#
#build calc_FracBound_difference
gcc ../TECtools/internals/calc_FracBound_difference/calc_FracBound_difference.c ../TECtools/internals/utils/io_management.c ../TECtools/internals/utils/debug.c ../TECtools/internals/utils/gen_utils.c ../TECtools/internals/TECdisplay_mapper/TECdisplay_output_column_headers.c ../TECtools/internals/TECdisplay_navigator/parse_TECdisplay_out_line.c ../TECtools/internals/TECdisplay_mapper/map_reads/map_expected/mk_output_files.c -lm -o calc_FracBound_difference
#
#
#
#
#build id2variant
gcc ../TECtools/internals/id2variant/id2variant.c ../TECtools/internals/utils/io_management.c  ../TECtools/internals/seq_utils/isDNAbase.c ../TECtools/internals/seq_utils/isIUPACbase.c  ../TECtools/internals/seq_utils/is_dgnrt_mtch.c ../TECtools/internals/seq_utils/test_possible_pairs.c ../TECtools/internals/seq_utils/basemap.c ../TECtools/internals/TECdisplay_navigator/parse_reference.c ../TECtools/internals/TECdisplay_navigator/read_vbase.c -o id2variant
#
#
#
#
#build TECDisplay_Hnav
gcc ../TECtools/internals/TECdisplay_Hnav/TECdisplay_Hnav.c ../TECtools/internals/utils/io_management.c ../TECtools/internals/utils/debug.c ../TECtools/internals/seq_utils/basemap.c ../TECtools/internals/seq_utils/isDNAbase.c ../TECtools/internals/seq_utils/isIUPACbase.c ../TECtools/internals/seq_utils/is_dgnrt_mtch.c  ../TECtools/internals/seq_utils/test_possible_pairs.c ../TECtools/internals/TECdisplay_mapper/TECdisplay_output_column_headers.c ../TECtools/internals/TECdisplay_navigator/merge_values_files.c ../TECtools/internals/TECdisplay_navigator/parse_TECdisplay_out_line.c ../TECtools/internals/TECdisplay_navigator/parse_reference.c ../TECtools/internals/TECdisplay_navigator/parse_constraints.c  ../TECtools/internals/TECdisplay_navigator/read_vbase.c ../TECtools/internals/TECdisplay_navigator/search_4_vbase_match.c ../TECtools/internals/TECdisplay_Hnav/TECdisplay_Hnav_global_vars.c ../TECtools/internals/TECdisplay_Hnav/get_constraint_metadata.c ../TECtools/internals/TECdisplay_Hnav/calc_output_files.c ../TECtools/internals/TECdisplay_Hnav/Hfilter.c ../TECtools/internals/TECdisplay_Hnav/process_output_files.c  -lm -o TECdisplay_Hnav

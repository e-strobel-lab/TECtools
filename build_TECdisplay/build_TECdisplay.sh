#!/bin/sh
#
#  build_TECprobe.sh
#  
#
#  Created by Eric Strobel on 2/7/23.
#  
#
#build TECdisplay_mapper
gcc ../internals/TECdisplay_mapper/TECdisplay_mapper_main.c ../internals/utils/debug.c ../internals/utils/gen_utils.c ../internals/utils/io_management.c ../internals/seq_utils/isDNAbase.c ../internals/seq_utils/isIUPACbase.c ../internals/seq_utils/revcomp.c ../internals/seq_utils/seq2bin_hash.c ../internals/seq_utils/seq2bin_long.c ../internals/seq_utils/mapping_metrics.c ../internals/variant_maker/constant_seqs.c ../internals/variant_maker/vmt_suffix.c  ../internals/TECdisplay_mapper/map_reads/UNV/call_fastp_TDSPLY.c ../internals/TECdisplay_mapper/map_reads/UNV/prcs_chnl_TDSPLY.c ../internals/TECdisplay_mapper/TECdisplay_output_column_headers.c ../internals/TECdisplay_mapper/map_reads/map_reads.c  ../internals/TECdisplay_mapper/map_reads/map_expected/set_barcoded_compact_target.c  ../internals/TECdisplay_mapper/map_reads/map_expected/get_key.c ../internals/TECdisplay_mapper/map_reads/map_expected/set_trgt.c ../internals/TECdisplay_mapper/map_reads/map_expected/parse_vmt_trgts.c ../internals/TECdisplay_mapper/map_reads/map_expected/mk_output_files.c ../internals/TECdisplay_mapper/map_reads/map_expected/print_navigator_template.c  ../internals/TECdisplay_mapper/testdata_analysis/mk_TDSPLY_test_data.c ../internals/TECdisplay_mapper/testdata_analysis/assess_TDSPLY_test_data.c -lm -o TECdisplay_mapper
#
#
#
#
#build variant_maker
gcc ../internals/variant_maker/variant_maker.c ../internals/utils/debug.c ../internals/utils/gen_utils.c ../internals/utils/io_management.c ../internals/seq_utils/ispair.c ../internals/seq_utils/isIUPACbase.c ../internals/seq_utils/isDNAbase.c ../internals/seq_utils/seq2bin_hash.c ../internals/seq_utils/seq2bin_long.c ../internals/seq_utils/test_possible_pairs.c ../internals/seq_utils/basemap.c ../internals/variant_maker/constant_seqs.c ../internals/variant_maker/vmt_suffix.c ../internals/variant_maker/count_variants.c ../internals/variant_maker/read_varFile.c  ../internals/variant_maker/expand_variant_template.c ../internals/variant_maker/filter_variants.c ../internals/variant_maker/mk_variant_nm.c  ../internals/variant_maker/make_barcodes.c ../internals/variant_maker/read_bcFile.c ../internals/variant_maker/print_output.c -lm -o variant_maker
#
#
#
#
#build TECDisplay_navigator
gcc ../internals/TECdisplay_navigator/TECdisplay_navigator.c ../internals/utils/io_management.c ../internals/utils/debug.c ../internals/seq_utils/isDNAbase.c ../internals/seq_utils/isIUPACbase.c ../internals/seq_utils/ispair.c ../internals/seq_utils/is_dgnrt_mtch.c ../internals/seq_utils/test_possible_pairs.c ../internals/seq_utils/basemap.c ../internals/TECdisplay_mapper/TECdisplay_output_column_headers.c ../internals/TECdisplay_navigator/merge_values_files.c ../internals/TECdisplay_navigator/parse_TECdisplay_out_line.c ../internals/TECdisplay_navigator/parse_reference.c ../internals/TECdisplay_navigator/parse_constraints.c ../internals/TECdisplay_navigator/read_vbase.c ../internals/TECdisplay_navigator/search_4_vbase_match.c ../internals/TECdisplay_navigator/filter_values.c  -lm -o TECdisplay_navigator
#
#
#
#
#build merge_TECdisplay_replicates
gcc ../internals/merge_TECdisplay_replicates/merge_TECdisplay_replicates.c ../internals/utils/io_management.c ../internals/utils/debug.c ../internals/utils/gen_utils.c ../internals/TECdisplay_mapper/TECdisplay_output_column_headers.c ../internals/TECdisplay_navigator/parse_TECdisplay_out_line.c ../internals/TECdisplay_mapper/map_reads/map_expected/mk_output_files.c -lm -o merge_TECdisplay_replicates
#
#
#
#
#build calc_FracBound_difference
gcc ../internals/calc_FracBound_difference/calc_FracBound_difference.c ../internals/utils/io_management.c ../internals/utils/debug.c ../internals/utils/gen_utils.c ../internals/TECdisplay_mapper/TECdisplay_output_column_headers.c ../internals/TECdisplay_navigator/parse_TECdisplay_out_line.c ../internals/TECdisplay_mapper/map_reads/map_expected/mk_output_files.c -lm -o calc_FracBound_difference
#
#
#
#
#build id2variant
gcc ../internals/id2variant/id2variant.c ../internals/utils/io_management.c  ../internals/seq_utils/isDNAbase.c ../internals/seq_utils/isIUPACbase.c  ../internals/seq_utils/is_dgnrt_mtch.c ../internals/seq_utils/test_possible_pairs.c ../internals/seq_utils/basemap.c ../internals/TECdisplay_navigator/parse_reference.c ../internals/TECdisplay_navigator/read_vbase.c -o id2variant
#
#
#
#
#build TECDisplay_Hnav
gcc ../internals/TECdisplay_Hnav/TECdisplay_Hnav.c ../internals/utils/io_management.c ../internals/utils/debug.c ../internals/seq_utils/basemap.c ../internals/seq_utils/isDNAbase.c ../internals/seq_utils/isIUPACbase.c ../internals/seq_utils/is_dgnrt_mtch.c  ../internals/seq_utils/test_possible_pairs.c ../internals/TECdisplay_mapper/TECdisplay_output_column_headers.c ../internals/TECdisplay_navigator/merge_values_files.c ../internals/TECdisplay_navigator/parse_TECdisplay_out_line.c ../internals/TECdisplay_navigator/parse_reference.c ../internals/TECdisplay_navigator/parse_constraints.c  ../internals/TECdisplay_navigator/read_vbase.c ../internals/TECdisplay_navigator/search_4_vbase_match.c ../internals/TECdisplay_Hnav/TECdisplay_Hnav_global_vars.c ../internals/TECdisplay_Hnav/get_constraint_metadata.c ../internals/TECdisplay_Hnav/calc_output_files.c ../internals/TECdisplay_Hnav/Hfilter.c ../internals/TECdisplay_Hnav/process_output_files.c  -lm -o TECdisplay_Hnav

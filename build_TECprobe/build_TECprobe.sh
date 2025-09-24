#!/bin/sh
#
#  build_TECprobe.sh
#  
#
#  Created by Eric Strobel on 2/7/23.
#  
#
#build cotrans_preprocessor
gcc -o cotrans_preprocessor ../internals/cotrans_preprocessor/cotrans_preprocessor_main.c ../internals/utils/debug.c ../internals/utils/gen_utils.c ../internals/utils/io_management.c ../internals/seq_utils/parse_fasta.c ../internals/seq_utils/mk_fasta.c ../internals/seq_utils/seq2bin_hash.c ../internals/seq_utils/seq2bin_long.c ../internals/seq_utils/isDNAbase.c ../internals/seq_utils/isIUPACbase.c ../internals/seq_utils/revcomp.c ../internals/seq_utils/appnd_att.c ../internals/variant_maker/constant_seqs.c ../internals/variant_maker/read_bcFile.c ../internals/cotrans_preprocessor/MLT_trgt_gen/mk_3pEnd_testdata.c ../internals/cotrans_preprocessor/MLT_trgt_gen/mk_MLT_trgts.c ../internals/cotrans_preprocessor/MLT_trgt_gen/mk_3pEnd_trgts.c ../internals/cotrans_preprocessor/MLT_trgt_gen/mk_intermed_trgts.c ../internals/cotrans_preprocessor/MLT_trgt_gen/printDels.c ../internals/cotrans_preprocessor/MLT_trgt_gen/printInsrts.c ../internals/cotrans_preprocessor/MLT_trgt_gen/printSubs.c ../internals/cotrans_preprocessor/MLT_trgt_gen/printQC_3pEnd_trgts.c ../internals/cotrans_preprocessor/SGL_trgt_gen/mk_SGL_target.c ../internals/cotrans_preprocessor/MUX_trgt_gen/mk_MUX_trgts.c ../internals/cotrans_preprocessor/MUX_trgt_gen/mk_MUX_testdata.c ../internals/cotrans_preprocessor/prcs_rds/UNV/call_fastp.c ../internals/cotrans_preprocessor/prcs_rds/UNV/bypass_fastp.c ../internals/cotrans_preprocessor/prcs_rds/UNV/prcs_chnl.c ../internals/cotrans_preprocessor/prcs_rds/UNV/print_splitting_metrics.c ../internals/cotrans_preprocessor/prcs_rds/UNV/mk_config.c ../internals/cotrans_preprocessor/prcs_rds/MLT/mk_smooth_script.c ../internals/cotrans_preprocessor/prcs_rds/MLT/parse_3pEnd_trgts.c  ../internals/cotrans_preprocessor/prcs_rds/MLT/prcs_MLT_cotrans.c ../internals/cotrans_preprocessor/prcs_rds/MLT/testdata3pEnd_analysis.c ../internals/cotrans_preprocessor/prcs_rds/MUX/prcs_MUX_cotrans.c ../internals/cotrans_preprocessor/prcs_rds/MUX/testdataMUX_analysis.c ../internals/cotrans_preprocessor/run_script_gen/UNV/config_struct.c ../internals/cotrans_preprocessor/run_script_gen/UNV/mk_run_script.c  ../internals/cotrans_preprocessor/run_script_gen/UNV/check_config.c ../internals/cotrans_preprocessor/run_script_gen/UNV/parse_config.c ../internals/cotrans_preprocessor/run_script_gen/UNV/mk_run_nm.c  ../internals/cotrans_preprocessor/run_script_gen/MLT/print_MLT_SM2_script.c -lm
#
#
#
#
#build mkmtrx
gcc -o mkmtrx ../internals/mkmtrx/mkmtrx.c ../internals/utils/debug.c ../internals/utils/gen_utils.c ../internals/utils/io_management.c ../internals/cotrans_preprocessor/run_script_gen/UNV/config_struct.c ../internals/cotrans_preprocessor/run_script_gen/UNV/mk_run_nm.c  ../internals/process_TECprobe_profiles/UNV/parse_sample_name.c ../internals/mkmtrx/cotrans_mtrx.c ../internals/mkmtrx/find_nxt_dir_entry.c ../internals/mkmtrx/get_fastp_out.c ../internals/mkmtrx/get_profiles.c ../internals/mkmtrx/parse_log.c ../internals/mkmtrx/mk_output_files.c  ../internals/mkmtrx/mk_rdat.c
#
#
#
#
#build mtrx2cols
gcc -o mtrx2cols ../internals/mtrx2cols/mtrx2cols.c ../internals/utils/io_management.c ../internals/mkmtrx/cotrans_mtrx.c
#
#
#
#
#build process_TECprobeVL_profiles
gcc -o process_TECprobeVL_profiles ../internals/process_TECprobe_profiles/process_TECprobe_profiles.c ../internals/utils/debug.c ../internals/utils/io_management.c ../internals/utils/gen_utils.c ../internals/seq_utils/isRNAbase.c ../internals/process_TECprobe_profiles/UNV/store_SM2_profile.c ../internals/process_TECprobe_profiles/UNV/initialize_empty_profile.c ../internals/process_TECprobe_profiles/UNV/calculate_normalization_factor.c ../internals/process_TECprobe_profiles/UNV/read_analysis_directories.c ../internals/process_TECprobe_profiles/UNV/input_validation.c ../internals/process_TECprobe_profiles/UNV/make_output_directories.c  ../internals/process_TECprobe_profiles/UNV/normalize_reactivities.c ../internals/process_TECprobe_profiles/UNV/merge_profiles.c ../internals/process_TECprobe_profiles/UNV/parse_sample_name.c ../internals/process_TECprobe_profiles/UNV/generate_sample_name.c ../internals/process_TECprobe_profiles/UNV/print_merged_profiles.c ../internals/process_TECprobe_profiles/UNV/print_legacy_compiled_table.c  ../internals/cotrans_preprocessor/run_script_gen/UNV/config_struct.c ../internals/cotrans_preprocessor/run_script_gen/UNV/mk_run_nm.c -lm
#
#
#
#
#build assemble_TECprobeLM_data
gcc -o assemble_TECprobeLM_data ../internals/assemble_TECprobeLM_data/assemble_TECprobeLM_data.c  ../internals/utils/debug.c ../internals/utils/io_management.c  ../internals/assemble_TECprobeLM_data/store_ipt_name.c ../internals/assemble_TECprobeLM_data/validate_input.c ../internals/assemble_TECprobeLM_data/count_delims_2_col.c ../internals/assemble_TECprobeLM_data/get_value.c ../internals/assemble_TECprobeLM_data/print_input_filenames.c ../internals/assemble_TECprobeLM_data/print_reactivity_output.c ../internals/assemble_TECprobeLM_data/print_linebar_output.c ../internals/assemble_TECprobeLM_data/print_length_dist_output.c ../internals/cotrans_preprocessor/run_script_gen/UNV/config_struct.c ../internals/cotrans_preprocessor/run_script_gen/UNV/mk_run_nm.c ../internals/process_TECprobe_profiles/UNV/parse_sample_name.c  ../internals/mkmtrx/mk_rdat.c
#
#
#
#
#build draw_intermediates
gcc -o draw_intermediates ../internals/draw_intermediates/draw_intermediates.c ../internals/utils/io_management.c ../internals/mkmtrx/cotrans_mtrx.c ../internals/seq_utils/parse_fasta.c ../internals/seq_utils/isDNAbase.c
#
#
#
#

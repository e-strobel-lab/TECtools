#!/bin/sh
#
#  build_TECprobe.sh
#  
#
#  Created by Eric Strobel on 2/7/23.
#  
#
#build cotrans_preprocessor
gcc -Wall -o cotrans_preprocessor ../internals/cotrans_preprocessor/cotrans_preprocessor_main.c ../internals/utils/debug.c ../internals/utils/gen_utils.c ../internals/utils/io_management.c ../internals/seq_utils/parse_fasta.c ../internals/seq_utils/mk_fasta.c ../internals/seq_utils/seq2bin_hash.c ../internals/seq_utils/isDNAbase.c ../internals/seq_utils/isIUPACbase.c ../internals/seq_utils/revcomp.c ../internals/seq_utils/appnd_att.c ../internals/cotrans_preprocessor/MLT_trgt_gen/mk_3pEnd_testdata.c ../internals/cotrans_preprocessor/MLT_trgt_gen/mk_MLT_trgts.c ../internals/cotrans_preprocessor/MLT_trgt_gen/mk_3pEnd_trgts.c ../internals/cotrans_preprocessor/MLT_trgt_gen/mk_intermed_trgts.c ../internals/cotrans_preprocessor/MLT_trgt_gen/printDels.c ../internals/cotrans_preprocessor/MLT_trgt_gen/printInsrts.c ../internals/cotrans_preprocessor/MLT_trgt_gen/printSubs.c ../internals/cotrans_preprocessor/MLT_trgt_gen/printQC_3pEnd_trgts.c ../internals/cotrans_preprocessor/SGL_trgt_gen/mk_SGL_target.c ../internals/cotrans_preprocessor/prcs_rds/UNV/call_fastp.c ../internals/cotrans_preprocessor/prcs_rds/UNV/prcs_chnl.c ../internals/cotrans_preprocessor/prcs_rds/MLT/mk_MLT_config.c ../internals/cotrans_preprocessor/prcs_rds/MLT/mk_smooth_script.c ../internals/cotrans_preprocessor/prcs_rds/MLT/parse_3pEnd_trgts.c ../internals/cotrans_preprocessor/prcs_rds/MLT/printQC_prcsMLT.c ../internals/cotrans_preprocessor/prcs_rds/MLT/prcs_MLT_cotrans.c ../internals/cotrans_preprocessor/prcs_rds/MLT/testdata3pEnd_analysis.c ../internals/cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.c ../internals/cotrans_preprocessor/run_script_gen/MLT/check_MLT_config.c ../internals/cotrans_preprocessor/run_script_gen/MLT/mk_MLT_run_nm.c ../internals/cotrans_preprocessor/run_script_gen/MLT/mk_MLT_run_script.c ../internals/cotrans_preprocessor/run_script_gen/MLT/parse_MLT_config.c ../internals/cotrans_preprocessor/run_script_gen/MLT/print_MLT_SM2_script.c -lm
#
#
#
#
#build mkmtrx
gcc -Wall -o mkmtrx ../internals/mkmtrx/mkmtrx.c ../internals/utils/debug.c ../internals/utils/gen_utils.c ../internals/utils/io_management.c ../internals/cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.c ../internals/cotrans_preprocessor/run_script_gen/MLT/mk_MLT_run_nm.c  ../internals/process_TECprobe_profiles/VL/parse_VL_sample_name.c ../internals/mkmtrx/cotrans_mtrx.c ../internals/mkmtrx/find_nxt_dir_entry.c ../internals/mkmtrx/get_fastp_out.c ../internals/mkmtrx/get_profiles.c ../internals/mkmtrx/parse_log.c ../internals/mkmtrx/mk_output_files.c  ../internals/mkmtrx/mk_rdat.c
#
#
#
#
#build mtrx2cols
gcc -Wall -o mtrx2cols ../internals/mtrx2cols/mtrx2cols.c ../internals/utils/io_management.c ../internals/mkmtrx/cotrans_mtrx.c
#
#
#
#
#build process_TECprobeVL_profiles
gcc -Wall -o process_TECprobeVL_profiles ../internals/process_TECprobe_profiles/VL/process_TECprobeVL_profiles.c ../internals/utils/debug.c ../internals/utils/io_management.c ../internals/utils/gen_utils.c ../internals/seq_utils/isRNAbase.c ../internals/process_TECprobe_profiles/global/store_SM2_profile.c ../internals/process_TECprobe_profiles/global/initialize_empty_profile.c ../internals/process_TECprobe_profiles/global/calculate_normalization_factor.c ../internals/process_TECprobe_profiles/VL/read_VL_analysis_directories.c ../internals/process_TECprobe_profiles/VL/VL_input_validation.c ../internals/process_TECprobe_profiles/VL/make_VL_output_directories.c  ../internals/process_TECprobe_profiles/VL/normalize_VL_reactivities.c ../internals/process_TECprobe_profiles/VL/merge_VL_profiles.c ../internals/process_TECprobe_profiles/VL/parse_VL_sample_name.c ../internals/process_TECprobe_profiles/VL/generate_VL_sample_name.c ../internals/process_TECprobe_profiles/VL/print_merged_VL_profiles.c ../internals/process_TECprobe_profiles/VL/print_legacy_compiled_table.c  ../internals/cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.c ../internals/cotrans_preprocessor/run_script_gen/MLT/mk_MLT_run_nm.c -lm
#
#
#
#
#build assemble_TECprobeLM_data
gcc -Wall -o assemble_TECprobeLM_data ../internals/assemble_TECprobeLM_data/assemble_TECprobeLM_data.c  ../internals/utils/debug.c ../internals/utils/io_management.c  ../internals/assemble_TECprobeLM_data/store_ipt_name.c ../internals/assemble_TECprobeLM_data/validate_input.c ../internals/assemble_TECprobeLM_data/count_delims_2_col.c ../internals/assemble_TECprobeLM_data/get_value.c ../internals/assemble_TECprobeLM_data/print_input_filenames.c ../internals/assemble_TECprobeLM_data/print_reactivity_output.c ../internals/assemble_TECprobeLM_data/print_linebar_output.c ../internals/assemble_TECprobeLM_data/print_length_dist_output.c ../internals/cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.c ../internals/process_TECprobe_profiles/VL/parse_VL_sample_name.c ../internals/cotrans_preprocessor/run_script_gen/MLT/mk_MLT_run_nm.c ../internals/mkmtrx/mk_rdat.c
#
#
#
#
#build draw_intermediates
gcc -Wall -o draw_intermediates ../internals/draw_intermediates/draw_intermediates.c ../internals/utils/io_management.c ../internals/mkmtrx/cotrans_mtrx.c ../internals/seq_utils/parse_fasta.c ../internals/seq_utils/isDNAbase.c
#
#
#
#

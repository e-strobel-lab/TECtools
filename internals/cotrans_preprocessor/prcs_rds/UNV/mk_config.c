//
//  mk_config.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

#include "../../run_script_gen/MLT/mk_MLT_run_script.h"

#include "mk_config.h"

/* mk_config:
 generate configuration file for use make_shapemapper2_run_script executable
 this config file is used to generate a shell script that run commands for
 shapemapper2 analysis
*/

void mk_config(TPROBE_names * nm, target3p_params * trg_prms, int mode)
{
    extern const char VL_LM_config_header[36];
    extern const char SL_config_header[24];
    extern const char MUX_config_header[25];
    
    FILE * out_fp = NULL;
    
    if ((out_fp = fopen("./config.txt", "w")) == NULL) {
        printf("make_config: ERROR - could not generate config file. Aborting program...\n");
        abort();
    }
    
    //print header line indicating the format of the config file
    if (mode == MULTI) {
        fprintf(out_fp, "%s\n", VL_LM_config_header);
    } else if (mode == SINGLE) {
        fprintf(out_fp, "%s\n", SL_config_header);
    } else if (mode == MULTIPLEX) {
        fprintf(out_fp, "%s\n", MUX_config_header);
    }
    
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# ***** shapemapper2 path *****\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# set shapemapper2_path to the filepath\n");
    fprintf(out_fp, "# of the shapemapper2 executable\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "> shapemapper2_path=FILL_IN\n");
    fprintf(out_fp, "#\n# ---------------------------------------------------------\n");
    fprintf(out_fp, "# ***** input file location *****\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# set input_file_location to the filepath of the folder that\n");
    fprintf(out_fp, "# contains the split fastq files to be analyzed\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "> input_file_location=FILL_IN\n");
    fprintf(out_fp, "#\n# ---------------------------------------------------------\n");
    fprintf(out_fp, "# ***** output file location *****\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# set output_file_location to the filepath of the folder that\n");
    fprintf(out_fp, "# that will be used to store the shapemapper2 output files\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "> output_file_location=FILL_IN\n");
    fprintf(out_fp, "#\n# ---------------------------------------------------------\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# ***** input file prefixes (auto-populated) *****\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "> untreated_read1_prefix=%s_UNT\n", nm->smpl[READ1]);
    fprintf(out_fp, "> untreated_read2_prefix=%s_UNT\n", nm->smpl[READ2]);
    fprintf(out_fp, ">  modified_read1_prefix=%s_MOD\n", nm->smpl[READ1]);
    fprintf(out_fp, ">  modified_read2_prefix=%s_MOD\n", nm->smpl[READ2]);
    fprintf(out_fp, "#\n# ---------------------------------------------------------\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# ***** target files location *****\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# set target_files_location to the full path of the folder\n");
    fprintf(out_fp, "# that contains the target files to be used for alignment\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "> target_files_location=FILL_IN\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "#\n# ---------------------------------------------------------\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# ***** target files prefix *****\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# set target_files_prefix to the name that was used for target\n");
    fprintf(out_fp, "# generation. e.g., if the the target name is ZTP_NNN.fa, where\n");
    fprintf(out_fp, "# 'NNN' is the transcript length, target_files_prefix=ZTP\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "> target_files_prefix=FILL_IN\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "#\n# ---------------------------------------------------------\n");
    
    //only print target min and max length fields for VL/LM/SL experiments
    if (mode == MULTI || mode == SINGLE) {
        fprintf(out_fp, "#\n");
        fprintf(out_fp, "# ***** transcript length bounds (auto-populated) *****\n");
        fprintf(out_fp, "#\n");
        fprintf(out_fp, "> min_target_length=%03d\n", trg_prms->min);
        fprintf(out_fp, "> max_target_length=%03d\n", trg_prms->max);
        fprintf(out_fp, "#\n# ---------------------------------------------------------\n");
    }
    
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# ***** sample name *****\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# provide a name for the RNA. this name will be used when generating \n");
    fprintf(out_fp, "# shapemapper2 output files and should only include alphanumeric \n");
    fprintf(out_fp, "# characters or underscores\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "> name=FILL_IN\n");
    fprintf(out_fp, "#\n# ---------------------------------------------------------\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# ***** chemical probe conditions *****\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# specify the chemical probe that was used for use in the shapemapper2\n");
    fprintf(out_fp, "# name variable and whether probing was cotranscriptional or equilibrium\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# set cotranscriptional to TRUE for cotranscriptional or FALSE for equilibrium\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# suggested codes for common probes:\n");
    fprintf(out_fp, "#   benzoyl cyanide:  BzCN\n");
    fprintf(out_fp, "#   dimethyl sulfate: DMS\n");
    fprintf(out_fp, "#   1-methyl-7-nitroisatoic anhydride: 1M7\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "> cotranscriptional=FILL_IN\n");
    fprintf(out_fp, "> chemical_probe=FILL_IN\n");
    fprintf(out_fp, "#\n# ---------------------------------------------------------\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# ***** run identifiers *****\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# specify whether the sequencing data were concatenated\n");
    fprintf(out_fp, "# from multiple runs and add a sequencing run ID(s) to\n");
    fprintf(out_fp, "# the shapemapper2 name variable\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# set concatenated as TRUE or FALSE (e.g. concatenated=FALSE)\n");
    fprintf(out_fp, "#   if concatenated=FALSE, provide no more than one runID\n");
    fprintf(out_fp, "#   if concatenated=TRUE, at least two runIDs must be provided\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# the maximum number of additional fields is 8\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "> concatenated=FILL_IN\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "> runID=FILL_IN\n");
    fprintf(out_fp, "> runID=NULL\n");
    fprintf(out_fp, "> runID=NULL\n");
    fprintf(out_fp, "> runID=NULL\n");
    fprintf(out_fp, "> runID=NULL\n");
    fprintf(out_fp, "> runID=NULL\n");
    fprintf(out_fp, "> runID=NULL\n");
    fprintf(out_fp, "> runID=NULL\n");
    fprintf(out_fp, "#\n# ---------------------------------------------------------\n");
    
    //only print smoothing field for VL/LM experiments
    if (mode == MULTI) {
        fprintf(out_fp, "#\n");
        fprintf(out_fp, "# ***** smoothing *****\n");
        fprintf(out_fp, "#\n");
        fprintf(out_fp, "# will the analyses be performed on smoothed data?\n");
        fprintf(out_fp, "# set as TRUE or FALSE. (e.g. smoothing=FALSE)\n");
        fprintf(out_fp, "#\n");
        fprintf(out_fp, "> smoothing=FILL_IN\n");
        fprintf(out_fp, "#\n# ---------------------------------------------------------\n");
    }

    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# ***** ligand conditions (optional) *****\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# specify the ligand that was used and its concentration\n");
    fprintf(out_fp, "# for use in the shapemapper2 name variable\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# ligand_name is user-defined\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# the format for ligand_conc is NNNxM, where\n");
    fprintf(out_fp, "#   NNN is a three digit number, and\n");
    fprintf(out_fp, "#   x is the unit (e.g., m, u, or n for milli, micro, and nano\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "> ligand_name=NULL\n");
    fprintf(out_fp, "> ligand_conc=NULL\n");
    fprintf(out_fp, "#\n# ---------------------------------------------------------\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# ***** additional fields (optional) *****\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "# optional fields for user-defined information\n");
    fprintf(out_fp, "# for use in the shapemapper2 name variable\n");
    fprintf(out_fp, "# the maximum number of additional fields is 8\n");
    fprintf(out_fp, "#\n");
    fprintf(out_fp, "> field=NULL\n");
    fprintf(out_fp, "> field=NULL\n");
    fprintf(out_fp, "> field=NULL\n");
    fprintf(out_fp, "> field=NULL\n");
    fprintf(out_fp, "> field=NULL\n");
    fprintf(out_fp, "> field=NULL\n");
    fprintf(out_fp, "> field=NULL\n");
    fprintf(out_fp, "> field=NULL\n");
    fprintf(out_fp, "#\n# ---------------------------------------------------------\n#\n");
    
    if ((fclose(out_fp)) == EOF) {
        printf("mk_config: ERROR - error occurred when closing config file. Aborting program...\n");
        abort();
    }
}

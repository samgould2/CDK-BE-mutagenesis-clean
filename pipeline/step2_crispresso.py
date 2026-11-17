"""Step 2 of the sensor pipeline: run CRISPResso per guide.

Runs on the cluster, once per sample, after step 1 has split the reads by
guide. Produces the per-guide CRISPResso output that step 3 reduces.

    python3 step2_crispresso.py <guide_table.csv> <sample_name> [outdir]

The quantification window is the 20 nt protospacer plus 5 nt either side
(coordinates 12-32). A dummy HDR amplicon carrying the intended edit is passed
so CRISPResso reports the edited allele separately.

This file replaces six identical copies in the original cluster_scripts/ tree.
"""
import pandas as pd
import os
import sys
from pathlib import Path


#-------Loading in system args
# Parse command-line arguments to read in files of interest
if len(sys.argv) not in (3, 4):
    print("Usage: python3 step2_crispresso.py <input_df> <sample_name> [outdir]")
    sys.exit(1)

df = pd.read_csv(Path(sys.argv[1]))
sample_name = str(sys.argv[2])
# The six per-screen copies of this script were identical except for this
# directory ('crispresso' vs 'crispresso3'); it is an argument now.
outdir = str(sys.argv[3]) if len(sys.argv) == 4 else "./crispresso"

#------run it for all guides
for i, val in df.iterrows():
    guide_id = val['gRNA_id']
    wt = val['sensor_wt']
    edited = val['sensor_alt']
    proto = val['Protospacer'][1:] #19 bp (excluding G start)
    #naming according to unique id since i'm performing multiple crispresso analyses per guide
    name = guide_id

    #generating crispresso command
    f_path = f"./{sample_name}/{guide_id}.fastq"
    start = f"CRISPResso --fastq_r1 {f_path} "

    #qwc is 5nt upstream and 5nt downstream of 20 nt guide sequence
    #default alignment parameters (60%)
    p2 = f'--amplicon_seq {wt} --expected_hdr_amplicon_seq {edited} --guide_seq {proto} --quantification_window_coordinates 12-32,12-32 --exclude_bp_from_left 0 --exclude_bp_from_right 0 --base_editor_output --plot_window_size 5 --suppress_report --suppress_plots '
    
    out_path = f"-o {outdir}/{sample_name} -n {name}"

    command = start+ p2 + out_path

    #running the crispresso command
    #apparently this is bad practice (but subprocess module didn't work...)
    os.system(command)

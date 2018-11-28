#!/bin/sh -f

# nicrabgrp/cli/process.py --help
# shows help message of sub-commands

#nicrabgrp/cli/process.py show-parameters setenv/parameter/setup_parameter.yaml 
# shows parameters

#rm -rf out/crabgrp_v181127
#nicrabgrp/cli/process.py convert-radiofiles setenv/parameter/setup_parameter.yaml --outdir out/crabgrp_v181127
# convert radio text files to fits files 

#rm -rf out/crabgrp_v181127/2018097/xray
#nicrabgrp/cli/process.py prepare-xrayfiles setenv/parameter/setup_parameter.yaml --outdir out/crabgrp_v181127
# make barycentric corrected event files with a PULSE_PHASE column

#nicrabgrp/cli/process.py add-grpflag-to-xrayfiles setenv/parameter/setup_parameter.yaml --outdir out/crabgrp_v181127
# add grp flags to the file

#rm -rf out/crabgrp_v181128/combined

#nicrabgrp/cli/process.py plot-individual-profiles setenv/parameter/setup_parameter.yaml --outdir out/crabgrp_v181128
# plot individual profiles 

rm -rf out/crabgrp_v181128_accumulation
nicrabgrp/cli/process.py accumulated-significance setenv/parameter/setup_parameter.yaml --indir out/crabgrp_v181128 --outdir out/crabgrp_v181128_accumulation
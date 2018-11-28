#!/bin/sh -f 

#nicrabgrp/cli/xray.py generate-profile-fitsfile \
#	out/crabgrp_v181128/2017314/xray/ni1013010110_0mpu7_cl_nibary_phase_grpflag_0p3to12p0keV.evt \
#	--outfitsfile "out1.fits" 

#nicrabgrp/cli/xray.py generate-profile-fitsfile \
#	out/crabgrp_v181128/2017314/xray/ni1013010110_0mpu7_cl_nibary_phase_grpflag_0p3to12p0keV.evt \
#	out/crabgrp_v181128/2017364/xray/ni1013010122_0mpu7_cl_nibary_phase_grpflag_0p3to12p0keV.evt \
#	--outfitsfile "out2.fits" 

#nicrabgrp/cli/xray.py plot-profile-fitsfile	out1.fits --outpdf out1.pdf

nicrabgrp/cli/xray.py get-enhancement-significance out/crabgrp_v181128/2017314/xray/ni1013010110_0mpu7_cl_nibary_phase_grpflag_0p3to12p0keV_pls.fits 246,247,248,249 --nphase 250
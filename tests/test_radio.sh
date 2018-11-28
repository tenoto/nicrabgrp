#!/bin/sh -f

nicrabgrp/cli/radio.py --help 

nicrabgrp/cli/radio.py convert-grp-txt2fits \
	/Users/enoto/Dropbox/enoto/library/nicrabgrp/data/radio/original/v181127/allradioSNge4.8PHpm3/2017221_UsdS_IPGRPlistDE430SNge4.8PHpm3_updated.txt \
	--outfitsfile 2017221_UsdS_IPGRPlistDE430SNge4.8PHpm3_updated.fits

nicrabgrp/cli/radio.py convert-gti-txt2fits \
	/Users/enoto/Dropbox/enoto/library/nicrabgrp/data/radio/original/v181127/allradioSNge4.8PHpm3/2017221_UsdS_DE430NEW2_GTI.txt \
	--outfitsfile 2017221_UsdS_DE430NEW2_GTI.fits	
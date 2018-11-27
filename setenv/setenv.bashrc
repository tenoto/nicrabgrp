#!/bin/bash 

echo "#############################"
echo "#     Giant Radio Pulse     #"
echo "#############################"

export PATH_NICRABGRP=$(pwd)
export PATH=$PATH_NICRABGRP/nicrabgrp/cli:$PATH
export PYTHONPATH=$PATH_NICRABGRP:$PYTHONPATH

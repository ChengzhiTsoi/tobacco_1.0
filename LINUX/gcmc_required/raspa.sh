#!/bin/bash
if [ -r /opt/crc/Modules/current/init/bash ]; then
        source /opt/crc/Modules/current/init/bash
fi
# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi


cd /scratch365/ccai/ADTL/LINUX/INDEX

export HOME="/scratch365/ccai"
export RASPA_DIR=${HOME}/Research/RASPA/simulations
export DYLD_LIBRARY_PATH=${RASPA_DIR}/lib
export LD_LIBRARY_PATH=${RASPA_DIR}/lib:$LD_LIBRARY_PATH
$RASPA_DIR/bin/simulate -i simulation.input
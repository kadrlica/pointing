#!/usr/bin/env bash
#$Id: setup.sh 11004 2015-08-29 20:11:53Z kadrlica $

#These should be set in DECamObserver rc file
#setup SISPIlib
#setup ephem
#setup matplotlib

POINTING_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PATH=${POINTING_DIR}/bin:${PATH}
#export PYTHONPATH=${POINTING_DIR}/python/PlotPointings:${PYTHONPATH}
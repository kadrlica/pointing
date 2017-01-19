#!/usr/bin/env bash
#$Id: setup.sh 11004 2015-08-29 20:11:53Z kadrlica $

#These should be set in DECamObserver rc file
#setup SISPIlib
#setup ephem
#setup matplotlib

export PLOTPOINTINGS_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PATH=${PLOTPOINTINGS_DIR}/bin:${PATH}
export PYTHONPATH=${PLOTPOINTINGS_DIR}/python/PlotPointings:${PYTHONPATH}

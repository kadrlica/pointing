#$Id: setup.csh 11004 2015-08-29 20:11:53Z kadrlica $

# These should all be set in DECamObserver rc file
#setup SISPIlib
#setup ephem
#setup matplotlib

set sourced=($_)
if ("$sourced" != "") then
  set basedir = `dirname $sourced[2]`
else
  set basedir = `dirname $0`
endif

setenv PLOTPOINTINGS_DIR `cd $basedir/.. && pwd`
setenv PATH ${PLOTPOINTINGS_DIR}/bin:${PATH}
setenv PYTHONPATH ${PLOTPOINTINGS_DIR}/python/PlotPointings:${PYTHONPATH}

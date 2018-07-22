#!/usr/bin/env python
""" Where you at? """
import sys,os
import logging
from collections import OrderedDict as odict
from datetime import datetime,timedelta,tzinfo
import dateutil.parser

import mpl_toolkits.basemap as basemap
from matplotlib.patches import Ellipse, Circle
import matplotlib.patheffects as patheffects
from _tkinter import TclError

import numpy as np
import pylab as plt
import ephem

__author__  = "Alex Drlica-Wagner"
__email__   = "kadrlica@fnal.gov"
__version__ = "2.1.2"

MAXREF=5000 # Maximum number of refreshes
DECAM=1.1   # DECam radius (deg)

# Accurate DECam marker size depends on figsize and DPI
# This is a mess...
FIGSIZE=(10.5,8.5)
SCALE=np.sqrt((8.0*6.0)/(FIGSIZE[0]*FIGSIZE[1]))
DPI=80;

FILTERS = ['u','g','r','i','z','Y','VR']
BANDS = FILTERS + ['all']
COLORS = odict([
    ('none','black'),
    ('u','blue'),
    ('g','green'),
    ('r','red'),
    ('i','gold'),
    ('z','magenta'),
    ('Y','black'),
    ('VR','gray'),
])

# Allowed map projections
PROJ = odict([
    ('ortho'  , dict(projection='ortho',celestial=True)),
    ('moll'   , dict(projection='moll',celestial=True)),
    ('mol'    , dict(projection='moll',celestial=True)),
    ('ait'    , dict(projection='hammer',celestial=True)),
    ('mbt'    , dict(projection='mbtfpq',celestial=True)),
    ('mbtfpq' , dict(projection='mbtfpq',celestial=True)),
    ('mcbryde', dict(projection='mbtfpq',celestial=True)),
])

# Derived from telra,teldec of 10000 exposures
SN = odict([
    ('E1',(7.874, -43.010)),
    ('E2',(9.500, -43.999)),
    ('X1',(34.476, -4.931)),
    ('X2',(35.664,-6.413)),
    ('X3',(36.449, -4.601)),
    ('S1',(42.818, 0.000)),
    ('S2',(41.193, -0.991)),
    ('C1',(54.274, -27.113)),
    ('C2',(54.274, -29.090)),
    ('C3',(52.647, -28.101)),
])

SN_LABELS = odict([
    ('SN-E',(8,-41)),
    ('SN-X',(35,-12)),
    ('SN-S',(45,1)),
    ('SN-C',(55,-35)),
])

# The allowed footprint outlines
FOOTPRINTS = ['none','des','des-sn','smash','maglites','bliss','decals']

# CTIO location taken from:
#http://www.ctio.noao.edu/noao/content/Coordinates-Observatories-Cerro-Tololo-and-Cerro-Pachon
#http://arxiv.org/pdf/1210.1616v3.pdf
#(-30h 10m 10.73s, -70h 48m 23.52s, 2213m)

TEL_LON = -70.80653
TEL_LAT = -30.169647
TEL_HEIGHT = 2213

# Create the observatory object
CTIO = ephem.Observer()
CTIO.lon,CTIO.lat = str(TEL_LON),str(TEL_LAT)
CTIO.elevation = TEL_HEIGHT

def get_datadir():
    """ Path to data directory. """
    return os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

# Stupid timezone definition
ZERO = timedelta(0)
HOUR = timedelta(hours=1)
class UTC(tzinfo):
    """UTC"""
    def utcoffset(self, dt):
        return ZERO

    def tzname(self, dt):
        return "UTC"

    def dst(self, dt):
        return ZERO

def safe_proj(bmap,lon,lat,inverse=False):
    """ Remove points outside of projection 
    
    Parameters:
    -----------
    bmap : basemap
    lon  : longitude
    lat  : latitude
    inverse : inverse projection

    Returns:
    --------
    x,y : projected coordinates    
    """
    x,y = bmap(np.atleast_1d(lon),np.atleast_1d(lat),inverse=inverse)
    x[np.abs(x) > 1e29] = None
    y[np.abs(y) > 1e29] = None
    return x,y


def airmass_angle(x=1.4):
    """ Zenith angle for a given airmass limit """
    return 90.-np.degrees(np.arcsin(1./x))

def load_data(opts):
    """ Load the data (either from DB of file). 

    Parameters:
    -----------
    opts : command line options
    
    Returns:
    --------
    data : numpy recarray
    """
    since = parse_since(opts.since)
    propid = '%' if opts.propid is None else opts.propid
    dtype=[('expnum',int),('telra',float),('teldec',float),('filter',object)]

    if opts.infile is None:
        selection = ['id','telra','teldec','filter']
        #filter = "exposed = TRUE AND flavor LIKE '%s' AND date > '%s' AND propid LIKE '%s' ORDER BY id DESC"%(opts.flavor,since.isoformat(),propid)
        filter = "exposed = TRUE AND flavor SIMILAR TO '%s' AND date > '%s' AND propid LIKE '%s' ORDER BY id DESC"%(opts.flavor,since.isoformat(),propid)
        # Use the FNAL mirror to avoid overloading CTIO
        try: from database import Database
        except ImportError: from pointing.database import Database
        db = Database(dbname='db-'+opts.db)
        db.connect()
        query = "SELECT %s FROM exposure WHERE %s"%(','.join(selection),filter)
        #query = "SELECT id as expnum,telra as ra,teldec as dec,filter as band FROM exposure WHERE exposed = TRUE AND flavor LIKE 'object' and telra between 80 and 82 AND teldec between -71 and -69"
        data = db.execute(query)

        if len(data): ret = np.rec.array(data,dtype=dtype)
        else:         ret = np.rec.recarray(0,dtype=dtype)

        return ret
    else:
        return np.loadtxt(opts.infile,dtype=dtype)

def mjd(datetime):
    """ Modified Julian Date (MJD) """
    mjd_epoch = dateutil.parser.parse('1858-11-17T00:00:00Z')
    mjd_date = (datetime-mjd_epoch).total_seconds()/float(24*60*60)
    return mjd_date

def lmst(observatory):
    """ Calculate Local Mean Sidereal Time (LMST) """
    lmst = np.degrees(observatory.sidereal_time())
    logging.debug('Using pyephem for LMST: %.3f'%lmst)
    return lmst

def moon(datetime):
    """ Moon location 
    
    Parameters:
    -----------
    datetime : the datetime of moon location request
    
    Returns:
    --------
    (ra, dec), phase : moon parameters [(deg, deg), %]
    """
    moon = ephem.Moon()
    moon.compute(CTIO)
    moon_phase = moon.moon_phase * 100
    moon_ra,moon_dec = np.degrees([moon.ra,moon.dec])
    return (moon_ra, moon_dec),moon_phase

def boolean(string):
    """ Convert strings to booleans for argparse """
    string = string.lower()
    if string in ['0', 'f', 'false', 'no', 'off']:
        return False
    elif string in ['1', 't', 'true', 'yes', 'on']:
        return True
    else:
        raise ValueError()

def splash_screen():
    """ Splash text to print """
    splash = """Running Alex Drlica-Wagner's DECam pointing script..."""
    logging.info(splash)

def parse_utc(value):
    """ Parse isoformat 'utc' option string. """
    if value is None:
        utc = datetime.now(tz=UTC())
    elif isinstance(value,datetime):
        utc = value
    else:
        utc = dateutil.parser.parse(value,tzinfos={'UTC':UTC})
    logging.debug("UTC: %s"%utc.strftime('%Y-%m-%d %H:%M:%S'))
    return utc

def parse_since(value):
    """ Parse isoformat 'since' option string. """
    if value is None:
        since = datetime.now(tz=UTC()) - timedelta(hours=12)
    elif isinstance(value,datetime):
        since = value
    elif value.lower() in ['all','none','forever']:
        since = dateutil.parser.parse('2012-01-01 12:00',tzinfos={'UTC':UTC})
    else:
        since = dateutil.parser.parse(value,tzinfos={'UTC':UTC})
    logging.debug("Since: %s"%since.strftime('%Y-%m-%d %H:%M:%S'))
    return since

def draw_constellation(bmap,name):
    """ Draw a map of the constellations (work in progress). """
    from constellations import CONSTELLATIONS
    points = np.array(CONSTELLATIONS[name])

    drawtype = points[:,0]
    radeg = points[:,1] * 1.0 / 1800 * 15
    decdeg = points[:,2] * 1.0 / 60
    print radeg,decdeg
    verts = zip(safe_proj(bmap,radeg,decdeg))
    codes = [XEPHEM2PATH[c] for c in points[:,0]]
    print x,y

def draw_des(bmap,**kwargs):
    """
    Plot the DES wide-field footprint.

    Parameters:
    -----------
    bmap   : The basemap object
    kwargs : Various plotting arguments

    Returns:
    --------
    None
    """
    # Plot the wide-field survey footprint
    logging.debug("Plotting footprint: %s"%opts.footprint)
    #basedir = os.path.dirname(os.path.abspath(__file__))
    infile = os.path.join(get_datadir(),'des-round17-poly.txt')
    perim = np.loadtxt(infile,dtype=[('ra',float),('dec',float)])
    proj = safe_proj(bmap,perim['ra'],perim['dec'])
    bmap.plot(*proj,**kwargs)

def draw_des_sn(bmap,**kwargs):
    """
    Plot the DES supernova fields.

    Parameters:
    -----------
    bmap   : The basemap object
    kwargs : Various plotting arguments

    Returns:
    --------
    None
    """
    # Plot the SN fields
    logging.debug("Plotting DES supernova fields.")
    # Check that point inside boundary
    # Doesn't work for 'ait' and 'moll' projections
    fact = 0.99
    projection = kwargs.pop('projection',None)
    if projection in basemap._pseudocyl:
        # This was estimated by eye...
        rminor=9.00995e6; rmajor = 2*rminor
        boundary = Ellipse((rmajor,rminor),
                           2*(fact*rmajor),2*(fact*rminor))
    else:
        boundary = Ellipse((bmap.rmajor,bmap.rminor),
                           2*(fact*bmap.rmajor),2*(fact*bmap.rminor))
    for v in SN.values():
        if not boundary.contains_point(bmap(*v)):
            continue
        # This does the projection correctly, but fails at boundary
        bmap.tissot(v[0],v[1],DECAM,100,**kwargs)

    # The SN labels
    sntxt_kwargs = dict(zorder=kwargs['zorder'],fontsize=12,
                        bbox=dict(boxstyle='round,pad=0',fc='w',ec='none',
                                  alpha=0.25))
    for k,v in SN_LABELS.items():
        plt.gca().annotate(k,bmap(*v),**sntxt_kwargs)

def draw_smash(bmap,**kwargs):
    """ Draw the SMASH fields 

    Parameters:
    -----------
    bmap   : The basemap object
    kwargs : Various plotting arguments

    Returns:
    --------
    None
    """
    filename = os.path.join(get_datadir(),'smash_fields_final.txt')

    smash=np.genfromtxt(filename,dtype=[('ra',float),('dec',float)],usecols=[4,5])
    smash_x,smash_y = safe_proj(bmap,smash['ra'],smash['dec'])
    kwargs.update(dict(facecolor='none'))
    bmap.scatter(smash_x,smash_y,color='k',**kwargs)

def draw_maglites(bmap,**kwargs):
    """
    Plot the MagLiteS Phase-I footprint.

    Parameters:
    -----------
    bmap   : The basemap object
    kwargs : Various plotting arguments

    Returns:
    --------
    None
    """

    # Plot the wide-field survey footprint
    logging.debug("Plotting MagLiteS footprint")
    infile = os.path.join(get_datadir(),'maglites-poly.txt')
    perim = np.loadtxt(infile,dtype=[('ra',float),('dec',float)])
    proj = safe_proj(bmap,perim['ra'],perim['dec'])
    bmap.plot(*proj,**kwargs)

def draw_maglites2(bmap,**kwargs):
    """
    Plot the MagLiteS Phase-II footprint.

    Parameters:
    -----------
    bmap   : The basemap object
    kwargs : Various plotting arguments

    Returns:
    --------
    None
    """
    # Plot the wide-field survey footprint
    logging.debug("Plotting footprint: %s"%opts.footprint)
    infile = os.path.join(get_datadir(),'maglitesII-poly.txt')
    perim = np.loadtxt(infile,dtype=[('ra',float),('dec',float),('poly',int)])
    for p in np.unique(perim['poly']):
        sel = (perim['poly'] == p)
        proj = safe_proj(bmap,perim[sel]['ra'],perim[sel]['dec'])
        bmap.plot(*proj,**kwargs)

def draw_bliss(bmap,**kwargs):
    """
    Plot the BLISS wide-field footprint.

    Parameters:
    -----------
    bmap   : The basemap object
    kwargs : Various plotting arguments

    Returns:
    --------
    None
    """
    # Plot the wide-field survey footprint
    logging.debug("Plotting footprint: %s"%opts.footprint)
    infile = os.path.join(get_datadir(),'bliss-poly.txt')
    perim = np.loadtxt(infile,dtype=[('ra',float),('dec',float),('poly',int)])
    for p in np.unique(perim['poly']):
        sel = (perim['poly'] == p)
        proj = safe_proj(bmap,perim[sel]['ra'],perim[sel]['dec'])
        bmap.plot(*proj,**kwargs)

def draw_decals(bmap,**kwargs):
    """
    Plot the DECaLS wide-field footprint.

    Parameters:
    -----------
    bmap   : The basemap object
    kwargs : Various plotting arguments

    Returns:
    --------
    None
    """
    # Plot the wide-field survey footprint
    logging.debug("Plotting footprint: %s"%opts.footprint)
    infile = os.path.join(get_datadir(),'decals-poly.txt')
    perim = np.loadtxt(infile,dtype=[('ra',float),('dec',float),('poly',int)])
    for p in np.unique(perim['poly']):
        sel = (perim['poly'] == p)
        proj = safe_proj(bmap,perim[sel]['ra'],perim[sel]['dec'])
        bmap.plot(*proj,**kwargs)


def plot(opts):
    """ 
    Core plotting function. Creates the basemap, overplots all of the
    requested features, and returns the map object.

    Parameters:
    -----------
    opts : command line options
    
    Returns:
    --------
    m : the basemap object
    """
    utc = parse_utc(opts.utc)
    CTIO.date = utc
    since = parse_since(opts.since)

    # Grab the data
    data = load_data(opts)

    # Subselect the data
    sel = np.in1d(data['filter'],FILTERS)
    if opts.band in FILTERS:
        sel &= (data['filter'] == opts.band)
    data = data[sel]

    expnum,telra,teldec,band = data['expnum'],data['telra'],data['teldec'],data['filter']

    # Set the colors
    if opts.color:
        nexp = len(expnum)
        ncolors = len(COLORS)
        color_repeat = np.repeat(COLORS.keys(),nexp).reshape(ncolors,nexp)
        color_idx = np.argmax(band==color_repeat,axis=0)
        color = np.array(COLORS.values())[color_idx]
    else:
        color = COLORS['none']

    # Select the exposure of interest
    if opts.expnum:
        match = np.char.array(expnum).endswith(str(opts.expnum))
        if not match.any():
            msg = "Exposure matching %s not found"%opts.expnum
            raise ValueError(msg)
        idx = np.nonzero(match)[0][0]
    elif len(data)==0:
        idx = slice(None)
    else:
        idx = 0

    # Create the figure
    if plt.get_fignums():
        fig,ax = plt.gcf(),plt.gca()
    else:
        fig,ax = plt.subplots(figsize=FIGSIZE,dpi=DPI)
        fig.canvas.set_window_title("DECam Pointings")
    #fig,ax = plt.subplots()

    # Zenith position
    lon_zen=lmst(CTIO); lat_zen = TEL_LAT
    # Create the Basemap
    proj_kwargs = PROJ[opts.proj]
    # Centering position
    if proj_kwargs['projection'] in basemap._pseudocyl:
        ### This should work, but doesn't.
        ### Compare lon_0=-80.58345277606 to lon_0=-80.6 or lon_0=-80.5
        #lon_0=lon_zen-360*(lon_zen>180),lat_zen=0
        lon_0,lat_0 = 0,0
    else:
        lon_0,lat_0 = -lon_zen, lat_zen # Center position

    proj_kwargs.update(lon_0=lon_0,lat_0=lat_0)

    m = basemap.Basemap(**proj_kwargs)
    def format_coord(x,y):
        #Format matplotlib cursor to display RA, Dec
        lon,lat = safe_proj(m,x,y,inverse=True)
        lon += 360*(lon < 0)
        return 'ra=%1.3f, dec=%1.3f'%(lon,lat)
    plt.gca().format_coord = format_coord

    parallels = np.arange(-90.,120.,30.)
    m.drawparallels(parallels)
    meridians = np.arange(0.,420.,60.)
    m.drawmeridians(meridians)
    for mer in meridians[:-1]:
        plt.annotate(r'$%i^{\circ}$'%mer,m(mer,5),ha='center')
    plt.annotate('West',xy=(1.0,0.5),ha='left',xycoords='axes fraction')
    plt.annotate('East',xy=(0.0,0.5),ha='right',xycoords='axes fraction')

    # markersize defined at minimum distortion point
    if proj_kwargs['projection'] in basemap._pseudocyl:
        x1,y1=ax.transData.transform(m(lon_0,lat_0+DECAM))
        x2,y2=ax.transData.transform(m(lon_0,lat_0-DECAM))
    else:
        x1,y1=ax.transData.transform(m(lon_zen,lat_zen+DECAM))
        x2,y2=ax.transData.transform(m(lon_zen,lat_zen-DECAM))

    # Since markersize defined in "points" in scales with figsize/dpi
    size = SCALE * (y1-y2)**2

    # Scale the marker size to the size of an exposure
    exp_zorder = 10
    exp_kwargs = dict(s=size,marker='H',zorder=exp_zorder,edgecolor='k',lw=1)

    # Projected exposure locations
    x,y = safe_proj(m,telra,teldec)

    # Plot exposure of interest
    if len(data):
        logging.debug("Plotting exposure: %i (%3.2f,%3.2f)"%(expnum[idx],telra[idx],teldec[idx]))
        # Hacked path effect (fix if matplotlib is updated)
        m.scatter(x[idx],y[idx],color='w',**dict(exp_kwargs,edgecolor='w',s=70,lw=2))
        m.scatter(x[idx],y[idx],color=color,**dict(exp_kwargs,alpha=1.0,linewidth=2))

    # Once matplotlib is updated
    #x = m.scatter(x[idx],y[idx],color=color,**exp_kwargs)
    #ef = patheffects.withStroke(foreground="w", linewidth=3)
    #x.set_path_effects([ef])

    # Plot previous exposures
    nexp_kwargs = dict(exp_kwargs)
    nexp_kwargs.update(zorder=exp_zorder-1,alpha=0.2,edgecolor='none')#,lw=0)

    exp_slice = slice(None,opts.numexp)
    numexp = len(x[exp_slice])
    logging.debug("Plotting last %s exposures"%(numexp))
    m.scatter(x[exp_slice],y[exp_slice],color=color[exp_slice],**nexp_kwargs)

    # Plot zenith position & focal plane scale
    zen_x,zen_y = m(lon_zen,lat_zen)
    #zen_kwargs = dict(color='green',alpha=0.75,lw=1,zorder=0)
    zen_kwargs = dict(color='green',alpha=0.75,lw=1,zorder=1000)
    if opts.zenith:
        logging.debug("Plotting zenith: (%.2f,%.2f)"%(lon_zen,lat_zen))
        m.plot(zen_x,zen_y,'+',ms=10,**zen_kwargs)
        logging.debug("Plotting focal plane scale.")
        m.tissot(lon_zen, lat_zen, DECAM, 100, fc='none', **zen_kwargs)

        # To test exposure size
        #m.tissot(lon_zen, lat_zen, DECAM, 100, fc='none', **zen_kwargs)
        #m.scatter(*m(lon_zen,lat_zen),**nexp_kwargs)
        #m.tissot(0, 0, DECAM, 100, fc='none', **zen_kwargs)
        #m.scatter(*m(0,0),**nexp_kwargs)


    # Plot airmass circle
    if opts.airmass < 1:
        logging.warning("Airmass must be greater than one.")
        opts.airmass = np.nan
    else:
        logging.debug("Plotting airmass: %s"%opts.airmass)
        angle = airmass_angle(opts.airmass)
        m.tissot(lon_zen, lat_zen, angle, 100, fc='none',**zen_kwargs)

    # Moon location and phase
    (moon_ra,moon_dec),moon_phase = moon(utc)
    if opts.moon:
        logging.debug("Plotting moon: %i%%,(%.1f,%.1f)"%(moon_phase,moon_ra,moon_dec))
        moon_txt = '%i%%'%moon_phase
        #bbox = dict(boxstyle='circle,pad=0.4',fc='k',ec='k',alpha=0.25,lw=2)
        moon_kwargs = dict(zorder=exp_zorder-1,fontsize=11,va='center',ha='center',weight='bold')
        ax.annotate(moon_txt,m(moon_ra,moon_dec),**moon_kwargs)
        # Again old matplotlib making things difficult
        moon_kwargs2 = dict(facecolor='k',alpha=0.25,lw=2,s=2000)
        ax.scatter(*m(moon_ra,moon_dec),**moon_kwargs2)

    # Plot footprint(s)
    fp_zorder=exp_zorder-1
    fp_kwargs=dict(marker='o',mew=0,mfc='none',color='k',lw=2,zorder=fp_zorder)
    if 'none' in opts.footprint:
        opts.footprint = ['none']
    if 'des' in opts.footprint:
        des_kwargs = dict(fp_kwargs,color='b')
        draw_des(m,**des_kwargs)
    if 'des' in opts.footprint or 'des-sn' in opts.footprint:
        sn_kwargs = dict(facecolor='none',edgecolor='b',projection=proj_kwargs['projection'],zorder=fp_zorder)
        draw_des_sn(m,**sn_kwargs)
    if 'smash' in opts.footprint:
        smash_kwargs = dict(facecolor='none',**exp_kwargs)
        smash_kwargs.update(zorder=exp_zorder+1)
        draw_smash(m,**smash_kwargs)
    if 'maglites' in opts.footprint:
        maglites_kwargs = dict(fp_kwargs,color='r')
        draw_maglites(m,**maglites_kwargs)
        draw_maglites2(m,**maglites_kwargs)
    if 'bliss' in opts.footprint:
        bliss_kwargs = dict(fp_kwargs,color='m')
        draw_bliss(m,**bliss_kwargs)
    if 'decals' in opts.footprint:
        decals_kwargs = dict(fp_kwargs,color='darkorange')
        draw_decals(m,**decals_kwargs)

    # Annotate with some information
    if opts.legend:
        logging.debug("Adding info text.")
        bbox_props = dict(boxstyle='round', facecolor='white')
        textstr= "%s %s\n"%("UTC:",utc.strftime('%Y-%m-%d %H:%M:%S'))
        if len(data):
            textstr+="%s %i (%s)\n"%("Exposure:",expnum[idx],band[idx])
        textstr+="%s %i\n"%("Num. Exp.:",numexp)
        textstr+="%s (%.1f$^{\circ}$, %.1f$^{\circ}$)\n"%("Zenith:",lon_zen,lat_zen)
        textstr+="%s %s\n"%("Airmass:",np.nan_to_num(opts.airmass))
        textstr+="%s %i%% (%.1f$^{\circ}$, %.1f$^{\circ}$)\n"%("Moon:",moon_phase,moon_ra,moon_dec)
        textstr+="%s %s"%("Footprint:",', '.join(opts.footprint))

        ax.annotate(textstr, xy=(0.90,1.05), xycoords='axes fraction',
                    fontsize=10,ha='left',va='top', bbox=bbox_props)

    # Plot filter legend
    if opts.color:
        logging.debug("Adding filter legend.")
        leg_kwargs = dict(scatterpoints=1,fontsize=10,bbox_to_anchor=(0.08,0.20))
        handles, labels = [],[]
        for k in FILTERS:
            if k == 'VR' and not (band=='VR').any(): continue
            labels.append(k)
            handles.append(plt.scatter(None,None,color=COLORS[k],**exp_kwargs))
        plt.legend(handles,labels,**leg_kwargs)

    # Plot the version number
    vers_kwargs = dict(xy=(0.985,0.015),ha='right',va='bottom',
                       xycoords='figure fraction',size=8)
    plt.annotate('pointing v.%s'%__version__,**vers_kwargs)

    # Plot the author's name
    auth_kwargs = dict(xy=(0.015,0.015),ha='left',va='bottom',
                       xycoords='figure fraction',size=8)
    plt.annotate(u'\u00a9'+' %s'%__author__,**auth_kwargs)

    return m

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('expnum',nargs='?',type=int,default=None,
                        help="exposure number to plot")
    parser.add_argument('-a','--airmass',default=1.4,type=float,
                        help='draw airmass limit')
    parser.add_argument('-b','--band',default='all',choices=BANDS,
                        help='draw exposures in specific band')
    parser.add_argument('-c','--color',default=True,type=boolean,
                        help='color corresponding to filter')
    parser.add_argument('--db',default='ctio',choices=['ctio','fnal'],
                        help='database to query for exposures')
    parser.add_argument('-f','--footprint',action='append',choices=FOOTPRINTS,
                        help='footprint to draw')
    parser.add_argument('--flavor',default='object|standard',type=str,
                        help='exposure type [object,zero,dome flat,etc.]')
    parser.add_argument('-i','--infile',default=None,
                        help='list of exposures to draw')
    parser.add_argument('--legend',default=True,type=boolean,
                        help='draw figure legend')
    parser.add_argument('-m','--moon',default=True,type=boolean,
                        help='draw moon location and phase')
    parser.add_argument('-n','--numexp',default=None,type=int,
                        help='number of most recent exposures to plot')
    parser.add_argument('-o','--outfile',default=None,
                        help='output file for saving figure')
    parser.add_argument('--propid',default=None,
                        help='draw exposures from specific propid')
    parser.add_argument('--proj',default='ortho',choices=PROJ.keys(),
                        help='projection for plot')
    parser.add_argument('--refresh',nargs='?',default=None,const=60,type=int,
                        help="refresh interval for figure (seconds).")
    parser.add_argument('--since',default=None,
                        help="UTC for first exposure (defaults to 12 hours)")
    parser.add_argument('--utc',default=None,
                        help="UTC for zenith position (defaults to 'now')")
    parser.add_argument('-v','--verbose',action='store_true',
                        help='output verbosity')
    parser.add_argument('--version',action='version',
                        version='%(prog)s '+__version__)
    parser.add_argument('-z','--zenith',default=True,type=boolean,
                        help="draw zenith position")

    opts = parser.parse_args()

    # Set logging level
    logging.basicConfig(level=logging.DEBUG if opts.verbose else logging.INFO,
                        format='%(message)s',stream=sys.stdout)

    if not opts.footprint: opts.footprint = ['des']
    
    # Do the plotting
    m = plot(opts)

    # In interactive session
    if sys.flags.interactive: plt.ion()

    if opts.outfile:
        # Save the figure
        logging.debug("Saving figure to: %s"%opts.outfile)
        plt.savefig(opts.outfile,dpi=250)
    elif not opts.refresh:
        # Show plot
        plt.show()
    else:
        # Refresh the plot
        plt.show(block=False)
        for i in range(MAXREF): # safer than while loop
            try: 
                plt.pause(opts.refresh)
            except TclError:
                # Catch the TclError thrown when window closed
                break
            logging.debug("Refreshing plot...")
            plt.cla()
            m = plot(opts)
        if i == MAXREF:
            logging.info("Reached max refresh number.")

#!/usr/bin/python
"""
Routines for comparing a set of lineaments to an NSR stress field.

"""

from . import lineament
from . import satstress
from . import stressplot
from pylab import *
from matplotlib import colors, colorbar
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import random
import os

# These are global variables.  I know, bad.  But this file is eventually going
# to get split into two pieces.  One for general purpose stuff, and another
# which is just a script generating figures for the NSR/TPW paper and my
# thesis.

# In general where to figs and other things get saved?
outdir = 'output'

# Where to put the figures for slurping up?
figdir = os.path.join(outdir,'figs')

# Where to store the giant tangle of pickled lineaments?
lindir = os.path.join(outdir,'lins')

# Format to save the figures in.  If None, don't save:
#save_fmt = None
save_fmt = 'pdf'
#save_fmt = 'png'

###############################################################################
# ANALYSES ####################################################################
###############################################################################
# These functions perform the analysis that makes up this paper.  Run from
# scratch, they start from the digital map, do the lineament fits and the TPW
# search work, and make all the figures for the paper and thesis on this topic.
# If run at high resolution, this will take days (or even weeks) on a desktop
# computer, but being able to force it all through from the beginning means we
# know there's nothing suspicious and leftover hiding in the shadows.
###############################################################################

def fromscratch(satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
                linfile="input/GlobalLineaments",\
                N_spin=1000, N_jack=100, N_pnp=100,\
                nb_hi=180, nb_lo=36, nlins=0,\
                pnp_lon=np.radians(80), pnp_lat=np.radians(10),\
                init_doppel_res=0.0, doppel_res=0.1, num_subsegs=10): #{{{
    """
    A "glue" function that strings together all of the analysis routines, so
    they can be run at one go.
    """
    from time import strftime
    from time import time
    from datetime import timedelta
    import pickle

    # All of the following are done by calcfits():
    #  - Read in the stress field definition (planet, satellite, materials, etc)
    #  - Read in the mapped features.
    #  - Create synthetic features based on the mapped features.
    #  - Calculate NSR fits for the mapped and synthetic features
    print("Beginning run at "+strftime('%Y-%m-%d %H:%M:%S'))
    t0 = time()
    maplins, synthlins, crazylins, gclins, tpwlins = calcfits(satfile=satfile, linfile=linfile, nb=nb_hi,\
                                                              pnp_lon=pnp_lon, pnp_lat=pnp_lat, nlins=nlins)

    NSR = maplins[0].stresscalc

    # Now we perform a bunch of different re-orientations of the mapped
    # dataset, re-fit them, and save the results for later plotting:

    # Even global coverage that will give us a distribution of RMDs:
    pnp_2lons, pnp_2lats = fibonacci_sphere(N_pnp*2)
    # limit the search to one hemisphere because of periodicity:
    pnp_lons = pnp_2lons[np.where(pnp_2lats > 0)]
    pnp_lats = pnp_2lats[np.where(pnp_2lats > 0)]

    devnull = tpw_polesearch(satfile=satfile, linfile=linfile, nb=nb_lo,\
                             pnp_lons=pnp_lons, pnp_lats=pnp_lats, nlins=nlins,\
                             poledir='tpw_polesearch')

    # Now we need to calculate the distribution of possible activity history
    # RMD values for completely randomized and longitudinally randomized
    # features.  In order for this to be a valid comparison, all of the
    # calculations should takae place at the same resolution, which has to be
    # pretty low or it will take forever:
    print("re-fitting maplins at lower resolution for jackstraws...")
    devnull = [ lin.calc_nsrfits(nb=nb_lo, stresscalc=NSR,\
                                 init_doppel_res=init_doppel_res,\
                                 doppel_res=doppel_res,\
                                 num_subsegs=num_subsegs) for lin in maplins ]

    spin_map_acthists = jackstraws(maplins, N=N_spin, spin=True, tpw=False)
    safe_pickle(spin_map_acthists, name='spin_map_acthists', dir=outdir)

    jack_map_acthists = jackstraws(maplins, N=N_jack, nb=nb_lo, spin=True, tpw=True)
    safe_pickle(jack_map_acthists, name='jack_map_acthists', dir=outdir)
    t1 = time()
    print("Finished run at "+strftime('%Y-%m-%d %H:%M:%S'))
    print("Elapsed time: "+ str(timedelta(seconds=t1-t0)))

#}}}

def deltasynth(satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
               linfile="input/GlobalLineaments", nb=180, nlins=500,\
               init_doppel_res=0.0, doppel_res=0.1, num_subsegs=10): #{{{
    """
    Read in the mapped features and calculate their fits for a variety of
    different values of NSR Delta.  Save the results.

    """

    europa = satstress.Satellite(open(satfile,'r'))
    NSR = satstress.NSR(europa)
    nsr_periods = europa.nsr_period*(np.logspace(-2, 1, num=7, base=10.0)/NSR.Delta())
    nsr_stresscalc_orig = satstress.StressCalc([satstress.NSR(europa),])

    labels = []
    lins_list = []
    for P_nsr in nsr_periods:
        europa.nsr_period = P_nsr
        NSR = satstress.NSR(europa)
        nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])
        label = "synth_Delta_%.3g" % (NSR.Delta(),)
        labels.append(label)

        print("Synthesizing %d NSR lineaments with Delta = %g" % (nlins, NSR.Delta()) )
        lins = random_nsrlins(nsr_stresscalc=nsr_stresscalc, nlins=nlins)

        print("Calculating fits for %s..." % (label,) )
        for lin,N in zip(lins,range(len(lins))):
            if (mod(N+1,50)==0):
                print("    N=%d" % (N+1,) )
            lin.calc_nsrfits(nb=nb, stresscalc=nsr_stresscalc_orig, init_doppel_res=init_doppel_res, doppel_res=doppel_res, num_subsegs=num_subsegs)

        print("Saving %s" % (label,) )
        safe_pickle(lins, name=label, dir=lindir)

        lins_list.append(lins)

    return(lins_list, labels)
#}}}

def deltafits(satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
              linfile="input/GlobalLineaments", nb=180, nlins=0,\
              init_doppel_res=0.0, doppel_res=0.1, num_subsegs=10): #{{{
    """
    Read in the mapped features and calculate their fits for a variety of
    different values of NSR Delta.  Save the results.

    """

    europa = satstress.Satellite(open(satfile,'r'))
    NSR = satstress.NSR(europa)
    nsr_periods = europa.nsr_period*(np.logspace(-2, 2, num=9, base=10.0)/NSR.Delta())

    lins_list = []
    for P_nsr in nsr_periods:
        europa.nsr_period = P_nsr
        NSR = satstress.NSR(europa)
        nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])

        label = "Delta_%.3g" % (NSR.Delta(),)
        lins = lineament.shp2lins(linfile, stresscalc=nsr_stresscalc)

        if nlins > 0:
            lins=lins[:nlins]

        print("Calculating fits for %s..." % (label,) )
        for lin,N in zip(lins,range(len(lins))):
            if (mod(N+1,60)==0):
                print("    N=%d" % (N+1,) )
            lin.calc_nsrfits(nb=nb, stresscalc=nsr_stresscalc, init_doppel_res=init_doppel_res, doppel_res=doppel_res, num_subsegs=num_subsegs)
        # only save if we're doing the whole dataset
        if nlins == 0:
            print("Saving %s" % (label,) )
            safe_pickle(lins, name=label, dir=lindir)

        lins_list.append(lins)

    return(lins_list)
#}}}

def calcfits(satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
             linfile="input/GlobalLineaments",
             nb=180, pnp_lon=radians(80), pnp_lat=radians(10), nlins=0,\
             init_doppel_res=0.0, doppel_res=0.1, num_subsegs=10): #{{{
    """
    Starting from scratch, read in the mapped features, calculate their fits,
    and use them to generate TPW features with the given pole, and to subsample
    the two sets of synthetic features (great circle segments and perfect NSR
    lineaments)

    """
    print("Initializing satellite")
    europa = satstress.Satellite(open(satfile,'r'))
    NSR = satstress.StressCalc([satstress.NSR(europa),])

    print("Loading mapped features")
    maplins = lineament.shp2lins(linfile, stresscalc=NSR)

    if nlins > 0:
        maplins=maplins[:nlins]

    print("Transforming to pre-TPW coordinates")
    tpwlins = [ lin.poleshift(pnp_lon=pnp_lon, pnp_lat=pnp_lat) for lin in maplins ]

    print("Generating randomized mapped dataset")
    crazylins = make_crazy(maplins, spin=True, tpw=True)

    print("Generating synthetic lineament pools")
    protolinlens = array([lin.length for lin in maplins])
    print("  great circle segments")
    gclins = linsample(proto=maplins, pool=random_gclins(nlins=10*len(maplins), minlen=protolinlens.min(), maxlen=protolinlens.max()))
    print("  model NSR lineaments")
    synthlins = linsample(proto=maplins, pool=random_nsrlins(nsr_stresscalc=NSR, nlins=10*len(maplins), minlen=protolinlens.min(), maxlen=protolinlens.max()))

    lins_list = (maplins, synthlins, crazylins, gclins, tpwlins)
    labels = ("map_nsrfit", "synth_nsrfit", "crazy_nsrfit", "gc_nsrfit", "tpw_nsrfit")

    # now calculate the fits...
    for lins,label in zip(lins_list, labels):
        print("Calculating fits for %s..." % (label,) )
        for lin,N in zip(lins,range(len(lins))):
            if (mod(N+1,60)==0):
                print("    N=%d" % (N+1,) )
            lin.calc_nsrfits(nb=nb, stresscalc=NSR, init_doppel_res=init_doppel_res, doppel_res=doppel_res, num_subsegs=num_subsegs)
        # only save if we're doing the whole dataset
        if nlins == 0:
            print("Saving %s" % (label,) )
            safe_pickle(lins, name=label, dir=lindir)

    return lins_list
#}}}

def tpw_polesearch(satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
                   linfile="input/GlobalLineaments",\
                   poledir="tpw_polesearch",\
                   nb=36, nlins=0, pnp_lons=None, pnp_lats=None): #{{{
    """
    Transform the lineaments read in from linfile to np different possible
    paleo north poles, specified by pnp_lons and pnp_lats. Fit them to the NSR
    stress field at nb different values of backrotation.  Save the results for
    later viewing.

    For testing purposes, set nlins to the number of features to fit.  Fits and
    transformed features will not be saved.

    Returns the relative mean differences (RMDs) associated with the mapped
    features transformed to those positions.

    """
    import pickle
    from time import strftime

    if pnp_lons is None or pnp_lats is None:
        pnp_lons, pnp_lats = fibonacci_sphere(100)

    print("Initializing satellite")
    europa = satstress.Satellite(open(satfile,'r'))
    NSR = satstress.StressCalc([satstress.NSR(europa),])

    print("Loading mapped features")
    maplins = lineament.shp2lins(linfile, stresscalc=NSR)

    if nlins > 0:
        maplins = maplins[:nlins]

    # remove link to newest pre-existing directory if it exists:
    try:
        os.unlink(os.path.join(outdir,poledir))
    except OSError:
        pass

    newdir = poledir+'_'+strftime('%Y%m%d%H%M%S')
    os.mkdir(os.path.join(outdir,newdir))
    print "linking %s to %s" % (newdir, os.path.join(outdir,poledir))
    os.symlink(newdir, os.path.join(outdir,poledir))

    tpwlins_RMDs = []
    for pnp_lon, pnp_lat, N in zip(pnp_lons, pnp_lats, arange(len(pnp_lons))):
        print("Fitting paleopole %d / %d (lon=%f, lat=%f)" % (N+1,len(pnp_lons), degrees(pnp_lon), degrees(pnp_lat)) )
        tpwlins = [ lin.poleshift(pnp_lon=pnp_lon, pnp_lat=pnp_lat) for lin in maplins ]

        label = "tpw_lon%.6f_lat%.6f_nsrfit" % (degrees(pnp_lon), degrees(pnp_lat))
        devnull = [ lin.calc_nsrfits(nb=nb, stresscalc=NSR, init_doppel_res=0.0, doppel_res=0.1, num_subsegs=10) for lin in tpwlins ]
        tpwlins_RMDs.append(acthist_RMD(tpwlins))

        # only save the features if we're doing the whole dataset:
        if nlins == 0:
            pickle.dump(tpwlins, open(os.path.join(outdir,poledir,label+'.pkl'),'w'))
    
    return(tpwlins_RMDs)
#}}}

###############################################################################
#  Analysis Helper Functions:
###############################################################################

def random_nsrlins(nsr_stresscalc=None, nlins=1000, minlen=0.0, maxlen=1.25): #{{{
    """

    Create nlins lineament objects, resulting from tensile fracture under the
    NSR stress field, having lengths between minlen and maxlen (radians of
    arc), with their locations and (maximum) lengths randomly distributed.

    Because some lineaments will self-terminate upon reaching compressive
    zones, this will naturally result in a bias towards shorter lineaments.
    This bias can be compensated for later by subsampling the synthetic
    lineament population according to whatever criteria the user desires.

    """

    if nsr_stresscalc is None:
        satfile   = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
        europa = satstress.Satellite(open(satfile,'r'))
        nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])

    linlist = []

    N=0
    max_lengths = minlen+(rand(2*nlins)*(maxlen-minlen))
    init_lons, init_lats = lineament.random_lonlatpoints(2*nlins)
    while len(linlist) < nlins:
        seg_len=min(0.01, max_lengths[N]/10.0)
        newlin = lineament.lingen_nsr(nsr_stresscalc, init_lon=init_lons[N], init_lat=init_lats[N], max_length=max_lengths[N], prop_dir='both', seg_len=seg_len, num_subsegs=2)
        if newlin.length > minlen:
            linlist.append(newlin)
        N = N+1

    return linlist
#}}}

def random_gclins(nlins=1000, minlen=0.0, maxlen=1.25): #{{{
    """
    Create nlins lineament objects, whose paths approximate great circle
    segments on the surface of the satellite, with lengths ranging from
    minlen to maxlen (radians of arc), randomly distributed in length,
    orientation, and location on the surface.

    """

    # - pick a random endpoint v1 on the surface of the sphere
    init_lons, init_lats = lineament.random_lonlatpoints(nlins)
    # - pick a random azimuth, az
    azimuths = 2*pi*rand(nlins)
    # - pick a random length L between minlen and maxlen
    linlengths = (minlen+(rand(nlins)*(maxlen-minlen)))
    # - calculate the location of endpoint v2, L radians away
    #   from v1, along the initial heading az, from v1.
    fin_lons, fin_lats = lineament.spherical_reckon(init_lons, init_lats, azimuths, linlengths)
    # - use lingen_greatcircle() to calculate intermediary vertices.
    return([ lineament.lingen_greatcircle(init_lons[N], init_lats[N], fin_lons[N], fin_lats[N]) for N in range(nlins) ])

#}}}

def linresample_byN(lins, fraction=1.0): #{{{
    """
    Create a synthetic lineament map which is similar to the one represented by
    the set of lineaments 'lins' that was passed in.

    The new dataset is generated by a process which is equivalent to laying all
    of the lineaments in the input set end to end, and then choosing a random
    distance between 0 and the sum of the lengths of the input dataset, and
    then picking a copy of that lineament which that distance lies within in
    the series of lineaments laying end to end.  This process is repeated until
    the new dataset has a length that is greater than the total length of the
    input dataset, minus half the median length of a feature in that dataset.

    """

    return(array(lins)[np.random.randint(len(lins), size=np.int(fraction*len(lins)))])
#}}} end linresample

def linresample_byL(lins, fraction=1.0): #{{{
    """
    Create a synthetic lineament map which is similar to the one represented by
    the set of lineaments 'lins' that was passed in.

    The new dataset is generated by a process which is equivalent to laying all
    of the lineaments in the input set end to end, and then choosing a random
    distance between 0 and the sum of the lengths of the input dataset, and
    then picking a copy of that lineament which that distance lies within in
    the series of lineaments laying end to end.  This process is repeated until
    the new dataset has a length that is greater than the total length of the
    input dataset, minus half the median length of a feature in that dataset.

    """

    linlengths = [ lin.length for lin in lins ]

    # Create a lookup table in which each input lineament has a number of of
    # entries proportional to its length, with sufficient resolution that even
    # the shortest features have several entries.
    mincount   = 10 # even the shortest feature should get this many counts
    count_mult = mincount/np.min(linlengths)
    lookup_list = []

    for i in range(len(linlengths)):
        lookup_list += [i,]*np.int(mincount*linlengths[i])
    lookup_arr = np.array(lookup_list)

    resampled_lins = []
    total_length = 0.0
    while total_length < fraction*(np.sum(linlengths)-np.median(linlengths)):
        resampled_lins.append(lins[lookup_arr[np.random.randint(len(lookup_arr))]])
        total_length += resampled_lins[-1].length
    
    return(resampled_lins)
#}}} end linresample

def linsample(proto=None, pool=None, nbins=20, fraction=1.0): #{{{
    """
    Select a set of lineaments from a pool with a similar length distribution
    to that of a prototype set.

    nbins controls the number of discrete, evenly spaced length bins which are
    to be approximated in replicating the length distribution.  The lower edge
    of the shortest bin is the length of the shortest feature in the prototype
    set, and similarly the high edge of the longest bin is the length of the
    longest feature in the prototype set.

    """

    # define the length bins to replicate:
    proto_lengths = array([ lin.length for lin in proto ])
    proto = array(proto)
    pool_lengths  = array([ lin.length for lin in pool  ])
    pool = array(pool)

    bins = linspace(min(proto_lengths)-0.00001, max(proto_lengths)+0.00001, num=nbins+1)
    sample_lins = []

    for i in range(nbins):
        # Organize the prototypes and pool into bins by length
        proto_binned_lins    = proto[where( logical_and(proto_lengths >= bins[i], proto_lengths < bins[i+1]) )]
        proto_binned_lengths = array([lin.length for lin in proto_binned_lins])

        pool_binned_lins    = pool[where(  logical_and(pool_lengths >= bins[i], pool_lengths < bins[i+1]) )]
        pool_binned_lengths = array([lin.length for lin in pool_binned_lins])

        if len(pool_binned_lins) > 0:
            sample_binned_lins = pool_binned_lins[ np.random.randint(0, high=len(pool_binned_lins), size=int(fraction*len(proto_binned_lins))) ]
            sample_binned_lengths = array([lin.length for lin in sample_binned_lins])
            sample_lins.append(sample_binned_lins)

        #print("%3.d: %7.3g %7.3g " % (i+1, np.mean(proto_binned_lengths), np.mean(sample_binned_lengths)) )

    return(concatenate(sample_lins))
#}}}

def safe_pickle(thing, name=None, dir=None): #{{{ TODO: find all safe_pickle and switch them
    """
    Pickle something without overwriting an old thing having the same name, by
    renaming the pickle file with a datetime string, and using a symlink to
    point at the most recent version.

    """
    from time import strftime
    import pickle
    import os

    assert name is not None
    assert dir is not None

    outfile="%s_%s.pkl" % (name, strftime('%Y%m%d%H%M%S'))
    outpath = os.path.join(dir,outfile)
    linkpath = os.path.join(dir,name)
    pickle.dump(thing, open(outpath,'w'))

    try:
        os.unlink(linkpath)
    except OSError:
        pass

    os.symlink(outfile, linkpath)

#}}} end safe_pickle()

def load_lins(linpath): #{{{
    """
    Load a list of lineaments from either a pickle file, or a directory filled
    with pickle files (as is produced by calc_fits)

    """
    import pickle

    linlist = []
    if os.path.isdir(linpath):
        lindir = os.listdir(linpath)
        for linfilename in lindir:
            linfilepartname = os.path.join(linpath, linfilename)
            print "Loading", linfilepartname
            newpart = pickle.load(open(linfilepartname))
            print "  found %d lineaments" % (len(newpart,))
            linlist += newpart
    elif os.path.isfile(linpath):
        linlist = pickle.load(open(linpath))
    else:
        raise os.error("Path: \'%s\' is not a file or directory" % (linpath,) )

    return linlist
#}}}

def reload_nsrfits(update=False): #{{{
    """
    Loads the lineaments which I most often end up fooling around with and
    returns them as a list, ordered as follows:

    maplins    Europa's lineaments, as mapped
    synthlins  A perfect synthetic dataset
    crazylins  The mapped lineaments, with random locations/orientations
    gclins     A collection of great circle segments
    tpwlins    The mapped lineaments, transformed to PNP 80E 10N

    maplins, synthlins, crazylins, gclins, tpwlins = nsrhist.reload_nsrfits()

    if update is True, then these datasets are read in, updated to reflect the
    most recent version of the Lineament object, saved back to disk, and then
    returned.

    """

    maplins   = load_lins(os.path.join(lindir,'map_nsrfit'))
    synthlins = load_lins(os.path.join(lindir,'synth_nsrfit'))
    crazylins = load_lins(os.path.join(lindir,'crazy_nsrfit'))
    gclins    = load_lins(os.path.join(lindir,'gc_nsrfit'))
    tpwlins   = load_lins(os.path.join(lindir,'tpw_nsrfit'))

    if update is True:
        maplins   = lineament.update_lins(maplins)
        synthlins = lineament.update_lins(synthlins)
        crazylins = lineament.update_lins(crazylins)
        gclins    = lineament.update_lins(gclins)
        tpwlins   = lineament.update_lins(tpwlins)

        safe_pickle(synthlins, name='synth_nsrfit', dir=lindir)
        safe_pickle(maplins,   name='map_nsrfit',   dir=lindir)
        safe_pickle(crazylins, name='crazy_nsrfit', dir=lindir)
        safe_pickle(gclins,    name='gc_nsrfit',    dir=lindir)
        safe_pickle(tpwlins,   name='tpw_nsrfit',   dir=lindir)

    return(maplins, synthlins, crazylins, gclins, tpwlins)
#}}}

def load_tpw_RMDs(tpwdir="tpw_polesearch", fast=True): #{{{
    """
    Look in the specified output directory for sets of lineaments which have
    been transformed to a paleopole location, and fit to the NSR field.  Parse
    their filenames to glean the paleopole location, and return three lists:

    (pnp_lons, pnp_lats, pnp_RMDs)

    representing the locations of the paleopoles and the RMD of the inferred
    activity histories associated with them.

    """
    from re import match
    from os import listdir
    import os.path
    import pickle

    # Check to see if we've read these things in before and saved the
    # results.
    try:
        tpw_polesearch_results = pickle.load(open(os.path.join(outdir,tpwdir+'.pkl')))
        pnp_lons = tpw_polesearch_results['pnp_lon']
        pnp_lats = tpw_polesearch_results['pnp_lat']
        acthist_RMDs = tpw_polesearch_results['acthist_RMD']
    except IOError:
    # If we get here, then we were told to do it fast, but can't because
    # the results haven't been cached.  Do it slow instead:
        print("Could not read cached TPW polesearch results.  Loading from scratch.")
        fast=False

    if fast is True:
        # we've already read these things in once.... don't need to do it again.
        pass
    else:
        acthist_RMDs = []
        pnp_lons = []
        pnp_lats = []
        # for each filename in the output directory
        outfiles = os.listdir(os.path.join(outdir,tpwdir))
        for n,fn in zip(range(len(outfiles)), outfiles):
            # see if it matches our output filename format
            filematch = match(r'tpw_lon(?P<pnp_lon>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)_lat(?P<pnp_lat>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)_nsrfit.pkl$', fn)
            # if it does, extract information from the filename, load the lineaments:
            if filematch is not None:
                pnp_lon, pnp_lat = float(filematch.group('pnp_lon')), float(filematch.group('pnp_lat'))
                print("TPW dataset %d / %d: pnp_lon=%.6f, pnp_lat=%.6f" % (n+1, len(outfiles), pnp_lon, pnp_lat))
                tpwlins = load_lins(os.path.join(outdir,tpwdir,fn))
                pnp_lons.append(radians(pnp_lon))
                pnp_lats.append(radians(pnp_lat))
                acthist_RMDs.append(acthist_RMD(tpwlins))

        # Now that we've gone to all the trouble to load this thing in, let's
        # save it for next time:
        dtype = [('pnp_lon',float),('pnp_lat',float),('acthist_RMD',float)]
        pickle(array([ (pnp_lons[n], pnp_lats[n], acthist_RMDs[n]) for n in range(len(pnp_lons))], dtype=dtype), open(os.path.join(outdir,tpwdir+'.pkl'),'w'))

    return(array(pnp_lons), array(pnp_lats), acthist_RMDs)
#}}}

def calc_acthist(lins, dbar_max=0.125, norm_length=None): #{{{
    if norm_length is None:
        norm_length = sum([lin.length for lin in lins])
    return(array([ lin.length*lin.nsrfits(dbar_max=dbar_max) for lin in lins ]).sum(axis=0)/norm_length)

def acthist_amplitude(lins, dbar_max=0.125, norm_length=None):
    acthist = calc_acthist(lins, dbar_max=dbar_max, norm_length=norm_length)
    return(max(acthist) - min(acthist))

def acthist_meandiff(lins, dbar_max=0.125, norm_length=None):
    acthist = calc_acthist(lins, dbar_max=dbar_max, norm_length=norm_length)
    return(np.fabs(np.tile(acthist,len(acthist)) - np.repeat(acthist,len(acthist))).mean())

def acthist_RMD(lins, dbar_max=0.125, norm_length=None):
    acthist = calc_acthist(lins, dbar_max=dbar_max, norm_length=norm_length)
    return(np.fabs(np.tile(acthist,len(acthist)) - np.repeat(acthist,len(acthist))).mean()/acthist.mean())

#}}}

def jackstraws(lins, stresscalc=None, N=10, nb=36, spin=True, tpw=False): #{{{
    """
    Take the given set of lineaments, resample them and scatter them over the
    surface of Europa.  If they've been completely randomized, also calculate
    their activity history (at a resolution specified by nb).  Repeat this
    process N times and return the length N list of activity histories.

    If spin is true, randomize their longitudes.  If tpw is true randomize
    their location and orientations fully.

    """
    if stresscalc is None:
        stresscalc = lins[0].stresscalc

    acthist_list = []
    for n in range(N):
        if np.mod(n+1,100) == 0:
            print("Tossing jackstraws %d / %d" % (n+1, N,))
        crazylins = make_crazy(linresample_byN(lins), tpw=tpw, spin=spin)

        if tpw is True:
            print("Fitting jackstraws %d / %d" % (n+1, N,))
            [ lin.calc_nsrfits(nb=nb, stresscalc=stresscalc) for lin in crazylins ]

        acthist_list.append(calc_acthist(crazylins))

    return(acthist_list)
#}}}

def nsrfit_Q(lins, use_stress=True, dbar_max=0.125): #{{{
    """
    Takes a list of lineaments which have had their NSR fits calculated and
    returns a the sum of the lineament lengths multiplied by their greatest
    value of f_{nsr}(b), divided by the overall length of the dataset.

    """

    return(sum([lin.length*lin.best_fit(use_stress=use_stress, dbar_max=dbar_max)[0] for lin in lins])/sum([lin.length for lin in lins]))

#}}}

def fibonacci_sphere(N, halfsphere=False): #{{{
    """
    Generate N points on the surface of the sphere, fairly evenly spaced from
    one another.  If halfsphere is True, return only points with longitudes
    between 0 and 180 degrees.

    """

    # if we're only doing half the sphere, need to double N first to end up
    # with the right number of points
    if halfsphere is True:
        N = 2*N

    inc = pi * (3 - sqrt(5))
    off = 2. / N
    k = arange(0,N)
    y = k*off - 1. + 0.5*off
    r = sqrt(1 - y*y)
    phi = k * inc
    x = cos(phi)*r
    z = sin(phi)*r
    theta = arctan2(sqrt(x**2+y**2),z)
    lons = arctan2(y,x)
    lats = pi/2-theta

    # Only retain those points with longtiudes between 0 and 180 degrees.
    if halfsphere is True:
        lats = lats[where(mod(lons,2*pi) < pi)]
        lons = lons[where(mod(lons,2*pi) < pi)]

    return lons, lats
#}}}

def make_crazy(lins, tpw=True, spin=True): #{{{
    """
    Take a list of lineaments and randomize their locations and orientations on
    the sphere by applying random TPW and NSR shifts to them.

    If tpw is True, re-orient and move the features on the surface arbitrarily.
    This results in the loss of any previously calculated fit information.

    if spin is True, simply shift them in longitude.  If only this
    transformation is performed, the fit information is preserved.

    """

    newlins = []
    if spin is True:
        rand_bs = 2*pi*np.random.random(len(lins))
        for lin, b in zip(lins, rand_bs):
            newlins.append(lin.lonshift(b))
        lins = newlins

    if tpw is True:
        newlins = []
        rand_lons, rand_lats = lineament.random_lonlatpoints(len(lins))
        for lin, pnp_lon, pnp_lat in zip(lins, rand_lons, rand_lats):
            newlins.append(lin.poleshift(pnp_lon=pnp_lon, pnp_lat=pnp_lat))

    return(newlins)
#}}}

def synthpeak(lins, N=100, peak_frac=0.4, scale=np.radians(15)): #{{{
    """
    Create a N synthetic datasets, each composed of peak_frac parts de-rotated
    features, subsampled from lins, and (1-peak_frac) parts longitudinally
    randomized resampling of lins Return the list of lists of lineaments.

    """

    map_acthist = calc_acthist(lins)
    mu = lins[0].bs[int(where(np.abs(map_acthist-np.max(map_acthist)) < 1e-6)[0])]
    maps = []
    for i in range(N):
        random_part = make_crazy(linresample_byN(lins, fraction=(1.0-peak_frac)),tpw=False)
        derotated_part = linresample_byN(lins, fraction=peak_frac)
        derotated_part = [ lin.lonshift(-b+lin.best_fit()[1]) for lin,b in zip(derotated_part,np.random.normal(loc=mu, scale=scale, size=len(derotated_part))) ]
        maps.append(concatenate([random_part,derotated_part]))

    return(mu, maps)
#}}}

def linlatstats(lins): #{{{
    """
    returns a triple (linlons, linlats, linlens) corresponding to the
    longitudes and latitudes of the midpoings of the segments making up the
    lineaments lins, and their lengths (in radians of arc).

    """

    linlats = np.array([])
    linlons = np.array([])
    linlens = np.array([])
    for lin in lins:
        linlens = np.concatenate([linlens,lin.seg_lengths()])
        mp_lons, mp_lats = lin.seg_midpoints()
        mp_lons = np.mod(mp_lons,2.0*np.pi)
        linlons = np.concatenate([linlons,mp_lons])
        linlats = np.concatenate([linlats,mp_lats])

    return(linlons, linlats, linlens)
    #}}}

###############################################################################
# PLOTTING
###############################################################################

def makefigs(dbar_max=0.125, numhists=100, all=False, stress=False, maps=False, lindensity=False, examples=False, hists=False, stats=False, tpw=False): #{{{
    """
    A wrapper that takes input parameters for the figures, and calls all of the
    individual figure plotting routines, below.

    """
    from osgeo import gdal
    from os import path

    # If we want to re-make all of the figures:
    if all is True:
        stress     = True # works w/ fromscratch()
        maps       = True # works w/ fromscratch()
        lindensity = True # works w/ fromscratch()
        examples   = True # works w/ fromscratch()
        hists      = True # works w/ fromscratch()
        stats      = True # works w/ fromscratch() but R^2(min(Dbar),S) w lat > 30 deg is very small...
        tpw        = True # works w/ fromscratch()

    maplins, synthlins, crazylins, gclins, tpwlins = reload_nsrfits()

    #if stress is True: #{{{2
        #for Delta_nsr in np.logspace(-2,2,9):
        #    StressFieldNSR(Delta_nsr=Delta_nsr)
        #for orbital_pos in np.linspace(0,360,24,endpoint=False):
        #    StressFieldDiurnal(orbital_pos=orbital_pos, scale=2e6)

    #}}}2

    if maps is True: #{{{2
        print("Plotting Mapped Lineaments")
        FitMap(maplins, nbins=18, titlestr="Global Lineaments as Mapped", dbar_max=dbar_max, outfile='FitMap_Mapped')

        print("Plotting Mapped Lineaments, relative to stress field")
        FitMap(maplins, nbins=18, titlestr="Mapped Lineaments Backrotated for Best Fit", dbar_max=dbar_max, derotate=True, showbad=False, outfile='FitMap_Derotated')

        print("Plotting Randomized Mapped Lineaments")
        FitMap(crazylins, nbins=18, titlestr="Randomized Global Lineaments", dbar_max=dbar_max, outfile='FitMap_Crazy')

        print("Plotting Derotated Randomized Lineaments")
        FitMap(crazylins, nbins=18, titlestr="Derotated Randomized Lineaments", dbar_max=dbar_max, outfile='FitMap_CrazyDerotated', showbad=False, derotate=True)

        print("Plotting Pre-TPW Lineaments, fit to NSR stresses")
        FitMap(tpwlins, nbins=18, titlestr="Mapped Lineaments before TPW w/ PNP=80E10N", dbar_max=dbar_max, outfile='FitMap_PreTPW')

        print("Plotting random great circle segments")
        FitMap(gclins, nbins=18, titlestr="Random Great Circle Segments", dbar_max=dbar_max, outfile='FitMap_GreatCircle')

        print("Plotting Derotated Great Circle Segments")
        FitMap(gclins, nbins=18, titlestr="Random Great Circle Segments Backrotated for Best Fit", dbar_max=dbar_max, derotate=True, showbad=False, outfile='FitMap_DerotatedGC')

        print("Plotting synthetic NSR lineaments")
        FitMap(synthlins, titlestr="Perfect Synthetic NSR Lineaments for b=0", dbar_max=dbar_max, outfile='FitMap_SyntheticNSR')
    #}}}2

    if lindensity is True: #{{{2
        # This shows the density of the mapped features, and a bunch of information
        # about how that relates to the spacecraft observations.
        LinDensityMap(maplins, maxdist=250, label="Density of Mapped Features")
        LinLatStatsCompare(maplins)
    #}}}2

    if examples is True: #{{{2
        # List of canonical features from the map: {{{3
        good_nsr1    = maplins[1]   # 1500 km long, symmetric arcuate lineament in the N. hemisphere
        good_nsr2    = maplins[53]  #  700 km long, asymmetric arcuate lineament in the N. hemisphere
        good_nsr3    = maplins[108] # 1800 km long, nearly symmetric arcuate lineament in the S. hemisphere
        bad_nsr1     = maplins[20]  #  553 km long, north-south lineament in the S. hemisphere
        bad_nsr2     = maplins[6]   # 1222 km long, diagonal, grazing the equator
        bad_nsr3     = maplins[36]  # 1175 km long, diagonal, equator crossing
        bad_nsr4     = maplins[54]  # 1036 km long, north-south, equator crossing
        bad_nsr5     = maplins[70]  # 1061 km long, one of the SCDs, crossing the equator
        bad_nsr6     = maplins[120] #  640 km, N-S equator crossing
        bad_nsr7     = maplins[122] # 1300 km, N-S equator crossing
        cycloid1     = maplins[112] # 1132 km long cycloid, 4-5 arcs, near 30S.  Low dbar, high delta
        cycloid2     = maplins[137] #  458 km long cycloid, 5 arcs, near 60S.  Low dbar, high delta
        cycloid3     = maplins[148] #  776 km long cycloid, 5 arcs, 45-70N.
        cycloid4     = maplins[155] # 1711 km long semi-cycloidal, actually passed the test
        cycloid5     = maplins[159] # 1527 km long semi-cycloidal, actually passed the test
        sinuous1     = maplins[23]  # 1334 km long, sinuous from 30N to 75N
        sinuous2     = maplins[136] # 1189 km long, sinuous from 30S to 70S
        short_hilat1 = maplins[11]  #  450 km, east-west, just above 30N... fits perfectly
        short_hilat2 = maplins[78]  #  200 km, diagonal, 50-60N, fits perfectly
        short_hilat3 = maplins[11]  #  183 km, east-west, just above 30N... fits perfectly
        short_hilat4 = maplins[160] #  200 km, diagonal, 75S... fits very well
        short_lolat1 = maplins[26]  #  500 km, diagonal, between 0N and 30N, fits perfectly
        short_lolat2 = maplins[41]  #  177 km, diagonal, between 0N and 30N, does not fit because of dbar
        short_lolat3 = maplins[43]  #  197 km, diagonal, between 0N and 30N, fits very well
        dbar_fail1   = maplins[7]   # 1073 km, arcuate, dbar_min doesn't quite line up with delta_min, and so it gets lost.
        dbar_fail2   = maplins[62]  #  500 km, east-west almost at the equator.  Large number of "good fits" possible.
        dbar_fail3   = maplins[63]  #  500 km, east-west almost at the equator.  Large number of "good fits" possible.
        dbar_fail4   = maplins[115] # 1262 km, looks good, but dbar just doesn't quite get low enough... Grr.
        #}}}3

        eg_lins = [ good_nsr1, bad_nsr2, cycloid1, sinuous1, sinuous2, short_hilat1, short_lolat2, dbar_fail2 ]

        eg_labels = [ "Good NSR", "Bad NSR", "Cycloial", "Sinuous", "Sinuous",\
                      "Short straight mid-latitude", "Short straight low-latitude",\
                      "Poorly constrained best fit" ]

        FitCurveExamples(eg_lins, labels=eg_labels, dbar_max=dbar_max)

    #}}}2

    if hists is True: #{{{2
        # These three figures demonstrate robustness of the metric
        ActHist_ByLength(maplins, dbar_max=dbar_max, norm_by_all=True)
        ActHist_BySin(maplins, dbar_max=dbar_max, norm_by_all=True)
        ActHist_ByDbar(maplins, norm_by_all=False)

        # Statistics of the dataset, and an attempt at synthesizing our map
        print("Calculating activity histories for %d resampled maps" % (numhists,) )
        ActHist_MapStats(maplins, N=numhists)

        print("Calculating activity histories for %d spin cycle maps" % (numhists,) )
        ActHist_SpinCycleStats(maplins, N=numhists)

        print("Calculating activity histories for %d partially derotated maps" % (numhists,) )
        ActHist_PeakStats(maplins, N=numhists, scale=np.radians(15))

        print("Calculating activity histories for %d synthetic map recreations" % (numhists,) )
        ActHist_Synthesized(maplins, N=numhists, peak_frac=0.4, scale=np.radians(15))

        # TODO:
        #print("Calculating activity histories for %d jackstraws maps" % (numhists,) )
        #ActHist_JackStrawsStats(maplins, N=numhists)

        print("Calculating activity histories for %d resampled great circle maps" % (numhists,) )
        ActHist_GreatCircleStats(gclins, maplins, N=numhists)

        #print("Calculating activity histories for %d derotated maps" % (numhists,) )
        #ActHist_Derotated(maplins, synthlins, crazylins, N=numhists)

        print("Plotting distribution of activity history uniformity")
        ActHist_RMDProb(maplins, tpwlins)
    #}}}2
 
    if stats is True: #{{{2
        LinLengthDist(maplins, label="Lineaments as Mapped",             outfile='LinLenDist_Mapped')
        LinLengthDist(gclins,  label="Random Great Circle Segments",     outfile='LinLenDist_GC')
        LinLengthDist(tpwlins, label="pre-TPW (pole=80E10N) Lineaments", outfile='LinLenDist_TPW')

        DbarLengthCorr(gclins)
        DbarSinuosityCorr(maplins)

        QvDbarMax(synthlins, gclins, crazylins, maplins, nq=101, use_stress=True, outfile='QvDbarMax_Wstress')
        QvDbarMax(synthlins, gclins, crazylins, maplins, nq=101, use_stress=False, outfile='QvDbarMax_NoStress')
    #}}}2

    if tpw is True: #{{{2
        # Plots having to do with the search for a best TPW pole...
        ActHist_RMDMap()

    #}}}2

# end makefigs }}}

def StressFieldNSR(scalar_field='tens', cmap=cm.gray, plot_tens=True, plot_comp=True, plot_lesser=True, plot_greater=True, Delta_nsr=None, scale=None): #{{{

    min_lon = np.radians(  0.0)
    max_lon = np.radians(180.0)
    min_lat = np.radians(-90.0)
    max_lat = np.radians( 90.0)

    scalar_nlons = 300
    scalar_nlats = scalar_nlons*(max_lat-min_lat)/(max_lon-min_lon)

    satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = satstress.Satellite(open(satfile,'r'))
    NSR = satstress.NSR(europa)
    nsr_stresscalc = satstress.StressCalc([NSR,])

    # We want to be able to pick a Delta value and have it shown
    if Delta_nsr is not None:
        europa.nsr_period *= Delta_nsr/NSR.Delta()
        NSR = satstress.NSR(europa)
        nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])
    else:
        Delta_nsr = NSR.Delta()

    the_fig = figure(figsize=(16,12))
    map_ax  = the_fig.add_subplot(1,1,1)
    map_ax.set_title("Tidal Stresses due to NSR ($\Delta_{nsr}=$%.3g)" % (Delta_nsr,))

    grid_basemap = Basemap(llcrnrlon = np.degrees(min_lon),\
                             llcrnrlat = np.degrees(min_lat),\
                             urcrnrlon = np.degrees(max_lon),\
                             urcrnrlat = np.degrees(max_lat),\
                                    ax = map_ax)

    field_data = stressplot.scalar_grid(stresscalc=nsr_stresscalc, nlons=scalar_nlons,\
                                        nlats=scalar_nlats, min_lon=min_lon,\
                                        max_lon=max_lon, min_lat=min_lat,\
                                        max_lat=max_lat, field=scalar_field,\
                                        cmap=cmap, basemap_ax=grid_basemap)
    vector_nlons = 13
    vector_nlats = 13

    vector_lons  = np.linspace(min_lon, max_lon, vector_nlons)
    vector_lats  = np.linspace(min_lat, max_lat, vector_nlats)
    vector_mesh_lons, vector_mesh_lats = np.meshgrid(vector_lons, vector_lats)
    vector_mesh_lons = np.ravel(vector_mesh_lons)
    vector_mesh_lats = np.ravel(vector_mesh_lats)

    if scale is None:
        nsr_scale = max(np.fabs(nsr_stresscalc.mean_global_stressmag()))*50
    else:
        nsr_scale = scale

    stressplot.vector_points(stresscalc=nsr_stresscalc, lons=vector_mesh_lons, lats=vector_mesh_lats,\
                             time_t=0.0, basemap_ax=grid_basemap, arrow_width=0.002,\
                             plot_tens=plot_tens, plot_comp=plot_comp,\
                             plot_lesser=plot_lesser, plot_greater=plot_greater, scale=nsr_scale)

    grid_basemap.drawmeridians(np.degrees(np.linspace(min_lon, max_lon, 13)), labels=[1,0,0,1], linewidth=0.1, color='gray')
    grid_basemap.drawparallels(np.degrees(np.linspace(min_lat, max_lat, 13)), labels=[1,0,0,1], linewidth=0.1, color='gray')
    grid_basemap.drawmapboundary()

    # Need some kind of scale bar for the scalar field in the background:
    cb_ax,kw = colorbar.make_axes(map_ax, orientation='vertical')
    colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=colors.Normalize(vmin=np.min(field_data), vmax=np.max(field_data)), orientation='vertical', format='%.2g')
    cb_ax.set_ylabel(scalar_field)
    draw()

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,'StressFieldNSR_Delta%.3g.%s' % (Delta_nsr,save_fmt)))

#}}} end StressFieldNSR

def StressFieldDiurnal(scalar_field='tens', cmap=cm.gray, plot_tens=True, plot_comp=True, plot_lesser=True, plot_greater=True, orbital_pos=0.0, scale=None): #{{{

    min_lon = np.radians(  0.0)
    max_lon = np.radians(180.0)
    min_lat = np.radians(-90.0)
    max_lat = np.radians( 90.0)

    scalar_nlons = 300
    scalar_nlats = scalar_nlons*(max_lat-min_lat)/(max_lon-min_lon)

    satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = satstress.Satellite(open(satfile,'r'))
    Diurnal = satstress.Diurnal(europa)
    diurnal_stresscalc = satstress.StressCalc([Diurnal,])

    the_fig = figure(figsize=(16,12))
    map_ax  = the_fig.add_subplot(1,1,1)
    map_ax.set_title(r'Diurnal Tidal Stresses ($nt=%.3g^\circ$)' % (orbital_pos,))

    grid_basemap = Basemap(llcrnrlon = np.degrees(min_lon),\
                             llcrnrlat = np.degrees(min_lat),\
                             urcrnrlon = np.degrees(max_lon),\
                             urcrnrlat = np.degrees(max_lat),\
                                    ax = map_ax)

    field_data = stressplot.scalar_grid(stresscalc=diurnal_stresscalc, nlons=scalar_nlons,\
                                        nlats=scalar_nlats, min_lon=min_lon,\
                                        max_lon=max_lon, min_lat=min_lat,\
                                        max_lat=max_lat, field=scalar_field,\
                                        time_t=(orbital_pos/360.0)*europa.orbit_period(),\
                                        cmap=cmap, basemap_ax=grid_basemap)
    vector_nlons = 13
    vector_nlats = 13

    vector_lons  = np.linspace(min_lon, max_lon, vector_nlons)
    vector_lats  = np.linspace(min_lat, max_lat, vector_nlats)
    vector_mesh_lons, vector_mesh_lats = np.meshgrid(vector_lons, vector_lats)
    vector_mesh_lons = np.ravel(vector_mesh_lons)
    vector_mesh_lats = np.ravel(vector_mesh_lats)

    if scale is None:
        diurnal_scale = max(np.fabs(diurnal_stresscalc.mean_global_stressmag()))*50
    else:
        diurnal_scale = scale

    stressplot.vector_points(stresscalc=diurnal_stresscalc, lons=vector_mesh_lons, lats=vector_mesh_lats,\
                             basemap_ax=grid_basemap, arrow_width=0.002,\
                             plot_tens=plot_tens, plot_comp=plot_comp,\
                             time_t=(orbital_pos/360.0)*europa.orbit_period(),\
                             plot_lesser=plot_lesser, plot_greater=plot_greater, scale=diurnal_scale)

    grid_basemap.drawmeridians(np.degrees(np.linspace(min_lon, max_lon, 13)), labels=[1,0,0,1], linewidth=0.1, color='gray')
    grid_basemap.drawparallels(np.degrees(np.linspace(min_lat, max_lat, 13)), labels=[1,0,0,1], linewidth=0.1, color='gray')
    grid_basemap.drawmapboundary()

    # Need some kind of scale bar for the scalar field in the background:
    cb_ax,kw = colorbar.make_axes(map_ax, orientation='vertical')
    colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=colors.Normalize(vmin=np.min(field_data),vmax=np.max(field_data)), orientation='vertical', format='%.2g')
    cb_ax.set_ylabel(scalar_field)
    plt.draw()

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,'StressFieldDiurnal_nt%.3g.%s' % (orbital_pos,save_fmt)))

#}}} end StressFieldNSR

def LinLengthDist(lins, label="", outfile=None): #{{{
    """
    Plot the overall distribution of feature lengths, color coding retained and
    rejected features differently.

    """
    # Need to generate two arrays of values pertaining to the combined set of
    # good and bad lineaments.
    #  - one of floats: lineament lengths
    #  - one of booleans: True/False depending on whether it was kept/rejected
    # Then I need to sort them both in order of lineament length.

    # The plot I want to draw is a curve, the Y-value of which is lineament
    # length, and the x-value is just N, the number of the lineament, ordered by length.
    # For each point in the curve, if the lineament was kept, it should be black beneath
    # and if it was rejected, it should be white.

    sat_radius_km = lins[0].stresscalc.stresses[0].satellite.radius()/1000

    # generate a list containing the lineament lengths in km:
    linlengths = [ l.length*sat_radius_km for l in lins ]

    # and a list of boolean values saying whether or not it passed:
    best_fit = [ l.best_fit()[0] for l in lins ]

    len_fits = array([ (linlengths[n],best_fit[n]) for n in range(len(lins)) ], dtype=[('length',float),('best_fit',float)])
    len_fits.sort(order='length')

    the_fig = figure(figsize=(9,3))
    plot_ax = the_fig.add_subplot(1,1,1)
    plot_ax.plot(len_fits['length'], color='black', linewidth=2)
    plot_ax.fill_between(np.arange(len(lins)), len_fits['length'], color='gray')
    plot_ax.set_title(label + ' Length Distribution')
    plot_ax.set_xlabel('N')
    plot_ax.set_ylabel('lineament length [km]')
    plot_ax.set_xlim(0,len(lins)-1)
    plot_ax.grid()

    if save_fmt is not None and outfile is not None:
        the_fig.savefig(os.path.join(figdir,outfile+'.'+save_fmt))
#}}}

def ActHist_ByLength(lins, dbar_max=0.125, norm_by_all=True): #{{{
    sat_radius_km = lins[0].stresscalc.stresses[0].satellite.radius()/1000
    labels_by_length = []
    lins_by_length = []
    lengths = [0,150,300,500,1000,2500]
    numlens = len(lengths)-1
    for N in range(numlens):
        lins_by_length.append([])
        labels_by_length.append("%d km < L < %d km" % (lengths[N], lengths[N+1]) )
        for lin in lins:
            lin_length_km = lin.length*sat_radius_km
            if (lin_length_km > lengths[N]) and (lin_length_km < lengths[N+1]):
                lins_by_length[N].append(lin)
                continue

    activity_history(lins_by_length, dbar_max=dbar_max, labels=labels_by_length, norm_by_all=norm_by_all, outfile='ActHist_ByLength')
#}}}

def ActHist_BySin(lins, dbar_max=0.125, norm_by_all=True): #{{{
    labels_by_sin = []
    lins_by_sin = []
    sins = [1.0, 1.02, 1.04, 1.06, 1.08, 1.10, 1.20]
    numsins = len(sins)-1
    for N in range(numsins):
        lins_by_sin.append([])
        labels_by_sin.append("%.3g < S < %.3g" % (sins[N], sins[N+1]) )
        for lin in lins:
            S = lin.sinuosity()
            if (S > sins[N]) and (S < sins[N+1]):
                lins_by_sin[N].append(lin)
                continue

    activity_history(lins_by_sin, dbar_max=dbar_max, labels=labels_by_sin, norm_by_all=norm_by_all, outfile='ActHist_BySin')
#}}}

def ActHist_ByDbar(maplins, norm_by_all=False): #{{{
    """
    Show that the time variability is not sensitive to the particular value
    chosen for dbar_max.

    """
    the_fig = figure(figsize=(12,8))

    dbars = (0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25)
    colors = cm.jet(linspace(0,1,len(dbars)))

    for dbar_max,N in zip(dbars,range(len(dbars))):
        activity_history([maplins,], labels=[r'$\bar{D}_{max}=%.3g$'%(dbar_max,),], norm_by_all=norm_by_all, the_fig=the_fig, dbar_max=dbar_max, colors=[colors[N],], outfile="ActHist_ByDbar")
#}}}

def ActHist_MapStats(maplins, nbins=20, N=100): #{{{
    """
    Do a Monte Carlo subsampling of our lineament map and compare the activity
    history of the true map with the synthetic ones, to see how self-consistent
    it is.

    """

    # pseudo-maps in transparent gray
    maps = [ linresample_byN(maplins) for i in range(N) ]
    alphas = [10./N,]*N
    colors = ['gray',]*N
    labels = [None,]*(N-1)
    labels.append('Synthetic Subsamples (N=%d)' % (N,))

    # map minus E15 mapping swath in solid blue, as an example of geographic exclusion
    labels.append('Map minus E15 Swath')
    noE15lins = [ lin for lin in maplins if (degrees(mod(array(lin.lons),2*pi)) > 290).all() or (degrees(mod(array(lin.lons),2*pi)) < 240).all() ]
    maps.append(noE15lins)
    colors.append('blue')
    alphas.append(1.0)

    # Real map in solid black
    labels.append('Mapped Lineaments')
    maps.append(maplins)
    colors.append('black')
    alphas.append(1.0)

    the_fig = figure(figsize=(12,8))
    activity_history(maps, labels=labels, colors=colors, alphas=alphas, the_fig=the_fig, verbose=False, titlestr="Natural Variability in the Mapped Dataset", outfile="ActHist_MapStats")
#}}}

def ActHist_SpinCycleStats(maplins, N=100): #{{{
    """
    Do a Monte Carlo subsampling of our lineament map, randomizing the
    longitudes of the subsampled data, and then compare the activity history of
    the true map with the synthetic ones, to see how far outside of the normal
    range our data is.

    """

    maps = [ make_crazy(linresample_byN(maplins),tpw=False) for i in range(N) ]
    maps.append(maplins)

    # pseudo-maps in gray
    colors = ['gray',]*N
    # real map in black
    colors.append('black')

    # pseudo-maps nearly transparent
    alphas = [15./N,]*N
    # real map opaque
    alphas.append(1.0)

    labels = [None,]*(N-1)
    labels.append('Spun Subsamples (N=%d)' % (N,))
    labels.append('Mapped Lineaments')

    the_fig = figure(figsize=(12,8))
    activity_history(maps, labels=labels, colors=colors, alphas=alphas, the_fig=the_fig, verbose=False, outfile="ActHist_SpinCycleStats")
#}}}

def ActHist_GreatCircleStats(gclins, maplins, N=100): #{{{
    """
    Re-sample the synthetic great circles dataset so we can see how much spread
    there is in the resulting activity histories.

    """
    maps = [ make_crazy(linresample_byN(gclins),tpw=False) for i in range(N) ]
    maps.append(maplins)
    colors = ['gray',]*N
    colors.append('black')
    alphas = [15./N,]*N
    alphas.append(1.0)
    labels = [None,]*(N-1)
    labels.append('Re-sampled Great Circle Segments (N=%d)' % (N,))
    labels.append('Mapped Lineaments')

    the_fig = figure(figsize=(12,8))
    activity_history(maps, labels=labels, colors=colors, alphas=alphas, the_fig=the_fig, verbose=False, outfile="ActHist_GreatCircleStats")

#}}} end ActHist_GreatCricleStats()

def ActHist_PeakStats(maplins, N=100, scale=np.radians(15)): #{{{
    """
    Create a large number of synthetic datasets composed partly of a
    longitudinally randomized sampling of maplins, and partly from a de-rotated
    set of maplins, and see how their activity histories compare to that of the
    mapped features in their real locations.

    """

    map_acthist = calc_acthist(maplins)
    mu = int(where(np.abs(map_acthist-np.max(map_acthist)) < 1e-6)[0])*(np.pi/len(map_acthist))
    maps = []
    mu, maps = synthpeak(maplins, N=N, peak_frac=1.0, scale=scale)
    maps.append(maplins)

    # pseudo-maps in gray
    colors = ['gray',]*N
    # real map in black
    colors.append('black')

    # pseudo-maps nearly transparent
    alphas = [10./N,]*N
    # real map opaque
    alphas.append(1.0)

    labels = [None,]*(N-1)
    labels.append('Synthetic Peak (N=%d)' % (N,))
    labels.append('Mapped Lineaments')

    the_fig = figure(figsize=(12,8))
    titlestr='Peak subsampled from derotated map $\mu=%g^\circ$, $\sigma=%g^\circ$' % (degrees(mu),degrees(scale))
    activity_history(maps, labels=labels, colors=colors, alphas=alphas, the_fig=the_fig, verbose=False, titlestr=titlestr, outfile='ActHist_PeakStats')
#}}}

def ActHist_Synthesized(maplins, N=100, peak_frac=0.4, scale=np.radians(15)): #{{{
    """
    Create a large number of synthetic datasets composed partly of a
    longitudinally randomized sampling of maplins, and partly from a de-rotated
    set of maplins, and see how their activity histories compare to that of the
    mapped features in their real locations.

    """

    mu, maps = synthpeak(maplins, N=N, peak_frac=peak_frac, scale=scale)
    maps.append(maplins)

    # pseudo-maps in gray
    colors = ['gray',]*N
    # real map in black
    colors.append('black')

    # pseudo-maps nearly transparent
    alphas = [10./N,]*N
    # real map opaque
    alphas.append(1.0)

    labels = [None,]*(N-1)
    labels.append('%d%% Rand + %d%% Peak (N=%d)' % (int(100*(1-peak_frac)),int(100*peak_frac),N))
    labels.append('Mapped Lineaments')

    the_fig = figure(figsize=(12,8))
    titlestr='Mix of random and peaked features w/ $\mu=%g^\circ$, $\sigma=%g^\circ$' % (degrees(mu),degrees(scale))
    activity_history(maps, labels=labels, colors=colors, alphas=alphas, the_fig=the_fig, verbose=False, outfile='ActHist_Synthesized' )
#}}}

def FitMap(lins, titlestr="Features colored by fit", lin_cm=cm.jet, nbins=18, stresscentric=False, outfile=None, dbar_max=0.125, showbad=True, derotate=False): #{{{
    """
    Creates a global map of the lineaments, color coding them by what amount of
    backrotation (b) minimizes delta_rms(b) when compared to the NSR stress
    field.  Lineament width indicates the value of min(delta_rms(b)), with
    wider lineaments agreeing better with NSR stresses.

    Those lineaments which did not meet the criteria for inclusion in the
    analysis are plotted thin and black.

    """

    # We can't both derotate and show the bad features, since bad features
    # don't necessarily have a best fit backrotation
    if derotate is True:
        showbad = False

    lin_cm = colors.ListedColormap(lin_cm(linspace(0,1,nbins)))
    sat_radius_km = lins[0].stresscalc.stresses[0].satellite.radius()/1000

    the_fig = figure(figsize=(12,8))
    if stresscentric:
        llcrnrlon=0
        llcrnrlat=0
        urcrnrlon=180
        urcrnrlat=90
        lat_mirror=True
        gridspace=15
    else:
        llcrnrlon=0
        llcrnrlat=-90
        urcrnrlon=360
        urcrnrlat=90
        lat_mirror=False
        gridspace=30

    lon_cyc=abs(radians(llcrnrlon-urcrnrlon))
    linfitmap = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat)
    linfitmap.drawmeridians(range(llcrnrlon,urcrnrlon+1,gridspace), labels=[1,0,0,1])
    linfitmap.drawparallels(range(llcrnrlat,urcrnrlat+1,gridspace), labels=[1,0,0,1])
    map_ax = the_fig.axes[0]

    for lin in lins:
        if len(lin.bs) > 0:
            backrot = 0
            # load some information about the lineament's best fit:
            best_fit, best_b = lin.best_fit(dbar_max=dbar_max)

            if best_fit > 0.0:
                if derotate:
                    backrot = best_b
                # Map the color of the lineament to its best_b
                lin_color = lin_cm(int((lin_cm.N)*(0.5+(best_b/pi))))
                # use line width to indicate goodness of best_fit
                lin_width = 1.0 + 5.0*best_fit

            elif showbad:
                lin_width = 1.0
                lin_color = 'black'

            else:
                continue

            if backrot == 0:
                lin2plot = lin
            else:
                lin2plot = lineament.Lineament(lons=lin.lons+backrot, lats=lin.lats)

            newline, linfitmap = lineament.plotlinmap([lin2plot,], map=linfitmap, linewidth=lin_width, color=lin_color, lon_cyc=lon_cyc, lat_mirror=lat_mirror)

    map_ax.set_title(titlestr)
    cb_ax,kw = colorbar.make_axes(map_ax, orientation="horizontal", pad=0.05, shrink=0.5)
    colorbar.ColorbarBase(cb_ax, cmap=lin_cm, orientation="horizontal", norm=colors.BoundaryNorm(linspace(-90,90,nbins+1),nbins), format=r'%.0f$^\circ$')
    # Fix up the colorbar a bit:
    #cb_ax.invert_xaxis()
    cb_ax.set_xlabel("longitudinal shift b")

    the_fig.show()
    if outfile is not None and save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,outfile+'.'+save_fmt))

# end FitMap }}}

def ActHist_RMDProb(maplins, tpwlins): #{{{
    """
    Create a figure showing the distribution of activity history non-uniformity
    for three different classes of maps:

      1. mapped features, under TPW.
      2. resampled map with randomized longitudes (spun).
      3. resampled map completely randomized and re-fit (jackstraws).

    Additionally, indicate the non-uniformity of several individual maps:
      1. Features as mapped.
      2. Mapped features derotated (maximum possible non-uniformity).
      3. The putative TPW orientation (PNP=80E10N).

    """
    import pickle

    spin_acthists = pickle.load(open('output/spin_map_acthists'))
    spin_nb = len(spin_acthists[0])
    spin_RMDs = [ np.fabs(np.tile(ah,len(ah)) - np.repeat(ah,len(ah))).mean()/ah.mean() for ah in spin_acthists ]

    jack_acthists = pickle.load(open(os.path.join(outdir,'jack_map_acthists')))
    jack_nb = len(jack_acthists[0])
    jack_RMDs = [ np.fabs(np.tile(ah,len(ah)) - np.repeat(ah,len(ah))).mean()/ah.mean() for ah in jack_acthists ]

    pnp_lons, pnp_lats, pnp_RMDs = load_tpw_RMDs(tpwdir='tpw_polesearch', fast=True)

    # In order for this to be a fair comparison, we have to sample the activity
    # histories at the same resolution.  Note, we're also assuming here that nb
    # (the number of 'b' values) used in the tpw_polesearch is the same as
    # jack_nb... which is true if we've used nsrhist.fromscratch() to do the
    # calculation, and which makes sense, because re-fitting takes forever...
    assert jack_nb == spin_nb

    if len(maplins[0].bs) != spin_nb:
        downsample_factor = len(maplins[0].bs)/spin_nb
        print("downsampling maplins fits by a factor of %g" % (downsample_factor,) )
        lowres_maplins = lineament.update_lins(maplins)
        for lin in lowres_maplins:
            lin.bs = lin.bs[::downsample_factor]
            lin.nsrdbars = lin.nsrdbars[::downsample_factor]
            lin.nsrstresswts = lin.nsrstresswts[::downsample_factor]
    else:
        lowres_maplins = maplins
    map_RMD = acthist_RMD(lowres_maplins)

    if len(tpwlins[0].bs) != spin_nb:
        downsample_factor = len(tpwlins[0].bs)/spin_nb
        print("downsampling tpwlins fits by a factor of %g" % (downsample_factor,) )
        lowres_tpwlins = lineament.update_lins(tpwlins)
        for lin in lowres_tpwlins:
            lin.bs = lin.bs[::downsample_factor]
            lin.nsrdbars = lin.nsrdbars[::downsample_factor]
            lin.nsrstresswts = lin.nsrstresswts[::downsample_factor]
    else:
        lowres_tpwlins = tpwlins
    tpw80E10N_RMD = acthist_RMD(lowres_tpwlins)

    bestfits = [ lin.best_fit()[1] for lin in lowres_maplins ]
    derotated = [ lin.lonshift(b) for lin,b in zip(lowres_maplins, bestfits) ]
    derotated_RMD = acthist_RMD(derotated)
    print("greatest possible RMD w/ map: %g" % (derotated_RMD,))

    the_fig = figure(figsize=(12,8))
    hist_ax = the_fig.add_subplot(1,1,1)
    spin_counts, spin_bins, spin_patches = hist_ax.hist(spin_RMDs, bins=96, range=(0.0,0.6), color='red',   lw=2, normed=True, histtype='step', label='spin (N=%d)' % (len(spin_RMDs),))
    jack_counts, jack_bins, jack_patches = hist_ax.hist(jack_RMDs, bins=48,  range=(0.0,0.6), color='green', lw=2, normed=True, histtype='step', label='jack (N=%d)' % (len(jack_RMDs),))
    pnp_counts,  pnp_bins,  pnp_patches  = hist_ax.hist(pnp_RMDs,  bins=48,  range=(0.0,0.6), color='blue',  lw=2, normed=True, histtype='step', label='TPW (N=%d)'  % (len(pnp_RMDs),))
    hist_ax.set_xlim(0,0.6)
    hist_ax.set_xticks(linspace(0,0.6,13))

    # add a vertical line to show where the actual mapped features fall:
    hist_ax.axvline(x=median(spin_RMDs), linewidth=3, linestyle=':', color='red')
    hist_ax.annotate('median: %.3g' % (median(spin_RMDs),), xy=(median(spin_RMDs)-0.003, max(spin_counts)/5.0), rotation='vertical', ha='right', color='red')

    hist_ax.axvline(x=median(jack_RMDs), linewidth=3, linestyle=':', color='green')
    hist_ax.annotate('median: %.3g' % (median(jack_RMDs),), xy=(median(jack_RMDs)+0.003, max(jack_counts)/3.0), rotation=270, ha='left', color='green')

    hist_ax.axvline(x=tpw80E10N_RMD, linewidth=3, linestyle='--', color='gray')
    hist_ax.annotate('PNP 80E10N: %.3g' % (tpw80E10N_RMD,), xy=(tpw80E10N_RMD-0.003, 0.7*max(spin_counts)), rotation='vertical', ha='right', color='gray')

    hist_ax.axvline(x=map_RMD, linewidth=3, linestyle='--', color='black')
    hist_ax.annotate('mapped: %.3g' % (map_RMD,), xy=(map_RMD-0.003, 0.5*max(spin_counts)), rotation='vertical', ha='right', color='black')

    hist_ax.set_title(r'Activity history RMD distributions for several lineament populations')
    hist_ax.set_xlabel('RMD(H(b))')
    hist_ax.set_ylabel('N (normalized)')
    hist_ax.legend(loc='upper right')
    hist_ax.grid()

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,'ActHist_RMDProb.'+save_fmt))

#}}}

def ActHist_RMDMap(tpwdir="tpw_polesearch", maplins=None, ncols=50, cmap='jet', fast=True, npole=False): #{{{
    """
    Reads in a set of paleopole longitudes and latitudes, and the corresponding
    lineaments and their fits from tpwdir, and creates a 2-D plot showing how
    the relative mean difference (RMD) of the activity history varies with
    paleopole location.

    If fast is True, read in the values from a pickled cache, named
    'tpwdir'.pkl

    if maplins is not None, draw a contour on the map corresponding to the RMD
    of the inferred activity history of maplins.  Be careful that you make sure
    maplins has had its fits calculated at the same resolution as the TPW fits
    were done at (otherwise it's not a fair comparison).

    """

    pnp_lons, pnp_lats, acthist_RMDs = load_tpw_RMDs(tpwdir=tpwdir, fast=fast)

    if maplins is not None:
        map_acthist_RMD = acthist_RMD(maplins)
    else:
        map_acthist_RMD = 0.0

    levels = np.linspace(acthist_RMDs.min(), acthist_RMDs.max(), ncols)

    # replicate data to fill the whole surface (w/ symmetry)
    pnp_lons = np.concatenate([pnp_lons, pnp_lons+np.pi])
    pnp_lons = np.mod(pnp_lons,2*np.pi)-np.pi
    pnp_lats = np.concatenate([pnp_lats, -pnp_lats])
    acthist_RMDs = np.tile(acthist_RMDs,2)

    # Create a three-panel plot, one for the equatorial to mid-latitudes in a
    # simple cylindrical projection, and one each for the poles, in polar
    # stereographic underneath the other one.
    if npole is True:
        the_fig = figure(figsize=(9,12))
        np_ax = the_fig.add_subplot(1,1,1)
        cb_ax,kw = colorbar.make_axes(np_ax, fraction=0.1, aspect=20, orientation='horizontal')
    else:
        the_fig = figure(figsize=(12,12))
        lowlat_ax = the_fig.add_subplot(2,1,1)
        np_ax = the_fig.add_subplot(2,2,3)
        sp_ax = the_fig.add_subplot(2,2,4)
        cb_ax,kw = colorbar.make_axes(lowlat_ax, fraction=0.075, aspect=40, orientation='horizontal')

    # Ideally, all this interpolation would be done in spherical coordinates,
    # but that's not built-in.  Getting it right in the poles is important, and
    # can be done well by doing the interpolation in the map-projected space.
    # For the low latitude regions, we can get away with just setting up a
    # pseudo-periodic environment by copying the data a few times on either
    # side of the region we actually care about.

    # North Pole:
    np_map = Basemap(projection='npstere',boundinglat=40,lon_0=90, ax=np_ax)
    np_x_data, np_y_data = np_map(degrees(pnp_lons), degrees(pnp_lats))
    np_x_grid = linspace(np_map.xmin,np_map.xmax,100)
    np_y_grid = linspace(np_map.ymin,np_map.ymax,100)
    np_pnp_amps = griddata(np_x_data, np_y_data, acthist_RMDs, np_x_grid, np_y_grid)
    np_x_gridmesh, np_y_gridmesh = meshgrid(np_x_grid,np_y_grid)
    np_map.contourf(np_x_gridmesh, np_y_gridmesh, np_pnp_amps, levels, cmap=cm.get_cmap(cmap, ncols-1), linewidth=0)
    np_map.contour(np_x_gridmesh, np_y_gridmesh, np_pnp_amps, [map_acthist_RMD,], colors='black', linewidths=[2,])
    #np_map.scatter(np_x_data, np_y_data, s=1, alpha=0.5, color='black')
    np_map.drawparallels(arange(-60.,61.,30.), labels=[1,0,0,1], dashes=[2,5], linewidth=0.5)
    np_map.drawmeridians(arange(-180.,181.,30.), labels=[1,0,0,0], dashes=[2,5], linewidth=0.5)

    if npole is True:
        np_ax.set_title("Non-uniformity (RMD) of activity history H(b)")
        np_map.drawparallels(arange(-60.,61.,30.), labels=[0,1,0,0], dashes=[2,5], linewidth=0.5)
        np_map.drawmeridians(arange(-180.,181.,30.), labels=[0,1,0,0], dashes=[2,5], linewidth=0.5)

    if npole is False:
    # change some labels so things aren't cluttered...
        np_ax.set_title("North Pole")
    # South Pole:
        sp_map = Basemap(projection='spstere',boundinglat=-40,lon_0=270, ax=sp_ax)
        sp_x_data, sp_y_data = sp_map(degrees(pnp_lons), degrees(pnp_lats))
        sp_x_grid = linspace(sp_map.xmin,sp_map.xmax,100)
        sp_y_grid = linspace(sp_map.ymin,sp_map.ymax,100)
        sp_pnp_amps = griddata(sp_x_data, sp_y_data, acthist_RMDs, sp_x_grid, sp_y_grid)
        sp_x_gridmesh, sp_y_gridmesh = meshgrid(sp_x_grid,sp_y_grid)
        sp_map.contourf(sp_x_gridmesh, sp_y_gridmesh, sp_pnp_amps, levels, cmap=cm.get_cmap(cmap, ncols-1))
        sp_map.contour(sp_x_gridmesh, sp_y_gridmesh, sp_pnp_amps, [map_acthist_RMD,], colors='black', linewidths=[2,])
        #sp_map.scatter(sp_x_data, sp_y_data, s=2, alpha=0.5, color='black', linewidth=0)
        sp_map.drawparallels(arange(-60.,61.,30.), labels=[0,1,0,1], dashes=[2,5], linewidth=0.5)
        sp_map.drawmeridians(arange(-180.,181.,30.), labels=[0,1,0,0], dashes=[2,5], linewidth=0.5)
        sp_ax.set_title("South Pole")

    # Low Latitudes:
        lowlat_map = Basemap(llcrnrlon=-180,llcrnrlat=-60,urcrnrlon=180,urcrnrlat=60, ax=lowlat_ax)
        lowlat_x_data = pnp_lons
        lowlat_y_data = pnp_lats
        lowlat_acthist_RMDs = acthist_RMDs
        lowlat_x_grid = linspace(-pi,pi,360)
        lowlat_y_grid = linspace(-pi/2,pi/2,180)
        lowlat_pnp_amps = griddata(lowlat_x_data, lowlat_y_data, lowlat_acthist_RMDs, lowlat_x_grid, lowlat_y_grid)
        lowlat_x_gridmesh, lowlat_y_gridmesh = meshgrid(lowlat_x_grid,lowlat_y_grid)
        lowlat_map.contourf(degrees(lowlat_x_gridmesh), degrees(lowlat_y_gridmesh), lowlat_pnp_amps, levels, cmap=cm.get_cmap(cmap, ncols-1))
        lowlat_map.contour(degrees(lowlat_x_gridmesh),  degrees(lowlat_y_gridmesh), lowlat_pnp_amps, [map_acthist_RMD,], colors='black', linewidths=[2,])
        #lowlat_map.scatter(degrees(lowlat_x_data), degrees(lowlat_y_data), s=2, alpha=0.5, color='black', linewidth=0)
        lowlat_map.drawparallels(arange(-60.,61.,30.), labels=[1,1,0,1], dashes=[2,5], linewidth=0.5)
        lowlat_map.drawmeridians(arange(-180.,181.,30.), labels=[1,1,0,1], dashes=[2,5], linewidth=0.5)
        lowlat_ax.set_title("Equatorial and Mid-Latitudes")

    colorbar.ColorbarBase(cb_ax, cmap=cm.get_cmap(cmap, ncols-1), norm=colors.Normalize(vmin=levels[0],vmax=levels[-1]), orientation='horizontal', format='%.2f')
    cb_ax.set_xlabel("RMD(H(b))")
    plt.draw()

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,'ActHist_RMDMap.'+save_fmt))

#}}}

def DbarLengthCorr(lins): #{{{
    """
    Several scatter plots showing the correlation (or lack thereof) between
    quality of best fit and lineament length, in equatorial and non-equatorial
    regions.

    """

    the_fig = figure(figsize=(12,8))
    ax1 = the_fig.add_subplot(1,1,1)

    hilatlins = [ lin for lin in lins if degrees(min(fabs(lin.lats))) > 30 ]
    hilat_best_fits = array([ lin.nsrdbars.min() for lin in hilatlins ])
    hilat_lengths = array([ lin.length for lin in hilatlins ])*1561
    hilat_r2 = corrcoef(hilat_lengths, hilat_best_fits)[0,1]**2
    hilat_symbs = scatter(hilat_lengths, hilat_best_fits, c='blue', marker='s', label=r'$|\theta_{min}|>30^\circ R^2=%.3g$' % ( hilat_r2, ) )

    lolatlins = [ lin for lin in lins if degrees(min(fabs(lin.lats))) <= 30 ]
    lolat_best_fits = array([ lin.nsrdbars.min() for lin in lolatlins ])
    lolat_lengths = array([ lin.length for lin in lolatlins ])*1561
    lolat_r2 = corrcoef(lolat_lengths,lolat_best_fits)[0,1]**2
    lolat_symbs = scatter(lolat_lengths, lolat_best_fits, c='red', marker='^', label=r'$\|\theta_{min}|\leq30^\circ R^2=%.3g$' % ( lolat_r2, ) )

    ax1.legend(loc='lower right')
    ax1.set_xlabel('lineament length [km]')
    ax1.set_ylabel(r'$\min(\bar{D}(b))$')
    ax1.set_ylim(0,0.3)
    ax1.set_xlim( min((hilat_lengths.min(),lolat_lengths.min())),\
                  max((hilat_lengths.max(),lolat_lengths.max())) )
    ax1.grid(True)
    ax1.set_title('Effects of length and latitude on fit to NSR for great circle segments, N=%d'\
                   % (len(hilatlins)+len(lolatlins)) )

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,'DbarLengthCorr.'+save_fmt))

# end DbarLengthCorr}}}

def DbarSinuosityCorr(lins): #{{{
    """
    Scatter plot showing the correlation (or lack thereof) between quality of
    best fit and lineament sinuosity in equatorial and non-equatorial regions.

    """

    the_fig = figure(figsize=(12,8))
    ax1 = the_fig.add_subplot(1,1,1)

    hilatlins = [ lin for lin in lins if degrees(min(fabs(lin.lats))) > 30 ]
    hilat_best_fits = array([ lin.nsrdbars.min() for lin in hilatlins ])
    hilat_sins = array([ lin.sinuosity() for lin in hilatlins ])
    hilat_r2 = corrcoef(hilat_sins,hilat_best_fits)[0,1]**2
    hilat_symbs = scatter(hilat_sins, hilat_best_fits, c='blue', marker='s', label=r'$|\theta_{min}|>30^\circ R^2=%.3g$' % ( hilat_r2, ) )

    lolatlins = [ lin for lin in lins if degrees(min(fabs(lin.lats))) <= 30 ]
    lolat_best_fits= array([ lin.nsrdbars.min() for lin in lolatlins ])
    lolat_sins = array([ lin.sinuosity() for lin in lolatlins ])
    lolat_r2 = corrcoef(lolat_sins,lolat_best_fits)[0,1]**2
    lolat_symbs = scatter(lolat_sins, lolat_best_fits, c='red', marker='^', label=r'$\|\theta_{min}|\leq30^\circ R^2=%.3g$' % ( lolat_r2, ) )

    ax1.legend(loc='lower right')
    ax1.set_xlabel('sinuosity')
    ax1.set_ylabel(r'$\min(\bar{D}(b))$')
    ax1.set_xlim(1,1.1)
    ax1.set_ylim(0,0.3)
    ax1.grid(True)
    ax1.set_title('Effects of sinuosity and latitude on fit to NSR for all mapped features, N=%d'\
                   % (len(hilatlins)+len(lolatlins)) )

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,'DbarSinCorr.'+save_fmt))

# end DbarSinuosityCorr }}}

def LinDensityMap(lins, maxdist=500, N=0, label="", cmap=cm.jet): #{{{
    """
    Calculate an interpolated grid showing the density of lineaments on the
    surface, as reflected by the sum of the lengths of the lineament segments
    within the given distance maxdist, of N sample points on the surface of the
    satellite.

    Compare that lineament density to the resolution of coverage we have from
    the USGS mosaic that was used to make the lineament map, determining to
    what degree the resolution of coverage determines the density of mapped
    lineaments.

    """
    import matplotlib.colors as colors
    from osgeo import gdal
    import re

    # The number of points we need to sample the distribution at scales
    # proportional to the surface area that each point is sampling.  This
    # number works well for Europa anyway.
    if N==0:
        N = 5.0e8/(maxdist**2)

    #randlons, randlats = lineament.random_lonlatpoints(N)
    randlons, randlats = fibonacci_sphere(N)
    randlons = mod(randlons,2*pi)
    reglats = linspace(-90,90,180)
    reglons = linspace(0,360,360)

    # adding these corner points makes sure that the griddata goes to the edges of the map
    edge_N = 10
    toplats = array(pi/1.99).repeat(edge_N)
    toplons = linspace(-0.01,2.01*pi,edge_N)
    bottomlats = -toplats
    bottomlons = toplons
    westlats = linspace(-pi/1.98,pi/1.98,edge_N)
    westlons = array(-0.02).repeat(edge_N)
    eastlats = westlats
    eastlons = array(2.01*pi).repeat(edge_N)
    randlons = hstack([randlons, toplons, bottomlons, eastlons, westlons])
    randlats = hstack([randlats, toplats, bottomlats, eastlats, westlats])

    seglons = array([])
    seglats = array([])
    seglens = array([])
    for lin in lins:
        newlons, newlats = lin.seg_midpoints()
        newlens = lin.seg_lengths()
        seglons = concatenate([seglons, newlons])
        seglats = concatenate([seglats, newlats])
        seglens = concatenate([seglens, newlens])

    # For each point (lon,lat) defined by randlons, randlats, calculate the sum
    # of the lengths of the segments closer than dist radians of arc away:
    nsegs = len(seglens)
    lensums = array([])

    print("Calculating lineament density map with d=%d km and N=%d" % (maxdist, N) )
    for lon,lat in zip(randlons, randlats):
        lon_arr = array(lon)
        lat_arr = array(lat)
        newsum = sum(where(lineament.spherical_distance(lon_arr.repeat(nsegs),lat_arr.repeat(nsegs),seglons,seglats) < maxdist/1561.0, seglens, 0.0))
        lensums = hstack([lensums,newsum])

    # convert these values of radians per footprint, into m/km^2
    lindensity = lensums*1561*1000/(pi*maxdist**2)

    lindensity_grid = griddata(degrees(randlons), degrees(randlats), lindensity, reglons, reglats)

    # get rid of the out-of-bounds points we added so that the interpolated grid would stretch to the edges of the map
    randlons, randlats, lindensity = randlons[:N], randlats[:N], lindensity[:N]

    lindens_fig = figure(figsize=(12,8))
    lindens_ax = lindens_fig.add_subplot(1,1,1)
    lindens_ax.contourf(reglons, reglats, lindensity_grid, 64, cmap=cmap, linewidth=0)
    lindens_ax.scatter(degrees(randlons), degrees(randlats), marker='o', color='white', s=2, edgecolor='white', alpha=0.5, linewidth=0)
    lindens_ax.scatter(mod(degrees(seglons),360), degrees(seglats), marker='o', color='black', s=200*seglens, edgecolor='black', alpha=0.375, linewidth=0)

    lindens_ax.set_xlim(0,360)
    lindens_ax.set_xticks(linspace(0,360,13))
    lindens_ax.set_ylim(-90,90)
    lindens_ax.set_yticks(linspace(-90,90,7))
    lindens_ax.set_title(label+" N=%d, d=%g km" % (N,maxdist) )

    lindensnorm = colors.Normalize(vmin=0,vmax=lindensity.max())
    cb_ax,kw = colorbar.make_axes(lindens_ax, pad=0.05, orientation='horizontal', shrink=0.5)
    colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=lindensnorm, orientation='horizontal')
    cb_ax.set_xlabel(r"mapped lineament density [ m/km$^2$]")

    if save_fmt is not None:
       lindens_fig.savefig(os.path.join(figdir,'LinDensity_Grid.'+save_fmt))

    print("Reading in USGS mosaic raster data")
    # Read in the resolution/coverage map from USGS:
    # The data cube has 4 bands:
    # 1. resolution [km/px]
    # 2. emission angle [degrees]
    # 3. incidence angle [degrees]
    # 4. phase angle [degrees]
    EuropaGrid = gdal.Open("input/europa_coverage_simpcyl.cub")
    raster_numlons = EuropaGrid.RasterXSize
    raster_numlats = EuropaGrid.RasterYSize
    raster_lonidx = raster_numlons*(randlons/(2*pi))-1.0
    raster_lonidx = raster_lonidx.round().astype(int)
    raster_latidx = (raster_numlats*(pi/2-(randlats))/pi)-1.0
    raster_latidx = raster_latidx.round().astype(int)

    # changing resolution to [px/km]:
    resolution_raster = 1.0 / EuropaGrid.GetRasterBand(1).ReadAsArray()
    # mask the raster to remove extreme polar distortions, nodata values
    resolution_raster = ma.masked_outside(resolution_raster, 0.0, 6.0)

    emission_angle_raster = EuropaGrid.GetRasterBand(2).ReadAsArray()
    emission_angle_raster = ma.masked_where(emission_angle_raster < 0, emission_angle_raster)
    incidence_angle_raster = EuropaGrid.GetRasterBand(3).ReadAsArray()
    incidence_angle_raster = ma.masked_where(incidence_angle_raster < 0, incidence_angle_raster)
    phase_angle_raster    = EuropaGrid.GetRasterBand(4).ReadAsArray()
    phase_angle_raster    = ma.masked_where(phase_angle_raster < 0, phase_angle_raster)

    rasters   = [resolution_raster, emission_angle_raster, incidence_angle_raster, phase_angle_raster]
    rastnames = ['Resolution', 'Emission Angle', 'Incidence Angle', 'Phase Angle']
    rastunits = ['[km/px]', '[degrees]', '[degrees]', '[degrees]']
    rastfigs = [ figure(figsize=(9,13)), figure(figsize=(9,13)), figure(figsize=(9,13)), figure(figsize=(9,13)) ]
    rastoutfiles = [ 'LinDensity_Resolution', 'LinDensity_Emission', 'LinDensity_Incidence', 'LinDensity_Phase' ]

    for raster,rastname,rastunit,rastfig,outfile in zip(rasters, rastnames, rastunits, rastfigs, rastoutfiles):
        rast_ax = rastfig.add_subplot(2,1,1)
        rast_ax.imshow(raster, extent=(0,360,-90,90))
        rast_ax.scatter(mod(degrees(seglons),360), degrees(seglats), marker='o', color='black', s=200*seglens, edgecolor='black', alpha=0.375, linewidth=0)

        rast_ax.set_xlim(0,360)
        rast_ax.set_xticks(linspace(0,360,13))
        rast_ax.set_ylim(-90,90)
        rast_ax.set_yticks(linspace(-90,90,7))
        rast_ax.set_title("USGS Europa Mosaic: %s" % (rastname,) )
        rast_norm = colors.Normalize(vmin=0,vmax=raster.max())
        cb_ax,kw = colorbar.make_axes(rast_ax, pad=0.1, orientation='horizontal', shrink=0.5)
        colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=rast_norm, orientation='horizontal')
        cb_ax.set_xlabel(rastunit)

        # The array lindensity contains values of lineament density for the N
        # (lon,lat) points defined by (randlons, randlats).  I need a similar
        # set of values for the same (lon,lat) points, but corresponding to the
        # values stored within various rasters we're comparing to:
        randsamples = ma.masked_less(raster[raster_latidx, raster_lonidx],0)
        # Calculate the correlation between lindensity and the raster in question:
        lindensity_raster_corrcoef = ma.corrcoef(randsamples, lindensity)[0,1] 
        print(r"lindensity v. %s: $R^2$=%g" % (rastname, lindensity_raster_corrcoef**2) )

        # Make a scatter plot showing the correlation (or lack thereof):
        lin_rast_corr_ax = rastfig.add_subplot(2,1,2)
        lin_rast_corr_ax.scatter(randsamples, lindensity, s=10, linewidth=0, marker='o', color='black', alpha=0.375)
        lin_rast_corr_ax.set_xlabel("USGS Mosaic %s, %s" % (rastname, rastunit) )
        lin_rast_corr_ax.set_ylabel(r'Lineament density [m/km$^2$]')
        lin_rast_corr_ax.set_title(r'd=%g km, N=%d, $R^2$=%.4g' % (maxdist, N, lindensity_raster_corrcoef**2) )
        lin_rast_corr_ax.set_ylim(0,lindensity.max())
        lin_rast_corr_ax.set_xlim(0,randsamples.max())
        lin_rast_corr_ax.grid(True)

        if save_fmt is not None and outfile is not None:
            rastfig.savefig(os.path.join(figdir,outfile+'.'+save_fmt))
#}}} end LinDensityMap

def FitCurveExamples(lins, labels=[], dbar_max=0.125): #{{{
    # Create a full page plot, with the top half consisting of a map of the example lineaments,
    # with each one having a different color.  In the bottom half, show on individual subplots
    # the fitcurves for each of the features, color coded similarly.
    the_fig= figure(figsize=(9,12))

    colors = cm.jet(linspace(0,1,len(lins)))

    eg_ax1 = the_fig.add_axes((0.1,0.4625,0.4,0.1375))
    eg_ax2 = the_fig.add_axes((0.5,0.4625,0.4,0.1375))
    eg_ax3 = the_fig.add_axes((0.1,0.3250,0.4,0.1375))
    eg_ax4 = the_fig.add_axes((0.5,0.3250,0.4,0.1375))
    eg_ax5 = the_fig.add_axes((0.1,0.1875,0.4,0.1375))
    eg_ax6 = the_fig.add_axes((0.5,0.1875,0.4,0.1375))
    eg_ax7 = the_fig.add_axes((0.1,0.0500,0.4,0.1375))
    eg_ax8 = the_fig.add_axes((0.5,0.0500,0.4,0.1375))

    eg_axes = [eg_ax1,eg_ax2,eg_ax3,eg_ax4,eg_ax5,eg_ax6,eg_ax7,eg_ax8]

    # this is the map
    eg_map_ax = the_fig.add_axes((0.1,0.6,0.8,0.4))
    llcrnrlon=0
    llcrnrlat=0
    urcrnrlon=180
    urcrnrlat=90
    eg_map = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat, ax=eg_map_ax)
    eg_map.drawmeridians(range(llcrnrlon,urcrnrlon+1,15), labels=[1,1,0,1], linewidth=0.25, dashes=(2,5))
    eg_map.drawparallels(range(llcrnrlat,urcrnrlat+1,15), labels=[1,1,0,1], linewidth=0.25, dashes=(2,5))

    for lin,eg_ax,N in zip(lins,eg_axes,range(len(lins))):
        # Plot the lineament, and color it
        lineament.plotlinmap([lin,], map=eg_map, lon_cyc=radians(180), lat_mirror=True, color=colors[N], linewidth=2)
        fitcurve(lin, ax=eg_ax, color=colors[N], dbar_max=dbar_max)

    # clean up the massively multiple axes:
    ys_to_hide = [ the_fig.axes[N] for N in (1,3,5,7,9,11,13,15) ]
    [ ax.set_ylabel('') for ax in ys_to_hide ]
    [ setp(ax.get_yticklabels(),visible=False) for ax in ys_to_hide ]
    [ setp(ax.get_xticklabels(),visible=False) for ax in the_fig.axes[0:6] ]
    [ setp(ax.get_xticklabels(),visible=False) for ax in the_fig.axes[8:17] ]
    [ the_fig.axes[N].set_ylabel(r'$\bar{D}(b)$') for N in (0,2,4,6) ]
    [ the_fig.axes[N].set_ylabel(r'$H(b)$', rotation=270) for N in (10,12,14,16) ]

    eg_map_ax.set_title(r'Example Lineaments and Fit Probabilities')

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,'FitCurveExamples.'+save_fmt))
#}}}

def LinLatStatsCompare(maplins, synlins=None, bins=30, fold=False): #{{{

    the_fig = plt.figure(figsize=(12,6))
    hist_ax = the_fig.add_subplot(1,1,1)

    max_lat = np.pi/2
    if fold is True:
        min_lat = 0
    else:
        min_lat = -np.pi/2

    bin_edges   = linspace(min_lat,max_lat,bins+1)
    bin_centers = (bin_edges[1:]+bin_edges[:-1])/2.0
    print degrees(bin_centers)

    if fold is True:
        hist_ax.set_xlim(0,89)
        hist_ax.set_xticks(np.linspace(0,90,7))
        plot_lats = np.fabs(np.degrees(bin_centers))
    else:
        hist_ax.set_xlim(-89,89)
        hist_ax.set_xticks(np.linspace(-90,90,7))
        plot_lats = np.degrees(bin_centers)

    print("Getting uniform latitude distribution")
    gc_lons, gc_lats = lineament.random_lonlatpoints(1e7)
    if fold is True:
        gc_lats = np.fabs(gc_lats)
    gc_n, gc_bins = np.histogram(gc_lats, bins=bins, range=(min_lat,max_lat), normed=True, weights=1./np.cos(gc_lats))
    hist_ax.axhline(1.0, color='black', ls=':', lw=3)
    plt.draw()

    print("Getting lons and lats for mapped features")
    map_lons, map_lats, map_lens  = linlatstats(maplins)
    if fold is True:
        map_lats = np.fabs(map_lats)
    map_n, map_bins = np.histogram(map_lats, bins=bins, range=(min_lat, max_lat), normed=True, weights=map_lens/np.cos(map_lats))
    hist_ax.plot(plot_lats, map_n/gc_n.mean(), color='red', label='map', lw=2)
    plt.draw()

    print("Getting lons and lats for synthetic features")
    if synlins is None:
        synlins = [ lin for lin in load_lins('output/lins/synth_pool') if lin.length > 0 ]
    syn_lons, syn_lats, syn_lens = linlatstats(synlins)
    if fold is True:
        syn_lats = np.fabs(syn_lats)
    syn_n, syn_bins = np.histogram(syn_lats, bins=bins, range=(min_lat,max_lat), normed=True, weights=syn_lens/np.cos(syn_lats))
    hist_ax.plot(plot_lats, syn_n/gc_n.mean(), color='green', label='NSR', lw=2)
    plt.draw()

    hist_ax.grid(True)
    hist_ax.legend()
    hist_ax.set_xlabel("latitude")
    hist_ax.set_ylabel("lineament density (relative to uniform)")
    hist_ax.set_title("Longitudinally averaged lineament density")
    plt.draw()

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,'LinLatStatsCompare.'+save_fmt))

#}}}

def QvDbarMax(synthlins, gclins, crazylins, maplins, nq=26, use_stress=False, outfile='QvDbarMax'): #{{{
    """
    Create a plot exploring how valid our fitness function is.

    For each of the following maps, plot Q(map) as a function of Dbar_max:

    * Synthetic NSR
    * Synthetic NSR (lat > 30)

    * Great Circles
    * Great Circles (lat > 30)

    * Mapped features
    * Mapped features (lat > 30)

    * Jackstraws
    * Jackstraws (lat > 30)

    Hypothesis: The only reason the NSR features do better than the GCs is the
    lack of badly oriented features in the equatorial regions.

    Both synthetic NSR and Great Circle segments have essentially zero
    sinuosity.  The major difference between them is that GCs have some
    features in the equatorial regions that are oriented "wrong" to be due to
    NSR.  Both NSR and GCs should have similar Q() vs. Dbar_max curves when
    only high latitudes are considered, but GCs should be uniformly lower when
    the entire globe is included.

    Hypothesis: The mapped and the randomized (jackstraws) features are
    indistinguishable from each other, as far as fit quality is concerned,
    regardless of what value of Dbar_max is used.

    """

    dbar_maxes = np.linspace(0.0,0.25,nq)

    the_fig = plt.figure(figsize=(10,10))
    the_ax = the_fig.add_subplot(1,1,1)
    the_ax.grid()


    for lins,label,color in zip((synthlins, gclins, crazylins, maplins),\
                                ('Synthetic NSR', 'Great Circles', 'Randomized Map', 'Mapped'),\
                                ('red','green','blue','black')):
        hilat_lins = [ lin for lin in lins if np.degrees(np.fabs(lin.lats).min()) > 30 ]

        lins_Q = [ nsrfit_Q(lins, use_stress=use_stress, dbar_max=dbar_max) for dbar_max in dbar_maxes ]
        lins_Q_hilat = [ nsrfit_Q(hilat_lins, use_stress=use_stress, dbar_max=dbar_max) for dbar_max in dbar_maxes ]
        the_ax.plot(dbar_maxes, lins_Q, linewidth=3, label=label, color=color)
        the_ax.plot(dbar_maxes, lins_Q_hilat, linewidth=3, linestyle='dashed', label=label+r' (lat>30$^\circ$)', color=color)
        plt.draw()

    the_ax.set_xlabel(r'$\bar{D}_{max}$')
    the_ax.set_xticks(np.linspace(0,0.25,11))
    the_ax.set_ylabel('Q')
    the_ax.legend(loc='lower right')
    titlestr = r'Fit quality Q vs. $\bar{D}_{max}$ for several maps'
    if use_stress is True:
        titlestr += ' (stress weighted)'
        the_ax.set_yticks(np.linspace(0,1.2,13))
    else:
        titlestr += ' (unweighted)'
        the_ax.set_yticks(np.linspace(0,1.0,11))
    the_ax.set_title(titlestr)
    plt.draw()

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,outfile+'.'+save_fmt))

#}}}

###############################################################################
#    Plotting Helper Functions
###############################################################################

def show_extreme_doppels(lin, dbar_max=0.125): #{{{
    mp_lon, mp_lat = lin.bfgcseg_midpoint()
    seg_len = min(0.01, lin.length/10.0)
    best_fit = min(lin.nsrdbars)
    best_b = lin.bs[where(fabs(lin.nsrdbars - best_fit) < 1e-9)[0][0]]
    best_doppel = lineament.lingen_nsr(stresscalc=lin.stresscalc, init_lon=mp_lon+best_b, init_lat=mp_lat, prop_dir="both", max_length=lin.length, seg_len=seg_len)
    best_doppel.lons -= best_b

    worst_fit = max(lin.nsrdbars)
    worst_b = lin.bs[where(fabs(lin.nsrdbars - worst_fit) < 1e-9)[0][0]]
    worst_doppel = lineament.lingen_nsr(stresscalc=lin.stresscalc, init_lon=mp_lon+worst_b, init_lat=mp_lat, prop_dir="both", max_length=lin.length, seg_len=seg_len)
    worst_doppel.lons -= worst_b

    lines, map = lineament.plotlinmap([lin,], color='black', linewidth=2.0)
    map.scatter(degrees(lin.lons), degrees(lin.lats), color='black')

    lineament.plotlinmap([best_doppel,], color='green', map=map, linewidth=2)
    #map.scatter(degrees(best_doppel.lons), degrees(best_doppel.lats), color='green')
    #map.scatter(degrees(best_doppel.lons+2*pi), degrees(best_doppel.lats), color='green')
    #map.scatter(degrees(best_doppel.lons-2*pi), degrees(best_doppel.lats), color='green')
    gca().annotate('D=%.3g' % (best_fit,), xy=(degrees(best_doppel.lons[0]),      degrees(best_doppel.lats[0])) )
    gca().annotate('D=%.3g' % (best_fit,), xy=(degrees(best_doppel.lons[0]+2*pi), degrees(best_doppel.lats[0])) )
    gca().annotate('D=%.3g' % (best_fit,), xy=(degrees(best_doppel.lons[0]-2*pi), degrees(best_doppel.lats[0])) )

    lineament.plotlinmap([worst_doppel,], color='red', map=map, linewidth=2)
    #map.scatter(degrees(worst_doppel.lons), degrees(worst_doppel.lats), color='red')
    #map.scatter(degrees(worst_doppel.lons+2*pi), degrees(worst_doppel.lats), color='red')
    #map.scatter(degrees(worst_doppel.lons-2*pi), degrees(worst_doppel.lats), color='red')
    gca().annotate('D=%.3g' % (worst_fit,), xy=(degrees(worst_doppel.lons[0]),      degrees(worst_doppel.lats[0])) )
    gca().annotate('D=%.3g' % (worst_fit,), xy=(degrees(worst_doppel.lons[0]+2*pi), degrees(worst_doppel.lats[0])) )
    gca().annotate('D=%.3g' % (worst_fit,), xy=(degrees(worst_doppel.lons[0]-2*pi), degrees(worst_doppel.lats[0])) )

    return(best_doppel, worst_doppel)
#}}}

def fitcurve(lin, color='black', ax=None, dbar_max=0.125, show_dbar=True, show_w=True, show_fnsr=True, use_stress=True): #{{{
    """
    Plot delta and dbar as a function of b for lin, and show the value of the
    weighting function that results.

    """

    if ax is None:
        the_fig = figure(figsize=(9,6))
        dbar_ax = the_fig.add_subplot(1,1,1)
    else:
        dbar_ax = ax

    fit_ax  = twinx(ax=dbar_ax)

    # Plot the mapped lineament's fit curve:
    dbar_ax.set_xlabel('b')
    if show_dbar is True:
        dbar_ax.plot(degrees(lin.bs), lin.nsrdbars, ls='-', linewidth=2, color=color)
        dbar_ax.grid()
        dbar_ax.axis([-90,90,0,0.3])
        dbar_ax.set_yticks(linspace(0,0.3,5,endpoint=False))
        dbar_ax.set_ylabel('Dbar')

    if show_fnsr is True:
        fit_ax.fill_between(degrees(lin.bs), lin.nsrfits(dbar_max=dbar_max, use_stress=use_stress), 0, color=color, alpha=0.3)
        fit_ax_label = "H(b)"

    if show_w is True:
        fit_ax.plot(degrees(lin.bs), lin.nsrstresswts, ls='--', linewidth=2, color=color)

    if show_fnsr is True or show_w is True:
        fit_ax.set_yticks(linspace(0,1.5,5,endpoint=False))
        fit_ax.grid()
        fit_ax.axis([-90,90,0,1.5])
        fit_ax.set_xticks(linspace(-90,90,7))
        fit_ax.set_ylabel(fit_ax_label)

    dbar_ax.grid(True)
    fit_ax.grid(True)
#}}}

def activity_history(lins_list, the_fig=None, labels=[], colors=[], alphas=[], lin_cm=cm.jet, norm_by_all=False, dbar_max=0.125, verbose=True, titlestr="Apparent Lineament Activity History", outfile=None): #{{{
    """
    Plots apparent activity histories for one or several sets of lineaments,
    allowing visual comparison.

    lins_list is a list of lists of Lineament objects, assumed to have their
    fits already calculated.

    the_fig is the matplotlib figure in which the activity history should be
    drawn.  If None, a new plot is created.

    labels and colors are lists of strings and matplotlib color specifications
    that will be used to in legend creation and drawing (colors are
    automatically generated from the lin_cm colormap if colors is empty)

    if norm_by_sum is True, then the lineaments are treated as being subsets of
    one larger dataset, and their apparent contributions to the activity
    history are normalized by the cumulative length of the dataset, otherwise,
    they are assumed to be separate datasets, and each is normalized by its own
    cumulative length.

    """

    # Set hist_ax to the current axes, if none was supplied.
    if the_fig is None:
        the_fig = figure(figsize=(12,8))

    if len(the_fig.axes) == 0:
        hist_ax = the_fig.add_subplot(1,1,1)
    else:
        hist_ax = the_fig.axes[0]

    # Make sure that if we got labels and colors, there's the right number:
    if (len(labels) == 0):
        labels = [ 'Series '+str(N) for N in arange(len(lins_list))+1 ]
    try:
        assert(len(lins_list) == len(labels))
    except AssertionError:
        print "len(lins_list) = %g != len(labels) = %g" % (len(lins_list), len(labels))

    if (len(colors) == 0):
        colors = lin_cm(linspace(0,1,len(lins_list)))
    try:
        assert(len(lins_list) == len(colors))
    except AssertionError:
        print "len(lins_list) = %g != len(colors) = %g" % (len(lins_list), len(colors))

    if (len(alphas) == 0):
        alphas = np.ones(len(lins_list))
    try:
        assert(len(lins_list) == len(alphas))
    except AssertionError:
        print "len(lins_list) = %g != len(alphas) = %g" % (len(lins_list), len(alphas))

    # Depending on whether we're segmenting a single dataset (say, by lineament
    # length) or trying to compare different datasets, we may want to norm each
    # set of lineaments by its own length, or the total length of all the
    # lineaments.
    norm_lengths = [ sum([lin.length for lin in lins]) for lins in lins_list ]
    if norm_by_all is True:
        norm_lengths = [sum(norm_lengths),] * len(lins_list)

    for lins,label,color,alpha,norm_length in zip(lins_list,labels,colors,alphas,norm_lengths):
        if verbose is True:
            print("Calculating activity histories for %s" % (label,))
        acthist = calc_acthist(lins, dbar_max=dbar_max, norm_length=norm_length)
        x = np.concatenate([degrees(lins[0].bs)-180, degrees(lins[0].bs), degrees(lins[0].bs)+180])
        y = np.tile(acthist,3)
        hist_ax.plot(x, y, label=label, lw=3, c=color, alpha=alpha)

    hist_ax.set_ylabel("H(b)")
    hist_ax.set_xlabel("longitudinal translation b")
    hist_ax.grid(True)
    hist_ax.set_xticks(linspace(-90,90,19))
    degree_formatter = matplotlib.ticker.FormatStrFormatter(r'%.0f$^\circ$')
    hist_ax.xaxis.set_major_formatter(degree_formatter)
    hist_ax.set_xlim(-90,90)
    hist_ax.legend(loc="upper left")
    if norm_by_all is True:
        titlestr += ' (Normalized)'
    hist_ax.set_title(titlestr)

    # show the results, so we don't have some long awkward wait...
    draw()

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,outfile+'.'+save_fmt))
#}}}

###############################################################################
# Testing functions... #{{{
###############################################################################

def test_trimming(stresscalc, init_lon=0, init_lat=np.radians(60), max_length=0.5, prop_dir='east', seg_len=0.1, num_subsegs=10, color='red'): #{{{2
    """
    Testing lineament trimming (getting rid of overshoot on synthetic features)

    """
    perfect = lineament.lingen_nsr(stresscalc, init_lon=init_lon, init_lat=init_lat,\
                                   max_length=max_length, prop_dir=prop_dir,\
                                   seg_len=0.01, num_subsegs=100)
    testlin = lineament.lingen_nsr(stresscalc, init_lon=init_lon, init_lat=init_lat,\
                                   max_length=max_length, prop_dir=prop_dir,\
                                   seg_len=seg_len, num_subsegs=num_subsegs)
    the_fig = figure(figsize=(8,8))
    the_ax = the_fig.add_subplot(1,1,1)
    the_ax.plot(np.degrees(perfect.lons), np.degrees(perfect.lats), lw=2, color='black')
    the_ax.scatter(np.degrees(testlin.lons), np.degrees(testlin.lats), s=10, color=color)
    the_ax.grid(True)
    the_ax.set_title("target: %g, perfect: %g, rough: %g" % (max_length, perfect.length, testlin.length))
    draw()
#}}}2

def test_reckon_speed(N=1, M=100): #{{{2
    """
    A function for testing how long it takes to do spherical reckoning of
    multiple points.

    """

    lons     = 2*np.pi*np.random.random(size=N)
    lats     = np.pi*(np.random.random(size=N)-0.5)
    azimuths = 2*np.pi*np.random.random(size=N)
    dists    = np.pi*np.random.random(size=N)

    lons_out_list = []
    lats_out_list = []
    for lon,lat,az,d in zip(lons,lats,azimuths,dists):
        lons_out, lats_out = lineament.spherical_reckon(lon, lat, az, np.linspace(0,d,M))
        lons_out_list.append(lons_out)
        lats_out_list.append(lats_out)

    return(lons_out_list, lats_out_list)
#}}}2

def test_KDTree(linlib=None, libsize=30, N=10, d_max=1.0): #{{{2
    """
    A routine for testing how well this whole idea of finding an approximate
    doppelganger in the library works.

    """
    from scipy.spatial import KDTree

    radius = 1.0
    # generate a library if we didn't get one passed in:
    if linlib is None:
        linlib = lineament.lingen_nsr_library(nlats=libsize)

    # make a list of all the longitude and latitude points in the library
    lib_lons = concatenate([lin.lons for lin in linlib])
    lib_lats = concatenate([lin.lats for lin in linlib])
    lib_x, lib_y, lib_z = lineament.sphere2xyz(radius, pi/2-lib_lats, lib_lons)
    lib_kdt = KDTree(array([lib_x, lib_y, lib_z]).T)

    test_lons, test_lats = array(lineament.random_lonlatpoints(N))
    test_x, test_y, test_z = lineament.sphere2xyz(radius, pi/2.0-test_lats, test_lons)
    dists, near_idx = lib_kdt.query(array([test_x, test_y, test_z]).T, distance_upper_bound=d_max)

    near_idx = near_idx[where(dists<=d_max)]
    dists = dists[where(dists<=d_max)]
    return(dists, near_idx)

    #near_lons = mod(lib_lons[:,near_idx],2*pi)
    #near_lats = lib_lats[:,near_idx]
    #return(near_lons, near_lats)

    #the_fig=figure(figsize=(10,5))
    #map = the_fig.add_subplot(1,1,1)
    #map.scatter(degrees(mod(lib_lons,2*pi)),  degrees(lib_lats),  c='black', linewidth=0, s=3)
    #colors = cm.jet(linspace(0,1,len(test_lons)))
    #map.scatter(degrees(mod(test_lons,2*pi)), degrees(test_lats), c=colors,  linewidth=0)
    #map.scatter(degrees(mod(near_lons,2*pi)), degrees(near_lats), c=colors,  linewidth=0, alpha=0.7)
    #map.set_xlim([0,360])
    #map.set_ylim([-90,90])
    #return(linlib, near_lons, near_lats)

    #return(nearest_pts)
#}}}2

def test_fastfit(libsize=90, linlib=None, d_max=0.01): #{{{2

    print("Loading and updating mapped features")
    maplins = load_lins(os.path.join(lindir,'map_nsrfit'))
    maplins = lineament.update_lins(maplins)
    lz = maplins[0]

    if linlib is None:
        print("Generating lineament library with N=%d" % (libsize,) )
        linlib = lineament.lingen_nsr_library(nlats=libsize)

    print("Fitting to NSR directly")
    lz.calc_nsrfits()
    plot(degrees(lz.bs), lz.nsrdbars, linewidth=2, color='black')
    print("Fitting to NSR using linlib")
    lz.calc_nsrfits(doppel_library=linlib)
    plot(degrees(lz.bs), lz.nsrdbars, linewidth=2, color='red')

#}}}2

def test_fitres(lin, dbar_max=0.125, d_max=None, inter_res=20.0, nb=180): #{{{2
    """
    A function for testing how the resolution of synthetic features
    (doppelgangers) affects the accuracy of our fit metric.

    """
    europa = satstress.Satellite(open("input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",'r'))
    NSR = satstress.StressCalc([satstress.NSR(europa),])

    the_fig = figure(figsize=(10,10))
    map_ax = the_fig.add_subplot(2,1,1)
    ax1 = the_fig.add_subplot(2,1,2)

    map_ax.set_aspect('equal')
    map_ax.plot(np.degrees(lin.lons), np.degrees(lin.lats), lw=2, color='black')
    map_ax.scatter(np.degrees(lin.lons), np.degrees(lin.lats), color='black', s=5)
    init_lon, init_lat = lin.bfgcseg_midpoint()
    map_ax.grid()

    if d_max is None:
        d_max = np.median(lin.seg_lengths())

    print("Performing low-resolution fits (d_max=%g)" % (d_max,))
    lowres = lin.lonshift(radians(60))
    lowres.calc_nsrfits(stresscalc=NSR, d_max=d_max, nb=nb)
    ax1.plot(np.degrees(lowres.bs), lowres.nsrfits(use_stress=False), color='red', lw=2, label="$d_{max}=%g$" % (d_max,))
    ax1.fill_between(np.degrees(lowres.bs), lowres.nsrfits(use_stress=False), 0, color='red', alpha=0.3)
    ax1.set_ylim(0,1)
    ax1.set_yticks(linspace(0,1.0,11))
    ax1.grid()

    lowres_init_lon, lowres_init_lat = lowres.bfgcseg_midpoint()
    lowres_best_dbar, lowres_best_b = lowres.best_fit(use_stress=False)
    lowres_best_doppel = lineament.lingen_nsr(stresscalc=NSR, init_lon=lowres_init_lon+lowres_best_b, init_lat=lowres_init_lat, prop_dir="both", max_length=lowres.length, seg_len=d_max)
    lrbd = lowres_best_doppel.lonshift(-lowres_best_b-radians(60))
    map_ax.plot(np.degrees(lrbd.lons), np.degrees(lrbd.lats), lw=2, color='red')
    map_ax.scatter(np.degrees(lrbd.lons), np.degrees(lrbd.lats), color='red', s=5)

    the_fig.show()

    print("Performing high-resolution fits (d_max=%g)" % (d_max/inter_res,))
    highres = lin.lonshift(radians(120))
    highres.calc_nsrfits(stresscalc=NSR, d_max=d_max/inter_res, nb=nb)
    ax1.plot(np.degrees(highres.bs), highres.nsrfits(use_stress=False), color='green', lw=2, label="$d_{max}=%g$" % (d_max/inter_res,))
    ax1.fill_between(np.degrees(highres.bs), highres.nsrfits(use_stress=False), 0, color='green', alpha=0.3)
    ax1.set_ylim(0,1)
    ax2 = ax1.twinx()
    ax2.set_yticks(linspace(0,1.0,11))
    ax1.set_title("L=%g np=%d <l>=%g" % (lin.length, len(lin.lons), np.median(lin.seg_lengths()) ) )
    ax1.legend()

    highres_init_lon, highres_init_lat = highres.bfgcseg_midpoint()
    highres_best_dbar, highres_best_b = highres.best_fit(use_stress=False)
    highres_best_doppel = lineament.lingen_nsr(stresscalc=NSR, init_lon=highres_init_lon+highres_best_b, init_lat=highres_init_lat, prop_dir="both", max_length=highres.length, seg_len=d_max/inter_res)
    hrbd = highres_best_doppel.lonshift(-highres_best_b-radians(120))
    map_ax.plot(np.degrees(hrbd.lons), np.degrees(hrbd.lats), lw=2, color='green')
    map_ax.scatter(np.degrees(hrbd.lons), np.degrees(hrbd.lats), color='green', s=5)
    map_ax.scatter([np.mod(np.degrees(init_lon),360),], [np.degrees(init_lat),], color='black', marker='x', s=50)

    the_fig.show()
#}}}2

def test_initdist(lins): #{{{2
    """
    Look at the correlation between the distance between the fracture
    initiation point and how good a lineament's best fit is.

    """

    init_lons, init_lats = array([ lin.bfgcseg_midpoint() for lin in lins ]).transpose()

    mp_lons_list = []
    mp_lats_list = []
    for lin in lins:
        mp_lons, mp_lats = lin.seg_midpoints()
        mp_lons_list.append(mp_lons)
        mp_lats_list.append(mp_lats)

    lin_lengths = array([ lin.length for lin in lins ])
    d_mins = array([ lineament.spherical_distance(init_lon, init_lat, mp_lons, mp_lats).min() for init_lon,init_lat,mp_lons,mp_lats in zip(init_lons,init_lats,mp_lons_list,mp_lats_list) ])
    d_mins = d_mins/lin_lengths

    best_fits, best_bs = array([ lin.best_fit(use_stress=False) for lin in lins ]).transpose()

    the_fig = figure(figsize=(10,10))
    ax1 = the_fig.add_subplot(1,1,1)
    colors = cm.jet(linspace(0,1,11))
    ax1.scatter(d_mins,best_fits, s=[100*lin.length for lin in lins],alpha=0.5,lw=0,color='black')
    ax1.set_xlabel("Normalized dist. from crack init. pt. to prototype feature")
    ax1.set_ylabel("Quality of best fit ($0 \leq \max(f_{nsr}(b)) \leq 1$)")
    ax1.grid()
    ax1.set_title("symbol size proportional to prototype length ($R^2=%g$)" % (corrcoef(d_mins, best_fits)[0][1],))
    the_fig.show()
#}}}2

def test_mhd_by_lat(lin): #{{{2

    from scipy.optimize import brent

    ep1_lon, ep1_lat, ep2_lon, ep2_lat = lin.bfgcseg_endpoints()
    mp_lon, mp_lat = lineament.spherical_midpoint(ep1_lon, ep1_lat, ep2_lon, ep2_lat)
    max_length = lineament.spherical_distance(ep1_lon, ep1_lat, ep2_lon, ep2_lat)

    # We want to be able to compare the optimized results against the old
    # method, so we need to see what we got before:
    best_dbar, best_b = lin.best_dbar()

    brent_lats = []
    brent_dbars = []
    brent_iters = []
    brent_funcalls = []
    nsegs = linspace(2,100,99)
    for nseg in nsegs:
        brent_lat, brent_dbar, brent_iter, brent_funcall  = brent(lineament.mhd_by_lat, args=(mp_lon, lin.stresscalc, nseg, lin, max_length, best_b), full_output=True)
        brent_lats.append(brent_lat)
        brent_dbars.append(brent_dbar)
        brent_iters.append(brent_iter)
        brent_funcalls.append(brent_funcall)

    the_fig = plt.figure(figsize=(10,10))
    subplots_adjust(hspace=0.001)
    dbar_ax = the_fig.add_subplot(211)
    lat_ax = the_fig.add_subplot(212, sharex=dbar_ax)

    dbar_ax.plot(nsegs, brent_dbars)
    lat_ax.plot(nsegs, degrees(brent_lats))
    lat_ax.axhline(degrees(mp_lat))
    setp(dbar_ax.get_xticklabels(), visible=False)
    setp(lat_ax.get_yticklabels()[-1], visible=False)
    #dbar_ax.set_xlabel('nsegs')
    #ax.set_xticks(linspace(0,100,11))

    #ax.axvline(degrees(mp_lat), color='black')
    #ax.axhline(best_dbar, color='black')

    #ax.axvline(degrees(brent_lat), color='green')
    #ax.axhline(brent_dbar, color='green')

    #ax.set_ylim(0,0.25)
    dbar_ax.grid(True)
    lat_ax.grid(True)
 #}}}2 end test_mhd_by_lat

def test_brent_init(lin, nb=180, init_dop_res=[], dop_res=[]): #{{{2
    """
    Takes a lineament, and assumes it has its fits calculated using the old
    bfgc_midpoint initiation point.  Plots those fits, and then re-calculates
    them using the new optimized initiation point, for comparison.

    """
    from time import time

    the_fig = plt.figure()
    ax = the_fig.add_subplot(111)
    ax.set_ylim(0.0,1.001)
    ax.grid()
    draw()

    best_fits = []
    best_bs = []
    timing = []

    for idr,dr,n in zip(init_dop_res, dop_res, range(len(dop_res))):
        print("Fitting w/ IDR=%g, DR=%g" % (idr,dr))
        testlin = lineament.Lineament(lons=lin.lons, lats=lin.lats, stresscalc=lin.stresscalc)
        t0 = time()
        testlin.calc_nsrfits(nb=nb, init_doppel_res=idr, doppel_res=dr)
        t1 = time()
        dt = t1-t0
        timing.append(dt)
        best_fit, best_b = testlin.best_fit(use_stress=False)
        best_fits.append(best_fit)
        best_bs.append(best_b)
        print("    best fit: %g\n      best b: %g deg\n        time: %d s" % (best_fit, degrees(best_b), int(dt)))
        ax.plot(np.degrees(testlin.bs), testlin.nsrfits(use_stress=False), color=cm.gray(float(n)/len(dop_res)))
        draw()
#}}}2 end test_brent_init

#}}} end testing functions
###############################################################################

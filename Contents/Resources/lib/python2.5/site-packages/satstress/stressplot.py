#!/usr/bin/python
"""
Create map-projected plots of the tidal stresses on icy satellites.

This module allows the display of tidal stresses calculated using SatStress as
either scalar or vector fields, either as regularly gridded points on the
surface of the satellite, or arbitrary lists of points selected by the user.

It includes several test routines which also serve as examples of how to use
the module.

"""

import sys
import numpy as np
import satstress as ss
from matplotlib import pyplot as plt
from mpl_toolkits import basemap

def scalar_grid(stresscalc=None,\
                nlons=181, nlats=181,\
                min_lon=0.0, max_lon=np.pi,\
                min_lat=-np.pi/2, max_lat=np.pi/2,\
                time_t=0.0, field='tens', basemap_ax=None,\
                cmap=plt.cm.jet): #{{{
    """
    Display a rasterized scalar stress field defined by stresscalc, at the
    resolution specified by nlons, nlats, within the box bounded by min_lon,
    max_lon and min_lat, max_lat, at a time defined by time_t.  Which scalar
    field is displayed is controlled by field, which may be set to any of the
    following strings:

        'tens' -- magnitude of the more tensile principal component
        'comp' -- magnitude of the more compressive principal component
         'Ttt' -- meridional component of the stress tensor (north-south)
         'Tpt' -- off-diagonal component of the stress tensor
         'Tpp' -- east-west component of the stress tensor
    'w_stress' -- significance of the stresses, viz. magnitude, isotropy

    The time of the calculation, in seconds after periapse, is defined by
    time_t.

    The axes the stresses are plotted on is defined by basemap_ax.

    cmap defines the colormap used to color the scalar field.

    Returns the 2-dimensional array containing the plotted raster data, for use
    in making a scale bar, or other manipulation.

    """
    # TODO: need to add a scale bar to the RHS of the plot.

    calc_phis, calc_thetas = np.meshgrid(np.linspace(min_lon, max_lon, nlons), (np.pi/2.0)-np.linspace(min_lat, max_lat, nlats))
    calc_phis   = np.ravel(calc_phis)
    calc_thetas = np.ravel(calc_thetas)

    # some of the possible fields are easier to compute with the principal
    # components:
    if field=='tens' or field=='comp' or field=='w_stress' or field=='mean' or field=='diff':
        tens_mag, tens_az, comp_mag, comp_az = stresscalc.principal_components(calc_thetas, calc_phis, time_t)
        if nlons > 0 and nlats > 0:
            tens_mag = tens_mag.reshape(nlats, nlons)
            comp_mag = comp_mag.reshape(nlats, nlons)
            comp_az = comp_az.reshape(nlats, nlons)
            tens_az = tens_az.reshape(nlats, nlons)
        if field=='w_stress':
            w_stress = (tens_mag - comp_mag)/stresscalc.mean_global_stressdiff()

    # Or if people want to see the raw tensor components we can do that too:
    if field=='Ttt' or field=='Tpt' or field=='Tpp':
        Ttt, Tpt, Tpp = stresscalc.tensor(calc_thetas, calc_phis, time_t)
        if nlons > 0 and nlats > 0:
            Ttt = Ttt.reshape(nlats, nlons)
            Tpt = Tpt.reshape(nlats, nlons)
            Tpp = Tpp.reshape(nlats, nlons)
    
    # Now we need to display the results of our calculations, For a gridded
    # calculation, we can just show a raster with imshow()
    plot_field = None
    if field=='tens':
        plot_field=tens_mag
    if field=='comp':
        plot_field=comp_mag
    if field=='w_stress':
        plot_field=w_stress
    if field=='Ttt':
        plot_field=Ttt
    if field=='Tpt':
        plot_field=Tpt
    if field=='Tpp':
        plot_field=Tpp
    if field=='mean':
        plot_field=(tens_mag + comp_mag) / 2
    if field=='diff':
        plot_field=tens_mag - comp_mag
    basemap_ax.imshow(plot_field, cmap=cmap)
    return(plot_field)

#}}} end scalar

def scalar_points(stresscalc=None, lons=None, lats=None, time_t=0.0,\
                  field='tens', basemap_ax=None, cmap=plt.cm.jet, symb_sizes=50): #{{{
    """
    Display the magnitude of a scalar stress field defined by stresscalc at
    several discrete points defined by lons, lats, and a time defined by
    time_t.  Which scalar field is displayed is controlled by field, which may
    be set to any of the following strings:

        'tens' -- magnitude of the more tensile principal component
        'comp' -- magnitude of the more compressive principal component
         'Ttt' -- meridional component of the stress tensor (north-south)
         'Tpt' -- off-diagonal component of the stress tensor
         'Tpp' -- east-west component of the stress tensor
    'w_stress' -- significance of the stresses, viz. magnitude, isotropy

    The time of the calculation, in seconds after periapse, is defined by
    time_t.

    The axes the stresses are plotted on is defined by basemap_ax.

    cmap defines the colormap used to color the scalar field.

    """
    # TODO: add a scale bar to the RHS of the plot

    calc_phis   = lons
    calc_thetas = (np.pi/2.0) - lats

    # some of the possible fields are easier to compute with the principal
    # components:
    if field=='tens' or field=='comp' or field=='w_stress':
        tens_mag, tens_az, comp_mag, comp_az = stresscalc.principal_components(calc_thetas, calc_phis, time_t)
        if field=='w_stress':
            w_stress = (tens_mag - comp_mag)/stresscalc.mean_global_stressdiff()

    # Or if people want to see the raw tensor components we can do that too:
    if field=='Ttt' or field=='Tpt' or field=='Tpp':
        Ttt, Tpt, Tpp = stresscalc.tensor(calc_thetas, calc_phis, time_t)

    if field=='tens':
        field_plot = tens_mag
    if field=='comp':
        field_plot = comp_mag
    if field=='w_stress':
        field_plot = w_stress
    if field=='Ttt':
        field_plot = Ttt
    if field=='Tpt':
        field_plot = Tpt
    if field=='Tpp':
        field_plot = Tpp

    field_norm = plt.Normalize(vmin=np.min(field_plot),vmax=np.max(field_plot))

    basemap_ax.scatter(np.degrees(lons), np.degrees(lats), lw=0, s=symb_sizes, c=field_plot, cmap=cmap, norm=field_norm)

#}}}

def vector_points(stresscalc=None, lons=None, lats=None, time_t=0.0,\
           plot_tens=True, plot_comp=True, plot_greater=True, plot_lesser=True,\
           basemap_ax=None, lonshift=0, w_stress=False,\
           scale=1e8, scale_arr=None, arrow_width=0.0015): #{{{
    """
    Display the principal components of the tidal stresses defined by the input
    stresscalc object at the points defined by lons and lats, which are one
    dimensional arrays of equal length, and a time defined by time_t, in
    seconds after periapse.

    The stress vectors are plotted on the map axes defined by basemap_ax.

    By default all the principal components are plotted, but if you wish to see
    only the more or less tensile (less or more compressive) or only those
    principal components which are absolutely compressive or tensile, you may
    exclude some subset of the vectors using the following flags:

       plot_tens  --  if True, plot all tensile stresses.
       plot_comp  --  if True, plot all compressive stresses.
    plot_greater  --  if True, plot the greater (more tensile) principal component
     plot_lesser  --  if True, plot the lesser (more compressive) principal component

    lonshift is a longitudinal displacement added to lons when the stresses are
    calculated, useful in creating plots of lineaments at their current
    location, compared to stresses that they would have experienced at their
    apparent location of formation (i.e. those stresses which best match the
    feature).  For instance, if you wished to show only those stresses which are 
    the more tensile, and which are actually tensile, you would need to set
    the flags: plot_comp=False, plot_lesser=False.

    If w_stress is true, the lengths of the arrows which are used to represent
    the stresses are scaled according to how significant their location within
    the stress field is, i.e. large stresses and anisotropic stresses will be
    more prominent than small stresses and isotropic stresses.

    scale determines the overall size of the arrows representing the stresses.
    A smaller scale means bigger arrows.

    scale_arr is an array of the same length as lons and lats, which is used to
    scale the lengths of the vectors.  Useful in showing the relative
    importance of different segments of a feature having non-uniform lengths.

    arrow_width is passed in to np.quiver(), and is the width of the arrow
    shaft, as a proportion of the width of the overall plot.

    """

    calc_phis   = lons
    calc_thetas = (np.pi/2.0)-lats

    Ttt, Tpt, Tpp = stresscalc.tensor(calc_thetas, calc_phis+lonshift, time_t)

    Tau = np.array([[Ttt,Tpt],[Tpt,Tpp]])
    eigensystems = [ np.linalg.eig(Tau[:,:,N]) for N in range(len(Tau[0,0,:])) ]
    evals = np.array([ e[0] for e in eigensystems ])
    evecs = np.array([ e[1] for e in eigensystems ])

    eigval_A = evals[:,0]
    ex_A     = evecs[:,0,1]
    ey_A     = evecs[:,0,0]

    eigval_B = evals[:,1]
    ex_B     = evecs[:,1,1]
    ey_B     = evecs[:,1,0]

    mag1 = np.where(eigval_A >  eigval_B, eigval_A, eigval_B)
    ex1  = np.where(eigval_A >  eigval_B, ex_A, ex_B)
    ey1  = np.where(eigval_A >  eigval_B, ey_A, ey_B)

    mag2 = np.where(eigval_A <= eigval_B, eigval_A, eigval_B)
    ex2  = np.where(eigval_A <= eigval_B, ex_A, ex_B)
    ey2  = np.where(eigval_A <= eigval_B, ey_A, ey_B)

    if np.shape(scale_arr) != np.shape(mag1):
        scale_arr = np.ones(np.shape(mag1))
    if np.shape(w_stress) == np.shape(mag1):
        scale_arr = scale_arr*(mag1 - mag2)/stresscalc.mean_global_stressdiff()

    mag1_comp = np.ma.masked_where(mag1 > 0, mag1)
    mag1_tens = np.ma.masked_where(mag1 < 0, mag1)

    mag2_comp = np.ma.masked_where(mag2 > 0, mag2)
    mag2_tens = np.ma.masked_where(mag2 < 0, mag2)

    # First plot all the compressional stresses
    if plot_comp is True:
        if plot_greater is True:
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas),  mag1_comp*ex1*scale_arr,  mag1_comp*ey1*scale_arr,\
                              lw=0., width=arrow_width, scale=scale, color='blue', pivot='tip')
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas), -mag1_comp*ex1*scale_arr, -mag1_comp*ey1*scale_arr,\
                              lw=0., width=arrow_width, scale=scale, color='blue', pivot='tip')

        if plot_lesser is True:
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas),  mag2_comp*ex2*scale_arr,  mag2_comp*ey2*scale_arr,\
                              lw=0., width=arrow_width, scale=scale, color='blue', pivot='tip')
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas), -mag2_comp*ex2*scale_arr, -mag2_comp*ey2*scale_arr,\
                              lw=0., width=arrow_width, scale=scale, color='blue', pivot='tip')

    # Now all the tensional stresses:
    if plot_tens is True:
        if plot_greater is True:
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas),  mag1_tens*ex1*scale_arr,  mag1_tens*ey1*scale_arr,\
                              lw=0., width=arrow_width, scale=scale, color='red', pivot='tail')
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas), -mag1_tens*ex1*scale_arr, -mag1_tens*ey1*scale_arr,\
                              lw=0., width=arrow_width, scale=scale, color='red', pivot='tail')

        if plot_lesser is True:
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas),  mag2_tens*ex2*scale_arr,  mag2_tens*ey2*scale_arr,\
                              lw=0., width=arrow_width, scale=scale, color='red', pivot='tail')
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas), -mag2_tens*ex2*scale_arr, -mag2_tens*ey2*scale_arr,\
                              lw=0., width=arrow_width, scale=scale, color='red', pivot='tail')

#}}} end vector_points


################################################################################
# Test routines and examples of use:
################################################################################

def test_vector_grid(plot_tens=True, plot_comp=True, plot_greater=True, plot_lesser=True, min_lon=0.0, max_lon=np.pi, nlons=13, min_lat=-np.pi/2.0, max_lat=np.pi/2.0, nlats=13, t=0.0): #{{{
    """
    An attempt to exercise the above routines...

    """

    satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = ss.Satellite(open(satfile,'r'))
    NSR = ss.StressCalc([ss.NSR(europa),])

    grid_min_lon = min_lon
    grid_max_lon = max_lon
    grid_min_lat = min_lat
    grid_max_lat = max_lat
    vector_grid_nlons = nlons
    vector_grid_nlats = nlats

    # Lon/Lat locations to plot the principal components:
    vector_grid_lons  = np.linspace(grid_min_lon, grid_max_lon, vector_grid_nlons)
    vector_grid_lats  = np.linspace(grid_min_lat, grid_max_lat, vector_grid_nlats)
    vector_mesh_lons, vector_mesh_lats = np.meshgrid(vector_grid_lons, vector_grid_lats)
                                                     
    vector_mesh_lons = np.ravel(vector_mesh_lons)
    vector_mesh_lats = np.ravel(vector_mesh_lats)

    vector_grid_fig = plt.figure(figsize=(10,10))
    vector_grid_ax  = vector_grid_fig.add_subplot(1,1,1)
    vector_grid_ax.set_title("Regularly gridded vectors")
    vector_grid_basemap = basemap.Basemap(llcrnrlon = np.degrees(grid_min_lon),\
                                          llcrnrlat = np.degrees(grid_min_lat),\
                                          urcrnrlon = np.degrees(grid_max_lon),\
                                          urcrnrlat = np.degrees(grid_max_lat),\
                                          ax = vector_grid_ax)

    vector_points(stresscalc=NSR, lons=vector_mesh_lons, lats=vector_mesh_lats, time_t=t,\
                  plot_greater=plot_greater, plot_lesser=plot_lesser,\
                  plot_comp=plot_comp, plot_tens=plot_tens, basemap_ax=vector_grid_basemap)

    vector_grid_basemap.drawmeridians(np.degrees(vector_grid_lons), labels=[1,0,0,1], linewidth=0.5, color='gray')
    vector_grid_basemap.drawparallels(np.degrees(vector_grid_lats), labels=[1,0,0,1], linewidth=0.5, color='gray')
    vector_grid_basemap.drawmapboundary()
#}}} end test_vector_grid

def test_vector_lin(lin, t=0.0, w_stress=False, w_length=False, scale=5e7, plot_tens=True, plot_comp=False, plot_greater=True, plot_lesser=False): #{{{
    
    import lineament

    satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = ss.Satellite(open(satfile,'r'))
    NSR = ss.StressCalc([ss.NSR(europa),])

    best_b = lin.best_fit()[1]
    mp_lons, mp_lats = lin.seg_midpoints()

    if w_length is True:
        scale_arr = len(mp_lons)*lin.seg_lengths()/lin.length
    else:
        scale_arr = np.ones(np.shape(mp_lons))

    min_lat = 0
    max_lat = np.pi/2.0
    min_lon = 0.0
    max_lon = np.pi

    vector_lin_fig = plt.figure()
    vector_lin_ax  = vector_lin_fig.add_subplot(1,1,1)
    vector_lin_ax.set_title("Vectors at arbitrary locations (e.g. along a feature)")

    vector_lin_basemap = basemap.Basemap(llcrnrlon = np.degrees(min_lon),\
                                         llcrnrlat = np.degrees(min_lat),\
                                         urcrnrlon = np.degrees(max_lon),\
                                         urcrnrlat = np.degrees(max_lat),\
                                         ax = vector_lin_ax)

    junk = lineament.plotlinmap([lin,], map=vector_lin_basemap, lat_mirror=True, lon_cyc=np.pi)
    vector_lin_basemap.drawmeridians(np.degrees(np.linspace(min_lon, max_lon,13)), labels=[1,0,0,1], linewidth=0.5, color='gray')
    vector_lin_basemap.drawparallels(np.degrees(np.linspace(min_lat, max_lat, 7)), labels=[1,0,0,1], linewidth=0.5, color='gray')
    vector_lin_basemap.drawmapboundary()

    vector_points(stresscalc=NSR, lons=np.mod(mp_lons,np.pi), lats=mp_lats, time_t=0.0, lonshift=best_b,\
                  plot_tens=plot_tens, plot_greater=plot_greater, plot_comp=plot_comp, plot_lesser=plot_lesser, basemap_ax=vector_lin_basemap, scale=scale, scale_arr=scale_arr, w_stress=w_stress)

#}}} end test_vector_lin

def test_scalar_grid(min_lon=0.0, max_lon=np.pi, min_lat=-np.pi/2.0, max_lat=np.pi/2.0,\
                     nlats=181, nlons=181, field='tens', cmap=plt.cm.gray_r): #{{{

    satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = ss.Satellite(open(satfile,'r'))
    NSR = ss.StressCalc([ss.NSR(europa),])

    the_fig = plt.figure(figsize=(10,10))
    map_ax  = the_fig.add_subplot(1,1,1)
    map_ax.set_title("Rasterized scalar stress field (%s)" % (field,) )
    scalar_basemap = basemap.Basemap(llcrnrlon = np.degrees(min_lon),\
                                     llcrnrlat = np.degrees(min_lat),\
                                     urcrnrlon = np.degrees(max_lon),\
                                     urcrnrlat = np.degrees(max_lat),\
                                     ax = map_ax)

    scalar_grid(stresscalc=NSR, nlons=nlons, nlats=nlats, min_lon=min_lon, max_lon=max_lon, min_lat=min_lat, max_lat=max_lat,\
                field=field, cmap=cmap, basemap_ax=scalar_basemap)

    scalar_basemap.drawmeridians(np.degrees(np.linspace(min_lon, max_lon, 7)), labels=[1,0,0,1], linewidth=0.5, color='gray')
    scalar_basemap.drawparallels(np.degrees(np.linspace(min_lat, max_lat, 7)), labels=[1,0,0,1], linewidth=0.5, color='gray')
    scalar_basemap.drawmapboundary()
#}}} end test_scalar_grid

def test_scalar_lin(lin, field='tens', cmap=plt.cm.jet): #{{{

    import lineament

    satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = ss.Satellite(open(satfile,'r'))
    NSR = ss.StressCalc([ss.NSR(europa),])

    mp_lons, mp_lats = lin.seg_midpoints()
    seg_lengths = lin.seg_lengths()

    min_lat = 0
    max_lat = np.pi/2.0
    min_lon = 0.0
    max_lon = np.pi

    the_fig = plt.figure(figsize=(10,10))
    map_ax  = the_fig.add_subplot(1,1,1)
    map_ax.set_title("Point by Point scalar stresses field (%s)" % (field,) )
    scalar_basemap = basemap.Basemap(llcrnrlon = np.degrees(min_lon),\
                                     llcrnrlat = np.degrees(min_lat),\
                                     urcrnrlon = np.degrees(max_lon),\
                                     urcrnrlat = np.degrees(max_lat),\
                                     ax = map_ax)

    scalar_points(stresscalc=NSR, lons=np.mod(mp_lons,np.pi), lats=np.mod(mp_lats,np.pi/2.0), symb_sizes=5000*seg_lengths, field=field, cmap=cmap, basemap_ax=scalar_basemap)
    junk = lineament.plotlinmap([lin,], map=scalar_basemap, lat_mirror=True, lon_cyc=np.pi)

    scalar_basemap.drawmeridians(np.degrees(np.linspace(min_lon, max_lon, 7)), labels=[1,0,0,1], linewidth=0.5, color='gray')
    scalar_basemap.drawparallels(np.degrees(np.linspace(min_lat, max_lat, 7)), labels=[1,0,0,1], linewidth=0.5, color='gray')
    scalar_basemap.drawmapboundary()
#}}} end test_scalar_points


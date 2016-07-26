from numpy import *
from pylab import *
from mpl_toolkits.basemap import Basemap
import pickle

# Open Source Geospatial libraries:
from osgeo import ogr

class Lineament(object): #{{{1
    """
    A one dimensional feature on the surface of a spherical satellite.

    """

    def __init__(self, lons=None, lats=None, stresscalc=None, bs=None, nsrdbars=None, nsrstresswts=None, gid=None): #{{{2
        """
        Create a lineament from a given list of (lon,lat) points.

        lons and lats are arrays of longitudes and latitudes, having equal
        length, defining the vertices making up a polyline.  East and North are
        taken to be positive.  Units are radians
        
        stresscalc is a L{satstress.StressCalc} object defining the field we
        are going to be comparing the lineament to.
        
        If bs, nsrdbars, and nsrstresswts are set and have the same length, then we have
        a complete NSR fit and we can go ahead and set those values for the
        feature.  Otherwise, set these attributes to None.

        """

        # make sure we've got good input points...
        assert (len(lons) == len(lats))
        assert (len(lons) > 0)

        self.lons = array(lons)
        self.lats = array(lats)
        self.length = self.calc_length()
        self.stresscalc = stresscalc
        self.gid = gid
        self.hashval = self.calc_hash()

        if bs is None or nsrdbars is None or nsrstresswts is None:
            self.bs = None
            self.nsrdbars = None
            self.nsrstresswts = None
        else:
            assert(len(bs) == len(nsrdbars) == len(nsrstresswts))
            self.bs = array(bs)
            self.nsrdbars = array(nsrdbars)
            self.nsrstresswts = array(nsrstresswts)
    #}}}2

    def calc_hash(self): #{{{2
        """
        In order to be able to use a Lineament as a node in a NetworkX graph,
        it needs to be hashable.  The hash value depends on the all the
        longitude and latitude values of the lineament, and on their adjacency,
        but not on which order they appear in.  For example, if the lineament
        is composed of 5 vertices: A,B,C,D,E, and it's hash is X, then
        reversing the order of the verticies: E,D,C,B,A will result in the same
        hash, X, but re-arranging their ordering: C,D,B,A,E would result in
        some other hash value.
        """

        # Concatenate lons and their reverse, to get same hash either direction
        hash_lons = np.append(self.lons, self.lons[::-1])
        # Periodicity of longitude should not affect hash value:
        hash_lons = np.mod(hash_lons,2.0*np.pi)
        # Needs to be an integer output ultimately, and we don't want any floating point
        # inaccuracies screwing up the equality values:
        hash_lons = (1e9*hash_lons).astype('int')

        hash_lats = np.append(self.lats, self.lats[::-1])
        hash_lats = (1e9*hash_lats).astype('int')

        hash_ns = (1e9*np.tile(arange(len(hash_lons)),2)).astype('int')

        hash=0
        for lon, lat, n in zip(hash_lons, hash_lats, hash_ns):
            hash = hash^(lon-lat+n)

        hash = hash^(np.int(1e9*self.length))

        return(hash)
    #}}}2

    def __hash__(self): #{{{2
        """
        Return the lineament's cached hash value (see calc_hash)

        """

        try:
            hashval = self.hashval
        except(AttributeError):
            self.hashval = self.calc_hash()

        return(self.hashval)

    #}}}2

    def __cmp__(self, other): #{{{2
        """
        Basic comparison operator for Lineaments, depends on their __hash__()
        values.  A lineament is defined by a series of lat/lon points, and will
        compare equal if the same series of lat/lon points are used for both
        features (to within one part in 10^9).  The order of the points may be
        reversed, and they objects will still compare equal.

        """

        return(self.__hash__()-other.__hash__())
    #}}}2

    def __str__(self): #{{{2
        """
        This string method only gives an ID number for the feature (it's hash
        value).  If you want a real string representation of the feature's
        geographic information, use Lineament.wkt() to get the "well known
        text" format.

        """

        return(str(hash(self)))
    #}}}2

    def wkt(self): #{{{2
        """
        Outputs the geographic data associated with the lineament using the
        so-called "well known text" (wkt) representation, in order to allow
        easy output and use with the GDAL/OGR open source geospatial libraries.


        """

        linestr_geom_wkt = 'LINESTRING('
        for lon, lat in zip(self.lons, self.lats):
            linestr_geom_wkt = linestr_geom_wkt + "%f %f, " % (lon, lat)
        # get rid of the trailing comma and space...
        linestr_geom_wkt = linestr_geom_wkt[:-2]
        # close the parentheses
        linestr_geom_wkt = linestr_geom_wkt +')'

        return(linestr_geom_wkt)
#}}}2

    def calc_length(self): #{{{2
        """
        Return the total length of the L{Lineament} in radians of arc, along
        the surface of the sphere.

        """
        if len(self.lons) > 1:
            return(sum(self.seg_lengths()))
        else:
            return(0.0)

    #}}}2 end calc_length
    
    def sinuosity(self): #{{{2
        """
        Return the sinuosity of the lineament, using spherical distances.

        Sinuosity is the ratio of the sum of the lengths of all the line
        segments making up the lineament, to the distance separating the
        endpoints.

        """
        # If the two endpoints are the same, we'll get a divide by zero error...
        assert( (self.lons[0],self.lats[0]) != (self.lons[-1], self.lats[-1]))
        return (self.length / spherical_distance(self.lons[0], self.lats[0], self.lons[-1], self.lats[-1]))

    #}}}2

    def midpoint(self): #{{{2
        """
        Returns the (lon,lat) point halfway along the lineament.

        """

        seg_lengths = self.seg_lengths()
        half_length = self.length/2.0

        cumulative_lengths = concatenate([ [0.0,], array([ seg_lengths[:N+1].sum() for N in range(len(seg_lengths)) ])])
        just_past_half = where(cumulative_lengths > half_length)[0][0]
        passed_by = cumulative_lengths[just_past_half] - half_length
        back_az = spherical_azimuth(self.lons[just_past_half], self.lats[just_past_half], self.lons[just_past_half-1], self.lats[just_past_half-1])
        mp_lon, mp_lat = spherical_reckon(self.lons[just_past_half], self.lats[just_past_half], back_az, passed_by)

        return(mp_lon, mp_lat)

    #}}}2

    def seg_midpoints(self): #{{{2
        """
        Return a list of the midpoints of each line segment in the Lineament.

        """
        if len(self.lons) == 1:
            mp_lons, mp_lats = self.lons, self.lats
        else:
            mp_lons,mp_lats = spherical_midpoint(self.lons[:-1], self.lats[:-1], self.lons[1:], self.lats[1:])

        return(mp_lons, mp_lats)
    #}}}2

    def seg_azimuths(self): #{{{2
        """
        Calculate the lineament azimuths (orientations) at the midpoints of the
        segments making up the lineament.  Because these segments are not
        directional, this will be an angle between 0 and pi.

        """

        mplons, mplats = self.seg_midpoints()
        return(mod(spherical_azimuth(mplons, mplats, self.lons[:-1], self.lats[:-1]), pi))
    # }}}2 end midpoint_azimuths

    def seg_lengths(self): #{{{2
        """
        Calculate the lengths (in radians of arc) of the line segments making
        up the lineament.

        """

        return(spherical_distance(self.lons[:-1], self.lats[:-1], self.lons[1:], self.lats[1:]))

    #}}}2 end seg_lengths

    def lonshift(self, b): #{{{2
        """
        Return the lineament shifted in longitude by 'b' radians on the
        surface, with its fits as they would have been if they'd been
        calculated from that location.

        """
        # Add b to each of self.lons
        shift_lons = self.lons+b

        if self.bs is not None:
            shift_bs = mod(self.bs-b,pi)
            nb = len(self.bs)
            dtype = [('bs',float),('nsrdbars',float),('nsrstresswts',float)]
            shift_fits = zeros(nb, dtype=dtype)
            shift_fits['bs'] = shift_bs
            shift_fits['nsrdbars'] = self.nsrdbars
            shift_fits['nsrstresswts'] = self.nsrstresswts
            shift_fits.sort(order='bs')
            new_bs = shift_fits['bs']
            new_nsrdbars = shift_fits['nsrdbars']
            new_nsrstresswts = shift_fits['nsrstresswts']
        else:
            new_bs = None
            new_nsrdbars = None
            new_nsrstresswts = None

        return(Lineament(lons=shift_lons, lats=self.lats, stresscalc=self.stresscalc,\
                         bs=new_bs, nsrdbars=new_nsrdbars, nsrstresswts=new_nsrstresswts))
    #}}}2

    def poleshift(self, pnp_lon=0.0, pnp_lat=pi/2): #{{{2
        """
        Return a Lineament object representing the location and orientation of
        self, when the point defined in modern lon,lat terms by pnp_lon pnp_lat
        was the north pole.

        """
        tpw_lons, tpw_lats = paleopole_transform(lon_in=self.lons, lat_in=self.lats, pnp_lon=pnp_lon, pnp_lat=pnp_lat)
        return(Lineament(lons=tpw_lons, lats=tpw_lats, stresscalc=self.stresscalc))
    #}}}2

    def d_min(self, linB): #{{{2
        """
        Return an array containing the distances from the midpoints of each
        segment in self to the nearest midpoint of a segment in linB.

        """
        
        return(d_min(self, linB))
    #}}}2

    def mhd(self, linB): #{{{2
        """
        MHD form A to B.

        """
        
        return(mhd(self, linB))
    #}}}2

    def doppelgen_midpoint_nsr(self, stresscalc=None, init_doppel_res=0.0, doppel_res=0.1, num_subsegs=10): #{{{2
        """
        Find a midpoint, representative of the lineament, and generate a synthetic
        NSR feature there which approximates the feature, for each value of
        backrotation (b, longitudinal shift) stored in self.bs.  Return those
        doppelgangers as a list.

        doppel_res determines the resolution of the doppelgangers.  A smaller
        number means lower resolution.  doppel_res=1.0 means the mean segment
        length is the same as the prototype lineament.  doppel_res=2.0 would
        mean the mean segment length of the prototype lineament is twice that
        of the doppelganger, etc.

        If init_doppel_res is <= 0.0 then initial doppelgangers are not
        generated, and the midpoint of the great circle segment is used as the
        initiation point instead, saving considerable computation.

        """

        mean_seg_len = (self.length/(len(self.lons)-1))
        doppel_seg_len = mean_seg_len/doppel_res
        if init_doppel_res > 0.0:
            init_doppel_seg_len = mean_seg_len/init_doppel_res
            init_lons, init_lats, bfgc_len = best_nsr_init_points(self, stresscalc, seg_len=init_doppel_seg_len, num_subsegs=num_subsegs)
        else:
            bfgc_ep1_lon, bfgc_ep1_lat, bfgc_ep2_lon, bfgc_ep2_lat, bfgc_mp_lon, bfgc_mp_lat, bfgc_len = self.bfgcseg()
            init_lons = self.bs + bfgc_mp_lon
            init_lats = array([bfgc_mp_lat,]).repeat(len(self.bs))

        # Now we need to generate doppelgangers, which are perfect synthetic
        # features that have resulted from the NSR stress field.
        doppels = [ lingen_nsr(stresscalc=stresscalc, init_lon=init_lon, init_lat=init_lat,\
                        max_length=bfgc_len, prop_dir="both", seg_len=doppel_seg_len, num_subsegs=num_subsegs) for init_lon,init_lat in zip(init_lons,init_lats) ]

        return(doppels)

    #}}}2 end doppelgen_midpoint_nsr

    def doppelgen_gcseg(self): #{{{2
        """
        Generate a great circle segment approximating the lineament.

        """
        init_lon, init_lat, fin_lon, fin_lat = self.bfgcseg_endpoints()
        return(lingen_greatcircle(init_lon, init_lat, fin_lon, fin_lat))
    #}}}2

    def bfgcseg_endpoints(self): #{{{2
        """
        Find the points on the lineament's best fit great circle which are
        closest to the endpoints of the lineament.

        """

        mplons, mplats = self.seg_midpoints()
        bfgcpole_lon, bfgcpole_lat = self.bfgc_pole()

        # find points on the great circle nearest to the lineament endpoints.
        # Go mod((pi/2)-(distance to gcpole),pi/2) radians away from the pole if closer than pi/2
        # and toward the pole if further away than pi/2.
        dist_to_bfgcpole1 = spherical_distance(self.lons[0], self.lats[0], bfgcpole_lon, bfgcpole_lat)
        az_to_bfgcpole1   = spherical_azimuth(self.lons[0], self.lats[0], bfgcpole_lon, bfgcpole_lat)

        # if further than pi/2 from the gc_pole, go mod(dist,pi/2) toward it.
        if(dist_to_bfgcpole1 > pi/2):
            bfgcseg_lon1, bfgcseg_lat1 = spherical_reckon(self.lons[0], self.lats[0], az_to_bfgcpole1, mod(dist_to_bfgcpole1, pi/2))
        # otherwise, go pi/2-dist away from it.
        else:
            bfgcseg_lon1, bfgcseg_lat1 = spherical_reckon(self.lons[0], self.lats[0], az_to_bfgcpole1+pi, pi/2-dist_to_bfgcpole1)

        dist_to_bfgcpole2 = spherical_distance(self.lons[-1], self.lats[-1], bfgcpole_lon, bfgcpole_lat)
        az_to_bfgcpole2   = spherical_azimuth(self.lons[-1], self.lats[-1], bfgcpole_lon, bfgcpole_lat)

        # if further than pi/2 from the gc_pole, go mod(dist,pi/2) toward it.
        if(dist_to_bfgcpole2 > pi/2):
            bfgcseg_lon2, bfgcseg_lat2 = spherical_reckon(self.lons[-1], self.lats[-1], az_to_bfgcpole2, mod(dist_to_bfgcpole2, pi/2))
        # otherwise, go pi/2-dist away from it.
        else:
            bfgcseg_lon2, bfgcseg_lat2 = spherical_reckon(self.lons[-1], self.lats[-1], az_to_bfgcpole2+pi, pi/2-dist_to_bfgcpole2)

        return(bfgcseg_lon1, bfgcseg_lat1, bfgcseg_lon2, bfgcseg_lat2)
    #}}}2

    def bfgcseg_midpoint(self): #{{{2
        """
        Find the lon,lat point which lies on the lineament's best fit great
        circle, and is halfway between points on the great circle which are
        closest to the endpoints of the lineament.  This point represents a
        kind of geometric center of the feature, which can be used to
        approximate its location for a variety of purposes, including
        generating doppelgangers.

        """

        bfgcseg_lon1, bfgcseg_lat1, bfgcseg_lon2, bfgcseg_lat2 = self.bfgcseg_endpoints()
        bfgc_mp_lon, bfgc_mp_lat = spherical_midpoint(bfgcseg_lon1, bfgcseg_lat1, bfgcseg_lon2, bfgcseg_lat2)
        return(bfgc_mp_lon, bfgc_mp_lat)

    #}}}2

    def bfgcseg(self): #{{{2
        """
        Find the best fit great circle segment representing the lineament and
        return some information about it as a tuple in the following order:

        (ep1_lon, ep1_lat, ep2_lon, ep2_lat, mp_lon, mp_lat, bfgcseg_length)

        where ep1_lon, ep1_lat, ep2_lon, ep2_lat, are the lons and lats of the
        two endpoints, mp_lon and mp_lat are the lon and lat of the midopint,
        and seg_length is the geodesic distance from one end to the other.

        """
        ep1_lon, ep1_lat, ep2_lon, ep2_lat = self.bfgcseg_endpoints()
        mp_lon, mp_lat = spherical_midpoint(ep1_lon, ep1_lat, ep2_lon, ep2_lat)
        bfgcseg_length = spherical_distance(ep1_lon, ep1_lat, ep2_lon, ep2_lat)

        return(ep1_lon, ep1_lat, ep2_lon, ep2_lat, mp_lon, mp_lat, bfgcseg_length)
    #}}}2

    def calc_nsrfits(self, nb=180, stresscalc=None, init_doppel_res=0.0, doppel_res=0.1, num_subsegs=10): #{{{2
        """
        For nb evenly spaced values of longitudinal translation, b, ranging
        from 0 to pi, calculate the fit metric (dbar) for the lineament,
        using the stresscalc object associated with the feature to generate
        the hypothesized synthetic features.  Use the stresscalc that was
        passed in if there is none associated with the feature, or 

        Store the resulting values of b and dbar in the lineament's
        attributes bs, and nsrdbars for later reference.

        This only works for a purely NSR stress field.

        Calculate a metric of the significance of the feature's location within
        the stress field, based on its magnitude and anisotropy relative to the
        global average, and store it for later use in weighting the feature's
        contribution to the overall results.  This is self.nsrstresswts

        Generate a synthetic lineament meant to replicate the feature (a
        doppelganger) by using the NSR stress field, and a starting point near
        the center of the feature.  Calculate the mean Hausdorff distance (MHD)
        between the prototype (self) and this doppelganger for every value of
        b, normalized by the feature's length, and store the results in
        self.nsrdbars.

        """

        # we have to have at least one stresscalc or this is pointless:
        assert(self.stresscalc is not None or stresscalc is not None)

        if stresscalc is None and self.stresscalc is not None:
            stresscalc = self.stresscalc
        elif self.stresscalc is None and stresscalc is not None:
            self.stresscalc = stresscalc

        # set the b values first, so that the fit metrics can refer to them.
        self.bs = linspace(-pi/2.0,pi/2.0,nb,endpoint=False)

        # Create vectors of non-stress and non-b dependent values:
        mp_lons, mp_lats = self.seg_midpoints()
        w_length = self.seg_lengths()/self.length

        # number of segments in this feature
        nsegs = len(w_length)

        # Create vectors containing all of the lons and lats at which stresses
        # will need to be calculated, based on the locations of the vertices
        # making up the lineament and the values of b at which we are doing
        # calculations
        calc_thetas = tile((pi/2)-mp_lats, nb)
        calc_phis = repeat(self.bs,nsegs) + tile(mp_lons,nb)

        # use SatStress to perform the stress calculations at those locations
        tens_mag, tens_az, comp_mag, comp_az = stresscalc.principal_components(calc_thetas, calc_phis, 0.0)

        # Create an (nsegs*nb) length array of w_stress values
        w_stress = (tens_mag - comp_mag)/stresscalc.mean_global_stressdiff()

        self.nsrstresswts = ((tile(w_length,nb)*w_stress)).reshape(nb,nsegs).sum(axis=1)

        doppels = self.doppelgen_midpoint_nsr(stresscalc, init_doppel_res=init_doppel_res, doppel_res=doppel_res, num_subsegs=num_subsegs)

        # Shift all of the doppelgangers such that they are superimposed upon
        # the prototype lineament:
        for doppel,b in zip(doppels,self.bs):
            doppel.lons -= b

        # for each point in self (the prototype) find the minimum distance to
        # any point in each doppelganger - this measures the similarity in
        # their shape.  It's normalized by self.length so as to be unitless,
        # and scaled linearly to the overall length of the feature.
        d_min = ravel(array([ self.d_min(doppel) for doppel in doppels ]))/self.length

        # note that this is the RMS minimum distance...
        self.nsrdbars = sqrt(((tile(w_length,nb)*(d_min**2))).reshape(nb,nsegs).sum(axis=1))

    #}}}2 end calc_nsrfits

    def nsrfits(self, dbar_max=0.125, use_stress=True): #{{{2
        """
        Return a heuristic metric of the probability that the feature fits the
        chosen stress field at each value of b, multiplied by the significance
        of such a fit, as measured by the lack of perturbability of the fit in
        this location within the stress field.

        """

        assert(self.bs is not None and self.nsrdbars is not None and self.nsrstresswts is not None)

        if use_stress is True:
            w_stress = self.nsrstresswts
        else:
            w_stress = 1.0

        return(w_stress*(1.0 - where(self.nsrdbars/dbar_max < 1.0, (self.nsrdbars/dbar_max), 1.0))**2)

    #}}}2

    def best_fit(self, dbar_max=0.125, use_stress=True): #{{{2
        """
        Return best fit, and the value of b at which it occurs.

        """

        nsrfits = self.nsrfits(use_stress=use_stress, dbar_max=dbar_max)
        best_fit = max(nsrfits)
        best_b = float(self.bs[where(fabs(nsrfits-best_fit) < 1e-9)[0][0]])

        return(best_fit, best_b)
    
    def best_dbar(self):
        """
        Return the lineament's smallest dbar, and the value of b at which it occurs

        """

        best_dbar = min(self.nsrdbars)
        best_b = float(self.bs[where(fabs(self.nsrdbars-best_dbar) < 1e-9)[0][0]])

        return(best_dbar, best_b)

    #}}}2

    def bfgc_pole(self): #{{{2
        """
        Return the (lon, lat) point on the sphere defining the great circle
        which minimizes the sum of the squares of the distance between it and
        the lineament.

        """

        weights = self.seg_lengths()/self.length
        mp_lons, mp_lats = self.seg_midpoints()

        # special case for features consisting of a single segment... in which
        # case all the eigenstuff is both broken and unnecessary:
        if len(mp_lons) == 1:
            best_pole_phi, bp_lat = spherical_reckon(mp_lons[0], mp_lats[0], spherical_azimuth(self.lons[0], self.lats[0], self.lons[-1], self.lats[-1])+pi/2.0, pi/2.0)
            best_pole_theta = pi/2.0 - bp_lat
        else:
            mp_x, mp_y, mp_z = sphere2xyz(weights, pi/2 - mp_lats, mp_lons)
            A = dot(array([mp_x, mp_y, mp_z]),array([mp_x, mp_y, mp_z]).transpose())
            eigenvals, eigenvecs = eig(A)
            best_pole_x, best_pole_y, best_pole_z = eigenvecs[:,where(fabs(eigenvals - min(eigenvals)) < 1e-9)[0][0] ]
            best_pole_r, best_pole_theta, best_pole_phi = xyz2sphere(best_pole_x, best_pole_y, best_pole_z)

        return(best_pole_phi, pi/2.0 - best_pole_theta)
    #}}}2 end bfgc_pole

    def plot(self, map, lon_cyc=2*pi, lat_mirror=False, color='black', alpha=1.0, linewidth=1.0): #{{{2
        """
        Plot the lineament on the provided Basemap object, map, with the
        associated color, alpha, and width.

        if lon_cyclic is non-zero, modulo the longitudes of the feature by the
        cyclic value, until any that can ever fall in the displayed area of the
        map do so, plotting the feature once for each distinct time it is the
        case.  The default value, 2*pi, means that features which cross the
        edge of the map will have that portion that runs out of the display
        plotted on the far side.  If lon_cyclic=pi, then all the features will
        be placed in the same non-repeating pi-wide portion of the stress
        field, for easy comparison therein.

        If lat_mirror is True, use the absolute values of the latitudes,
        instead of the latitudes, in order to map the lineament into the
        northern hemisphere.  This allows more compact comparison of lineaments
        in stress-equivalent locations.
        
        """

        plotted_lines = []

        # What does it mean to have a longitude discontinuity?  It means
        # that two adjacent points in the lineament have longitude values
        # which differ by more than they need to... that is, for two points
        # ptA, ptB, having longitudes lonA and lonB, there exists some integer
        # N, such that:
        # abs(lonA - (2*pi*N + lonB)) < abs(lonA - lonB) 
        # this can only be true if:
        # abs(lonA - lonB) > pi
        lons = fixlons(self.lons)

        # In order to create the illusion that the lineaments are wrapping
        # around in longitude, we plot one copy at each longitude where a
        # portion of them will show up on the map...
        nrange = unique1d(floor((lons-radians(map.llcrnrlon))/lon_cyc))
        for N in nrange:
            (x,y) = map(degrees(lons-N*lon_cyc), degrees(self.lats))
            plotted_lines.append(map.plot(x, y, color=color, alpha=alpha, linewidth=linewidth))

            # If we're mirroring across the equator, we need to plot all of the
            # lineaments with their anti-latitudes as well:
            if lat_mirror:
                (x,y) = map(degrees(lons-N*lon_cyc), -degrees(self.lats))
                plotted_lines.append(map.plot(x, y, color=color, alpha=alpha, linewidth=linewidth))

        return(plotted_lines)
    #}}}2

#}}}1 end of the Lineament class

################################################################################
# Helpers having to do with fit metrics or lineament generation.
################################################################################
def lingen_nsr(stresscalc, init_lon=None, init_lat=None, max_length=2*pi, prop_dir="both", seg_len=0.01, num_subsegs=10): # {{{
    """
    Generate a synthetic NSR feature, given a starting location, maximum
    length, propagation direction, and a step size on the surface.

    Assumes tensile fracture, perpendictular to the most tensile principal
    component of the stresses.

    """


    # Make a note of the fact that we've got to do the other half too
    if prop_dir == "both":
        max_length = max_length/2.0
        prop_dir = "east"
        done = False
    else:
        done = True

    lons = array([init_lon,])
    lats = array([init_lat,])
    # Calculate the stresses at the given time and initial location
    (tens_mag, tens_az, comp_mag, comp_az) = stresscalc.principal_components(theta=(pi/2.0)-lats[0], phi=lons[0], t=0.0)

    lin_length = 0.0
    ice_strength = stresscalc.stresses[0].satellite.layers[-1].tensile_str

    # This can't be vectorized because we don't know where we'll end up until
    # we get there.
    while lin_length < max_length and tens_mag > ice_strength:
        # Choose which direction to go based on prop_dir, knowing that comp_az
        # is always an angle clockwise from north, between 0 and pi
        if prop_dir == "east":
            prop_az = comp_az
        else:
            assert(prop_dir == "west")
            prop_az = comp_az + pi

        # in order to ensure that the distance between the points representing the feature
        # is unlikely to be the determining factor in how well it compares to another
        # feature, the great-circle segment making up the feature is broken down into many
        # pieces.  This is done by passing in a list of distances to spherical reckon, instead
        # of just a single distance, with the linspace(0,seg_len,num_subsegs).
        newlons, newlats = spherical_reckon(lons[-1], lats[-1], prop_az, linspace(0,seg_len,num_subsegs))

        # Make sure that our new longitude is within 2*pi in longitude of the
        # previous point in the feature, i.e. don't allow any big
        # discontinuities (this is a hack to deal with longitude cyclicity)
        while (abs(newlons[-1] - lons[-1]) > abs((newlons[-1] - 2*pi) - lons[-1])):
            newlons = newlons - 2*pi
        while (abs(newlons[-1] - lons[-1]) > abs((newlons[-1] + 2*pi) - lons[-1])):
            newlons = newlons + 2*pi
        
        lons = concatenate([lons,newlons])
        lats = concatenate([lats,newlats])

        lin_length += seg_len

        # Calculate the stresses at the new location
        (tens_mag, tens_az, comp_mag, comp_az) = stresscalc.principal_components(theta=(pi/2.0)-lats[-1], phi=lons[-1], t=0.0)

    # if we only got a single point, then we failed to initiate a fracture, and
    # should not even try to make the second part.
    if len(lons) > 1:
        # because seg_len may be a significant portion of the length of the
        # feature, it's easy to end up with something that's a bit too long,
        # and because we're interpolating vertices between the actual stress
        # calculations, it's easy to trim the feature down to be as close to
        # the target length as possible:
        first_part = Lineament(lons=lons, lats=lats, stresscalc=stresscalc)
        length_to_trim = first_part.length - max_length

        if length_to_trim > 0:
            subseg_len = seg_len/num_subsegs
            nv2t = np.int(length_to_trim/subseg_len)
            if nv2t > 0:
                lons = first_part.lons[:-nv2t]
                lats = first_part.lats[:-nv2t]

        if not done:
            second_part = lingen_nsr(stresscalc, init_lon=init_lon, init_lat=init_lat, max_length=max_length, prop_dir="west", seg_len=seg_len, num_subsegs=num_subsegs)
            # the second_part lat/lon arrays are in reversed order here to preserve
            # the overall directionality of the vertices (see Numpy fancy slicing).
            # also, we don't want to include the midpoint twice, so we only use up through the
            # next-to-last point in the reversed second part.
            lons = concatenate([second_part.lons[:0:-1], lons])
            lats = concatenate([second_part.lats[:0:-1], lats])

    return(Lineament(lons=lons, lats=lats, stresscalc=stresscalc))

#}}} end lingen_nsr

def lingen_nsr_library(nlats=36): #{{{
    """
    Create a regularaly spaced "grid" of synthetic NSR lineaments, for use in
    visualizing the expected NSR trajectories, and in fast fitting of
    lineaments to the stress field.

    """
    import satstress
    satfile   = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = satstress.Satellite(open(satfile,'r'))
    NSR = satstress.StressCalc([satstress.NSR(europa),])

    linlist = []
    for lat in linspace(0,pi/2,nlats+2)[1:-1]:
        linlist.append(lingen_nsr(NSR, init_lon=0.0, init_lat=lat, max_length=2*pi, prop_dir='both'))

    # Now mirror and shift that set of regular lineaments to cover the surface entirely:
    linlist += [Lineament(lons=lin.lons, lats=-lin.lats, stresscalc=lin.stresscalc) for lin in linlist]
    linlist += [Lineament(lons=lin.lons+pi, lats=lin.lats, stresscalc=lin.stresscalc) for lin in linlist]

    return linlist
#}}}

def lingen_greatcircle(init_lon, init_lat, fin_lon, fin_lat, seg_len=0.01): #{{{
    """
    Return a L{Lineament} object closely approximating the shortest great
    circle route between (lon1,lat1) and (lon2,lat2).

    """

    gc_length  = spherical_distance(init_lon, init_lat, fin_lon, fin_lat)
    # this forces each feature to have at least 10 segments
    num_segs = max(10, gc_length/seg_len)
    seg_dists = linspace(0,gc_length,num_segs)
    init_az    = spherical_azimuth(init_lon, init_lat, fin_lon, fin_lat)
    lons, lats = spherical_reckon(init_lon, init_lat, init_az, seg_dists)
    return(Lineament(lons=lons, lats=lats))

#}}} end lingen_greatcircle

def d_min(linA, linB): #{{{
    """
    For each point in linA, find the minimum distance to a point in linB.

    """
    mp_lonsA, mp_latsA = linA.seg_midpoints()
    mp_lonsB, mp_latsB = linB.seg_midpoints()

    # Create a matrix of lat/lon points such that pairwise calculations in
    # spherical_distance do the entire grid of possible distances between
    # points in linA and linB:
    distmatrix = spherical_distance(repeat(mp_lonsA,len(mp_lonsB)), repeat(mp_latsA,len(mp_latsB)),\
                                      tile(mp_lonsB,len(mp_lonsA)),   tile(mp_latsB,len(mp_latsA)))

    # reshape distmatrix such that each column corresponds to a point in linA and find the
    # minimum distances for each point in A:
    return(distmatrix.reshape((len(mp_lonsA), len(mp_lonsB))).min(axis=1))

#}}} end d_min

def mhd(linA, linB): #{{{
    """
    Calculate the MHD from linA to linB.  Duh.

    """

    return(sum(d_min(linA, linB)*(linA.seg_lengths()/linA.length))/linA.length)
#}}}

def find_nearest_lins(lins=None, lons=None, lats=None, d_max=0.01): #{{{
    """
    Given a list of lineaments (lins) and a set of points on the surface,
    defined by (lons,lats), return a list of lineaments chosen from lins, which
    are the features closest to the respective (lon,lat) points.  If no
    features within the list of lineaments are within d_max of a (lon,lat)
    point, for that point return a lineament having (lon,lat) as its only
    vertex (a zero-length lineament, indicating failed doppelganger creation).

    """
    from scipy.spatial import KDTree

    # first, we need to create a lookup table that keeps track of which
    # lineament from the library a given lat/lon point is in, so that when we
    # figure out where that nearest point is, we can choose the appropriate
    # library lneament as a doppelganger.

    lib_lons = concatenate([ lin.lons for lin in lins ])
    lib_lats = concatenate([ lin.lats for lin in lins ])
    lib_Ns    = concatenate([ repeat(N, len(lin.lons)) for N,lin in zip(range(len(lins)),lins) ])

    dtype = [('lons',float),('lats',float),('Ns',int)]
    lins_lookup = zeros(len(lib_lons), dtype=dtype)
    lins_lookup['lons'] = lib_lons
    lins_lookup['lats'] = lib_lats
    lins_lookup['Ns']   = lib_Ns

    # Now we create a KDTree to search, using all the lon/lat points making up
    # the library of features, but converted to cartesian x,y,z:
    radius = 1.0
    lib_x, lib_y, lib_z = sphere2xyz(radius, pi/2-lib_lats, lib_lons)
    lib_kdt = KDTree(array([lib_x, lib_y, lib_z]).T)
    x, y, z = sphere2xyz(radius, pi/2.0-lats, lons)
    dists, near_idx = lib_kdt.query(array([x, y, z]).T, distance_upper_bound=d_max)

    # For those input points which we successfully found nearest lon,lat
    # points, we return the feature which that nearest point was in, otherwise,
    # we return a zero length feature having as its only vertex the input
    # point:
    nearest_lins = []
    for idx,N in zip(near_idx,range(len(lons))):
        if idx < len(lib_lons):
            nearest_lins.append(lins[lins_lookup['Ns'][idx]])
        else:
            nearest_lins.append(Lineament(lons=array([lons[N],]), lats=array([lats[N],])))

    return(nearest_lins)
#}}}

def mhd_by_lat(init_lat, init_lon, stresscalc, seg_len, num_subsegs, lin, max_length, lonshift): #{{{
    """
    A function to be used with optimizers, for finding the best latitude at
    which to initiate a synthetic NSR feature, given a longitude (init_lon),
    the lineament it is meant to replicate (lin), the stress field to use
    (stresscalc) the number of segments the doppelgangers should have (nsegs),
    the maximum length for the synthetic features (max_length) and how far away
    in longitude from the lineament's true (mapped) location the synthesis is
    happening so that they can be shifted to the same longitude for shape
    comparison (lonshift)

    """

    # generate a doppelganger at the given location
    doppel = lingen_nsr(stresscalc=stresscalc, init_lon=init_lon+lonshift, init_lat=init_lat, seg_len=seg_len, num_subsegs=num_subsegs, max_length=max_length)
    doppel.lons -= lonshift

    # calculate and return the MHD from the prototype to it
    return(mhd(lin, doppel))

#}}} end mhd_by_lat()

def best_nsr_init_points(lin, stresscalc=None, seg_len=0.01, num_subsegs=10): #{{{
    """
    Given a lineament and a stresscalc object, find the best latitude at which
    to initiate tensile cracking when generating NSR doppelgangers.  Find one
    latitude for each value of lin.bs, and return both the lons and the lats
    found.

    Also returns max_length, which is the target length for the doppelgangers,
    based on the length of the best fit great circle segment representing lin.

    """
    from scipy.optimize import brent

    if stresscalc is None:
        stresscalc=lin.stresscalc

    # Figure out where to initiate the formation of the doppelganger:
    ep1_lon, ep1_lat, ep2_lon, ep2_lat = lin.bfgcseg_endpoints()
    mp_lon, mp_lat = spherical_midpoint(ep1_lon, ep1_lat, ep2_lon, ep2_lat)

    # Set the max_length of the doppelganger to be the length of the best
    # fit great circle, and not the prototype, since NSR features are almost
    # perfectly straight.
    max_length = spherical_distance(ep1_lon, ep1_lat, ep2_lon, ep2_lat)

    # The simplest thing to do now, in choosing the initiation point is
    # just to use (mp_lon, mp_lat), but that point won't always be very
    # close to the mapped feature.  The separation between the initiation
    # point and the mapped lineament ends up being a strong determinant of
    # how good of a fit can be attained ultimately for those features which
    # fit NSR well.  This isn't really acceptable.  One way to get around
    # this would be to use the longitude of the bfgc, and try several
    # different latitudes, choosing the one which generates the best
    # doppelganger.
    
    # Now we use the Brent scalar function optimizer to search to find the
    # right latitude for each of those longitudes.
    init_lats = []
    for b in lin.bs:
        init_lats.append(brent(mhd_by_lat, args=(mp_lon, stresscalc, seg_len, num_subsegs, lin, max_length, b), full_output=1)[0])

    return(lin.bs+mp_lon, init_lats, max_length)
#}}} end best_nsr_init_points

################################################################################
# Helpers having to do with spherical geometry
################################################################################
def spherical_distance(lon1, lat1, lon2, lat2): #{{{
    """
    Calculate the distance between two points on a sphere, in radians.

    """

    sph_len = 0.0
    sph_len += sin((lat1-lat2)/2)**2 
    sph_len += cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2)**2)
    sph_len =  2.0*arcsin(sqrt(sph_len))
    return(sph_len)
# }}} end spherical_distance

def spherical_azimuth(lon1, lat1, lon2, lat2): # {{{
    """
    Calculate the azimuth between one point and another.

    The returned azimuth angle is the initial heading (angle clockwise from
    north) to set out on along a great-circle route in order to travel from
    point 1 to point 2.

    """
    lon1 = 2.0*pi - lon1
    lon2 = 2.0*pi - lon2

    return(mod(arctan2(sin(lon1-lon2)*cos(lat2), cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon1-lon2)), 2.0*pi))

# }}} end spherical_azimuth

def spherical_reckon(lon1, lat1, az12, ang_dist): # {{{
    """
    Calculate current location given a starting point, heading, and distance
    traveled.

    """
    # This formula uses WEST positive longitude... so we need to adjust.
    lon1 = 2.0*pi - lon1

    lat2 = arcsin( sin(lat1)*cos(ang_dist) + cos(lat1)*sin(ang_dist)*cos(az12) )
    dlon = arctan2( sin(az12)*sin(ang_dist)*cos(lat1), cos(ang_dist) - sin(lat1)*sin(lat2) )
    lon2 = mod( lon1 - dlon + pi, 2.0*pi ) - pi

    # This formula used WEST positive longitude... so we need to re-adjust.
    lon2 = 2.0*pi - lon2
    return(array([lon2,lat2]))

# }}} end spherical_reckon

def spherical_midpoint(lon1, lat1, lon2, lat2): #{{{
    """
    Given two points on the surface of a sphere, return the point that lies
    halfway between them on a great circle route.

    """

    ang_dist = spherical_distance(lon1, lat1, lon2, lat2)
    az12     = spherical_azimuth(lon1, lat1, lon2, lat2)

    return(spherical_reckon(lon1, lat1, az12, ang_dist/2.0))

# }}}

def random_lonlatpoints(N): #{{{
    """
    Generate N evenly distributed random points on a sphere.

    """

    return(2*pi*random(N), (pi/2)-arccos(2*random(N)-1))
#}}}

def sphere2xyz(r, theta, phi): #{{{
    """
    Takes a point in spherical coordinates and returns its cartesian
    equivalent.

    r is the distance from the origin.

    theta is south positive co-latitude (zero at the north pole)

    phi is the east positive longitude.

    Assumes that the z-axis is the rotation axis.

    """

    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)

    return(x, y, z)
#}}}

def xyz2sphere(x, y, z): #{{{
    """
    Takes a Cartesian (x,y,z) point and returns the equivalent point in terms
    of spherical (r, theta, phi), where:

    r is the radial distance from the origin.

    theta is co-latitude.

    phi is east-positive longitude.

    And the z-axis is taken to be the same as the rotation axis.

    """

    r     = sqrt(x**2 + y**2 + z**2)
    theta = arctan2(sqrt(x**2 + y**2), z)
    phi   = arctan2(y, x)

    return(r, theta, phi)
#}}}

def paleopole_transform(pnp_lon, pnp_lat, lon_in, lat_in): #{{{
    """
    Transforms the location of a point on the surface of a sphere, defined by
    (lon_in,lat_in) in east-positive longitude, returning the (lon,lat)
    location that point would be located at if the point defined by
    (pnp_lon,pnp_lat) is moved directly north until it is at the north pole.

    """

    # First we convert the point to be transformed into Cartesian coords.
    # Since we only care about the lon,lat position in the end, we'll treat the
    # the body as a unit sphere.  Remember that sphere2xyz needs CO-latitude:
    colat_in = pi/2 - lat_in
    xyz_in = array(sphere2xyz(1.0, colat_in, lon_in))

    # Now, remember that what we're doing is bringing a wayward pole back to
    # the top of the sphere... which means the angles we're interested in are
    # actually -colat, -lon:
    alpha = pi/2 + pnp_lon 
    beta  = pi/2 - pnp_lat 
    gamma = 0

    # The definition of the X-Z rotation matrix:
    rot_mat = array([ [ cos(alpha)*cos(gamma)-sin(alpha)*cos(beta)*sin(gamma), -cos(alpha)*sin(gamma)-sin(alpha)*cos(beta)*cos(gamma),  sin(beta)*sin(alpha) ],\
                      [ sin(alpha)*cos(gamma)+cos(alpha)*cos(beta)*sin(gamma), -sin(alpha)*sin(gamma)+cos(alpha)*cos(beta)*cos(gamma), -sin(beta)*cos(alpha) ],\
                      [                 sin(beta)*sin(gamma),                                     sin(beta)*cos(gamma),                       cos(beta)      ] ])

    # Do the matrix multiplication using dot():
    xyz_out = dot(xyz_in.transpose(), rot_mat)

    # Transform back to spherical coordinates:
    r_out, theta_out, phi_out = xyz2sphere(xyz_out[:,0], xyz_out[:,1], xyz_out[:,2])

    # and back to lon/lat from colon/lat:
    lon_out = mod(phi_out + alpha,2*pi)
    lat_out = pi/2 - theta_out # co-lat to lat

    return(lon_out, lat_out)
#}}}

def fixlons(lons): #{{{
    """
    Takes a set of longitudes, and forces it to be within a continuous range
    (i.e. no skips at "international date line")

    """

    lons = mod(lons, 2*pi)
    while len(where(abs(lons[:-1]-lons[1:]) > pi)[0] > 0):
        lons[1:] = where(abs(lons[:-1]-lons[1:]) > pi, lons[1:]+2*pi*sign(lons[:-1]-lons[1:]), lons[1:])
    return(lons)
#}}}

def crop_circle(lins, clon=0, clat=0, maxdist=None): #{{{
    """
    Take a set of lineaments, lins, and remove any vertices which are more than
    maxdist radians away from the point clon, clat.  Recalculate all the
    lineament lengths, and return the resulting cropped features.  If none of
    a feature's points are within range, it is removed entirely.

    This may do strange things if you have high sinuosity features, as it does
    not split the features into separate entities if they leave and re-enter
    the area being considered.  Should be fine in general for NSR lineaments
    and great circle segments.

    """
    if maxdist is None:
        cropped_lins = lins
    else:
        cropped_lins = []
        newlins = update_lins(lins)
        for lin in newlins:
            keep_idx = np.where(spherical_distance(clon, clat, lin.lons, lin.lats) < maxdist)
            if len(keep_idx[0]) > 1:
                lin.lons = lin.lons[keep_idx]
                lin.lats = lin.lats[keep_idx]
                lin.length = lin.calc_length()
                cropped_lins.append(lin)

    return(cropped_lins)
            
#}}}

################################################################################
# Helpers for input/output/saving/updating: (TODO: portable output...)
################################################################################
def update_lins(lins): # {{{
    """
    Just a quick and easy way to remember what all has to be done in order to
    update a set of lineaments to include all the newest and bestest
    functionality from the module.

    """

    newlins = []
    linlen=len(lins)
    for lin,N in zip(lins,range(linlen)):
        newlin = Lineament(lons=lin.lons, lats=lin.lats, stresscalc=lin.stresscalc, bs=lin.bs, nsrdbars=lin.nsrdbars, nsrstresswts=lin.nsrstresswts)
        newlins.append(newlin)

    return(newlins)
#}}} end update_lins

def plotlinmap(lins, map=None, color='black', alpha=1.0, linewidth=1.0, lon_cyc=2*pi, lat_mirror=False): #{{{
    """
    Plot a map of the lineaments listed in 'lins'.  Plot it to 'map' if given.

    """
    if map is None:
        thefig = figure()
        map = Basemap()
        map.drawmeridians(range(-180,181,30), labels=[1,0,0,1])
        map.drawparallels(range(-90,91,30), labels=[1,0,0,1])

    wasinteractive = isinteractive()
    interactive(False)
    lines = []

    lines += [ lin.plot(map, color=color, alpha=alpha, linewidth=linewidth, lon_cyc=lon_cyc, lat_mirror=lat_mirror) for lin in lins ]

    interactive(wasinteractive)

    return([line for line in flatten(lines)], map)
#}}} end plotlinmap

def shp2lins(shapefile, stresscalc=None): #{{{
    """
    Create a list of L{Lineament} objects from an ESRI shapefile.

    The shapefile must contain one or more linear features, defined in
    geographic coordinates using decimal degrees.  The shapefile is read in
    using the GDAL/OGR geospatial library.

    If there is a field named 'gid' in the attribute table of the features,
    then each Lineament object's 'gid' value is set using it.  This is useful
    for pulling features back into GIS afterward, and having them be
    identifiable.

    """

    # This is the list of lineaments we are going to return eventually:
    linlist = []
    # First read in the shapefile as an OGR "data source"
    data_source = ogr.Open(shapefile, update=0)

    # OGR data sources can in general have many data layers, but ours will only have
    # one, containing linear features.  Get that layer:
    lineaments = data_source.GetLayer(0)

    # Within that layer, there will be potentially many features: individual
    # lineaments which should be extracted and turned into Lineament objects
    # independently.  Each one should be extracted and made into a Lineament object
    # independently:

    ogr_lin_feat = lineaments.GetNextFeature()

    while ogr_lin_feat is not None:
        ogr_lin_geom = ogr_lin_feat.GetGeometryRef()
        pointlist = array([ ogr_lin_geom.GetPoint(i)[0:2] for i in range(ogr_lin_geom.GetPointCount()) ])

        if(len(pointlist) > 0):
            newlats = radians(pointlist[:,1])
            newlons = radians(pointlist[:,0])
            try:
                gid = ogr_lin_feat.GetField(ogr_lin_feat.GetFieldIndex('gid'))
            except(ValueError):
                gid = None
            linlist.append(Lineament(lons=newlons, lats=newlats, stresscalc=stresscalc, gid=gid))

        ogr_lin_feat = lineaments.GetNextFeature()

    return(linlist)

# }}} end shp2lins

def lins2kml(lins=[], kmlfile=None): #{{{ TODO: WRITE IT!
    """
    Create a KML file from a list of L{Lineament} objects, defining their
    shapes.

    """
    pass
#}}} end lins2kml

def lins2shp(lins=[], shapefile=None): #{{{ TODO: WRITE IT!
    """
    Create a shapefile from a list of L{Lineament} objects, including additions
    to the attribute tables that have been calculated so that those attributes
    can be displayed easily within ArcGIS.

    """
    pass
#}}} end lins2shp

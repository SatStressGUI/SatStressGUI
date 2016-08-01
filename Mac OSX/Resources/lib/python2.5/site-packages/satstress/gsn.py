#!/usr/bin/python
"""
A Geological Superposition Network (GSN) is a data structure encoding both
geographic and stratigraphic information about features on the surface of a
planetary body.  It uses directed graphs (or "digraphs"), a construct from
discrete mathematics, to manipulate both the geography and superposition
relationships by storing each geologic feature as a "node" (or "vertex") in the
graph, and each cross-cutting relationship between features as an "edge", such
that if a feature A is seen to cross-cut another feature B, there exists a
directed edge in the graph from B to A, and in the language of graph theory, A
is said to be a successor to B.  The direction of the edge thus indicates the
arrow of time, with the implication that because A appears superimposed upon B,
B must be older than A.

Within each node of the GSN is stored the geographic information pertaining to
the geologic feature it represents.  In the applications below, the features
being analyzed are linear tectonic features, and so satstress Lineament objects
are used as nodes.

It is possible for two features to have more than one intersection, and so it
follows that in some instances there may be multiple edges between the same two
nodes.  Because of this, it is necessary for each edge to have additional data
associated with it to allow multiple edges between the same two nodes to be
distinguished from one another.  It is also useful to have the location of each
intersection stored, for display and processing purposes.  Thus, each edge also
has a latitude and longitude value.

Because in some instances the exact nature of the superposition relationship is
unclear, the ordering of each intersection is assigned an estimated probability
of its being correct.  When this value is 1, the ordering is completely clear.
When it is 0.5, no temporal information is discernible.

Because the constraints that the cross cutting relationships place on the
overall order of formation are in general fairly loose, it is usually not
practical to use a GSN to try and infer a particular chronology.  However, the
construct is still useful in answering certain kinds of questions, such as:

"What are the chances that one feature formed before or after another?"

which is really just a special case of:

"To what degree are the mapped stratigraphic relationships compatible with a
particular ordering inferred by some other method?"

A GSN may also be used to identify sets of stratigraphic relationships which
appear to be logically incompatible, and therefore indicate non-instantaneous
activity.  Such sets of intersections form loops, or "cycles" within the
directed graph.

"""

import networkx as nx
import pygraphviz as pgv
import pickle
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pylab
from mpl_toolkits import basemap
from osgeo import ogr
from osgeo import osr
from . import lineament
from . import satstress
from . import nsrhist

# define the coordinate system we're working in:
# This is lame, but for the moment, let's assume that we're on Europa:
IAU2000_Europa_SRS_WKT="""GEOGCS["Europa 2000",DATUM["D_Europa_2000",SPHEROID["Europa_2000_IAU_IAG",1564130.0,488.79062499999998]],PRIMEM["Greenwich",0],UNIT["Radian",1.000000000]]"""

e15_lin_shp="""/Users/zane/Desktop/GSN_Mapping/E15_LinsGSN_byGID"""
e15_cross_shp="""/Users/zane/Desktop/GSN_Mapping/E15_LinCrossClean"""

class GeoSupNet(nx.MultiDiGraph): #{{{
    """
    A Geological Superposition Network (GSN) encapsulates the cross-cutting
    relationships between a number of geologic features, as well as their
    geography.  It provides methods for resolving the stratigraphic
    relationships between the features, even when some of those relationships
    are ambiguous.  The class is particularly designed to deal with linear
    tectonic features of the type found on Jupiter's moon Europa and other
    geologically active icy bodies.

    """

    # In the event that we're actually generating this GSN from an ordered list
    # of features, we should store that information for later reference.
    linstack = []

    ##########################################################################
    # Translation:
    ##########################################################################

    def to_bidirected(self): #{{{2
        """
        Take a unidirectional GSN, in which all intersection confidences are
        greater than or equal to 0.5, and return the equivalent bidirected GSN,
        in which every edge with confidence P from A to B has a corresponding
        edge with confidence 1-P from B to A.

        """
        bidirected_gsn = GeoSupNet()
        forward_edges = self.edges(data=True, keys=True)
        back_edges = []
        for e in forward_edges:
            top_lin = e[0]
            bottom_lin = e[1]
            key = e[2]
            lon = e[3]['lon']
            lat = e[3]['lat']
            weight = -np.log(1.0-np.exp(-e[3]['weight']))
            edge_data = {'weight':weight,'lon':lon,'lat':lat}
            back_edges.append((bottom_lin,top_lin,key,edge_data))

        bidirected_gsn.add_edges_from(forward_edges)
        bidirected_gsn.add_edges_from(back_edges)
        return(bidirected_gsn)
    #}}}2

    def to_shp(lin_shp_file, cross_shp_file): #{{{2 # TODO: write it!
        """
        Output shapefiles representing the GSN, for display in a GIS
        application.

        """
        pass
    #}}}2

    def get_sub_GSNs(self, minconf=0.5, minsize=2): #{{{2
        """
        Return a list of all the connected components of the GSN having size
        greater than or equal to minsize, when only intersections with weights
        corresponding to confidences greater than minconf. If the original GSN
        has a linstack associated with it, add the subset of lineaments in each
        sub_GSN to its associated linstack in the same order in which they
        appear in the original GSN.

        """

        # Convert the GSN to an undirected graph
        # and find the connected subcomponents
        good_edges = [ e for e in self.edges(data=True, keys=True) if np.exp(-e[3]['weight']) >= minconf ]
        U = nx.MultiGraph()
        U.add_edges_from(good_edges)
        connected_nbunches = nx.connected_components(U)
        G = GeoSupNet()
        G.add_edges_from(good_edges)

        # re-constitute those connected bunches as a new list of GSNs:
        sub_GSNs = []
        for nbunch in connected_nbunches:
            if len(nbunch) >= minsize:
                sub_gsn = G.subgraph(nbunch)
                # Add the appropriate sub-linstack if it's available:
                sub_gsn.linstack = [ ]
                for lins in self.linstack:
                    sub_gsn.linstack.append( [ lin for lin in lins if lin in sub_gsn ] )
                sub_GSNs.append(sub_gsn)

        return(sub_GSNs)
    #}}}2 end get_sub_GSNs()

    ##########################################################################
    # GSN Manipulation
    ##########################################################################

    def stratigraphic_sort(self): #{{{2 # TODO: bi-directed, log-weights
        """
        Find upper and lower bounds on the stratigraphic location of each
        feature in the GSN.

        This is done by finding the sizes of the successor and predecessor
        digraphs.  Anything that isn't in one of those graphs is of
        indistinguishable age.

        """

        stratsorts = []
        for sub_gsn in self.get_sub_GSNs():
            stratsort = {}
            rev_gsn = sub_gsn.reverse()
            for lin in sub_gsn.nodes():
                ub = len(sub_gsn) - len(nx.single_source_shortest_path_length(sub_gsn, lin))
                lb = len(nx.single_source_shortest_path_length(rev_gsn, lin)) - 1
                stratsort[lin] = (lb, ub)
            stratsorts.append(stratsort)
        
        return(stratsorts)

    #}}}2 end stratigraphic_sort

    def load_fits(self, linfile): #{{{2
        """
        Insert pre-calculated NSR fits into the lineaments that are GSN nodes.

        """
        for fitlin in nsrhist.load_lins(linfile):
            for gsnlin in self.nodes():
                if(fitlin == gsnlin):
                    gsnlin.stresscalc = fitlin.stresscalc
                    gsnlin.bs = fitlin.bs
                    gsnlin.nsrdbars = fitlin.nsrdbars
                    gsnlin.nsrstresswts = fitlin.nsrstresswts
    #}}}2

    ##########################################################################
    # Metrics and Analysis:
    ##########################################################################

    def shortest_paths(self, minconf=0.5): #{{{2
        """
        For each ordered pair of nodes (A,B) in the GSN, find the shortest
        weighted path from A to B.  Return both the paths, and their
        confidences, as two dictionaries of dictionaries such that for example:

        conf[A][B]

        would be the confidence of the path from A to B, which are part of the
        largest sub-GSN, and:

        path[A][B]

        would be the corresponding path, represented as an ordered list of
        nodes.

        Paths with confidences less than minconf are excluded.

        """

        path = {}
        dist = {}
        conf = {}
        maxdist = -np.log(minconf)

        # Find the shortest path from every node to every other reachable node.
        # Store the results (both the distances and the paths) in two
        # dictionaries of dictionaries, keyed by source and target node:
        for source in self.nodes():
            dist[source], path[source] = nx.single_source_dijkstra(self, source, cutoff=maxdist)
            conf[source] = {}

        # Transform the "distance" dictionary into a "confidence" dictionary,
        # since that's the metric we're mostly using.
        for source in dist.keys():
            for target in dist[source].keys():
                conf[source][target] = np.exp(-dist[source][target])

        return(conf, path)
    #}}}2

    def enumerate_cycles(self, minconf=0.5): #{{{2
        """
        Find all pairs of nodes (A,B) for which there exist paths both from A
        to B and B to A with confidences greater than minconf.  Many of these
        paths will be isomorphic with each other (i.e. passing through all of
        the same nodes... but with different starting and ending points A and
        B).  Compare all of the paths, and return one from each isomorphic
        group, as a list of nodes.  Also return the confidences of the cycles.

        (confs, cycles)

        """

        cycles = []
        cycle_confs = []
        conf, path = self.shortest_paths(minconf=np.sqrt(minconf))

        for n1 in self.nodes():
            for n2 in self.nodes():
                try:
                    if conf[n1][n2] and conf[n2][n1] and n1 != n2:
                        cycles.append(path[n1][n2][:-1]+path[n2][n1][:-1])
                        cycle_confs.append(conf[n1][n2]*conf[n2][n1])

                except(KeyError):
                    continue

        # only include one copy of each isomorphic cycle
        unique_cycles = []
        unique_cycle_confs = []
        cycles_as_sets = []
        for i in range(len(cycles)):
            cyc_set = set(cycles[i])
            if cyc_set not in cycles_as_sets:
                cycles_as_sets.append(cyc_set)
                unique_cycles.append(cycles[i])
                unique_cycle_confs.append(cycle_confs[i])

        # sort the returned cycles first by confidence, and then by number of
        # nodes involved in the cycle, so that the most important, and easily 
        # understood cycles come first.
        dtype = [('cycle',object),('path_length',int),('conf',float)]
        cycle_sort = np.zeros(len(unique_cycles), dtype=dtype)
        cycle_sort['cycle'] = unique_cycles
        cycle_sort['path_length'] = np.array([ len(c) for c in unique_cycles ])
        cycle_sort['conf'] = -np.array(unique_cycle_confs)
        cycle_sort.sort(order=['conf','path_length'])

        return(-cycle_sort['conf'], cycle_sort['cycle'])

    #}}}2 end enumerate_cycles

    def net_relations(self, minconf=0.0, minsize=2): #{{{2
        """
        In a GSN with N features, here are N*(N-1) possible ordered binary
        relationships which could be defined.  Relationships which have
        confidences less than unity we consider to only be partially defined.
        If both the ordered relations (A,B) and (B,A) are defined, then there
        is ambiguity as to the actual ordering of the features, and the net
        definition is: |P(A,B)-P(B,A)|.  The sum of those terms, for all
        possible pairs (A,B) is what this function returns.
        
        """

        info_sums = []
        for sub_gsn in self.get_sub_GSNs(minconf=minconf, minsize=minsize):
            conf, path = sub_gsn.shortest_paths(minconf=minconf)

            info_sum = 0.0
            for source in conf.keys():
                for target in conf[source].keys():
                    try:
                        forward = conf[source][target]
                    except(KeyError):
                        forward = 0.0
                    try:
                        backward = conf[target][source]
                    except(KeyError):
                        backward = 0.0

                    info_sum += np.fabs(forward-backward)

            # need to divide the sum by two, because we'll catch each one twice
            # iterating over the whole set of nodes as above:
            info_sums.append(info_sum/2.0)
        
        return(info_sums)
    #}}}2

    def completeness(self, minconf=0.0, minsize=2): #{{{2
        """
        The completeness of a GSN is an indication of how many of the possible
        pairwise stratigraphic relationships between the mapped features are
        actually defined by their cross cutting relationships.  If both the
        relations (A,B) and (B,A) are defined, then the stratigraphic
        relationship between A and B has actually been made ambiguous, and so
        the magnitude of the difference between them is used as a metric.  I.e.
        if P(A,B)=0.9 and P(B,A)=0.1 then the net definition is 0.8, and that's
        the contribution to the overall completeness that comes from the pair
        (A,B).  The total number of relationships which could be defined is
        N choose 2 = N*(N-1)/2 if there are N features in the GSN.
        
        """

        sub_GSNs = self.get_sub_GSNs(minconf=minconf, minsize=minsize)
        gsn_lens = np.array([ len(sub_gsn) for sub_gsn in sub_GSNs ])

        net_rels = np.array(self.net_relations(minconf=minconf,minsize=minsize))

        return( net_rels/(gsn_lens*(gsn_lens-1.0)/2.0) )
    #}}}2

    def intersection_confidences(self): #{{{2
        """
        Return an array containing the confidence values for all of the
        intersections in the GSN.

        """
        weights = np.array([e[2]['weight'] for e in self.edges(data=True)])
        return(np.exp(-weights))
    #}}}2

    def shortest_path_confidences(self, minconf=0.0): #{{{2
        """
        Return an array containing the confidence values for all of the
        shortest paths calculated within the GSN.

        """
        confs, paths = self.shortest_paths(minconf=minconf)
        out_confs = []
        for source in self.nodes():
            for target in confs[source].keys():
                out_confs.append(confs[source][target])

        return(out_confs)
    #}}}2

    def agreement(self, linstack, minconf_path=0.0): #{{{2
        """
        Calculate a metric of how consistent the relationships defined by
        linstack are with the relationships defined by the GSN.

             agree
        ----------------
        agree + disagree

        where agree is the sum of the confidences of the paths from A to B for
        all paths (A,B) defined by both the GSN and the linstack, and disagree
        is the sum of the confidences of all the paths (A,B) defined in the GSN
        where the path (B,A) is defined in the linstack.

        """

        my_confs, my_paths = self.shortest_paths(minconf=minconf_path)

        linstack_pairs = linstack2pairs(linstack)

        agree = 0.0
        disagree = 0.0
        for n1 in my_confs.keys():
            for n2 in my_confs[n1].keys():
                if (n1,n2) in linstack_pairs:
                    agree += my_confs[n1][n2]
                if (n2,n1) in linstack_pairs:
                    disagree += my_confs[n1][n2]

        if agree+disagree == 0.0:
            return(1.0)
        else:
            return(agree/(agree+disagree))
    #}}}2

    def disagreement(self, confs): #{{{2
        """
        The complement to "agreement", a metric of how much in conflict two
        sets of binary relations are.

        """
        return(1-self.agreement(confs))
    #}}}2

    def reorder(linstack): #{{{2 TODO
        """
        Adjust the direction and confidences of the edges in the GSN to be
        consistent with the ordering of the features in linstack.

        - In order to be able to tell how remarkable a particular ordering is, i.e.
          whether it's agreement() is better than we'd expect from chance, we need
          to be able to know how well chance does.
        - Need to be able to re-order a linstack quickly, so we can calculate many
          different values of agreement() and find out what the distribution is like.
        - Doing it quickly means that we need not to re-calculate all of the locations
          of the intersections, we need to be able to just switch the direction of the
          intersections, if need be.
        - Several different kinds of random reordering are possible:
          a) Keep the same layers, but re-shuffle their ordering.
          b) Keep the same number of layers, and the same number of features in each
             one, but shuffle which features are in which layer.
          c) Keep the same number of layers, but randomly assign features to each one
          d) Get rid of layers altogether, and choose a complete ordering.
        """
    #}}}2

    def geometric_information_distribution(gsn, iterations=100, minsize=2): #{{{2 TODO
        """
        Calculate the relative information content of a GSN repeatedly, re-ordering
        its true order of formation (linstack) repeatedly, while holding its
        geometry constant.  Do this for each connected component separately, and
        return 2D array of the results, with one row for each sub GSN.

        """

        # TODO: totally broken: doesn't work with weighted edges, no longer have a
        # fast reordering function.  May or may not re-write.  Saved here for
        # reference.

        # Make a copy of the GSN so we don't alter the original:
        test_gsn = gsn2gsn(gsn)

        sub_GSNs = test_gsn.get_sub_GSNs(minsize=minsize, minconf=minconf)
        sub_comps = np.zeros((len(sub_GSNs),iterations))
        for sub_gsn,n in zip(sub_GSNs,range(len(sub_GSNs))):
            for i in range(iterations):
                pylab.shuffle(sub_gsn.linstack)
                sub_gsn.reorder(sub_gsn.linstack)
                sub_comps[n,i] = sub_gsn.completeness()[0]

        return(sub_comps)
    #}}}2

    def significance(self, linstack, num_iter=100, minconf_cross=0.7, minconf_path=0.0): #{{{2
        """
        Calculate the agreement between linstack and the GSN, using only those
        edges in the GSN with confidence greater than minconf_cross, and those
        paths with confidence greater than minconf_path.

        Then try to figure out whether or not that level of agreement is
        actually statistically significant, by calculating the agreement
        between the GSN and a suite of num_iter randomly ordered linstacks.

        Return a tuple (hyp_A, rand_As) where hyp_A is the agreement between
        the linstack passed in (the hypothesis) and rand_As is a list of the
        agreements between the randomized linstacks and the GSN.  These can be
        used to make a histogram, or calculate the percentile score of the
        hypothesis.

        """

        test_gsn = gsn2gsn(self, minconf=minconf_cross)
        print("Calculating agreement between hypothesis and GSN:")
        hyp_A = test_gsn.agreement(linstack, minconf_path=minconf_path)
        print("    hyp_A = %g" % (hyp_A,) )

        flat_linstack = []
        for lins in linstack:
            for lin in lins:
                flat_linstack.append([lin,])
        rand_A = []
        print("Calculating agreement between random linstacks and GSN:")
        while len(rand_A) < num_iter:
            pylab.shuffle(flat_linstack)
            rand_A.append(test_gsn.agreement(flat_linstack, minconf_path=minconf_path))
            print("    A = %g (%d / %d) " % (rand_A[-1], len(rand_A), num_iter) )

        return(hyp_A, rand_A)
    #}}}2

    def nsr_compare(self, nsr_stresscalc=None, nb=180, nbins=9, num_iter=1000): #{{{2
        """
        Calculate NSR fits for the features in the GSN, and using those that do fit
        NSR, create a linstack (having nbins layers) based on their best fit
        backrotation.  Calculate the agreement between that ordering and the
        ordering relationships implied by the GSN.  Also calculate the agreement
        between num_iters different random orderings of the features and the GSN
        for comparison.

        """
        lins = self.nodes()
        linstack = nsr_linstack(lins, nb=nb, nbins=nbins, nsr_stresscalc=nsr_stresscalc)

        hyp_A, rand_A = self.significance(linstack, num_iter=num_iter)

        return(hyp_A, rand_A)

    #}}}2

    ##########################################################################
    # Display Routines:
    ##########################################################################

    def draw_graph(self, ax=None, minsize=2, minconf=0.5, cmap=pylab.cm.jet_r): #{{{2
        """
        Draw graph view of the GSN.  Plot each sub-GSN in a separate figure.
        This is janky and should be replace with something that will plot using
        graphviz instead of matplotlib.

        minsize allows you to limit which connected components are plotted, and
        avoid showing a bunch of uninteresting size 2 and 3 graphs.

        """
        from matplotlib.colors import rgb2hex
        from matplotlib import rcParams
        import Image

        sub_GSNs = self.get_sub_GSNs(minsize=minsize, minconf=minconf)

        C = self.completeness(minconf=minconf)
        for sub_gsn,n in zip(sub_GSNs, range(len(sub_GSNs))):

            # generate an array of colors to use on the nodes:
            nc_array = np.array([ cmap(float(n)/len(sub_GSNs)) for i in range(len(sub_gsn)) ])

            sub_agsn = nx.to_agraph(sub_gsn)
            sub_agsn.graph_attr.update(rankdir  = 'BT',\
                                       nodesep  = '0.25',\
                                       ranksep  = '0.5',\
                                       rank     = 'source')
            sub_agsn.node_attr.update(color=rgb2hex(nc_array[0][0:3]), shape='point', width='0.25')
            sub_agsn.layout(prog='dot')
            sub_agsn.draw('output/graphs/test_agsn.png')

            png_out = Image.open('output/graphs/test_agsn.png')
            dpi = rcParams['figure.dpi']
            fig_width  = png_out.size[0]/dpi
            fig_height = png_out.size[1]/dpi
            the_fig = plt.figure(figsize=(fig_width, fig_height))
            gax = the_fig.add_subplot(111, frameon=False)
            im = plt.imshow(png_out, origin='lower')
            gax.set_xticks([])
            gax.set_yticks([])
            gax.set_xlabel('C=%.3g' % (C[n],))

                
    #}}}2 end draw_graph

    def draw_map(self, intersections=False, ax=None, cmap=pylab.cm.jet_r, minsize=2, minconf=0.5): #{{{2
        """
        Draw a map view of the GSN using Basemap and Matplotlib.

        """
        if ax is None:
            the_fig = plt.figure(figsize=(12,6))
            ax = the_fig.add_subplot(1,1,1)
            map_ax = basemap.Basemap(ax=ax)
            map_ax.drawmeridians(range(-180,181,30), labels=[1,0,0,1])
            map_ax.drawparallels(range(-90,91,30), labels=[1,0,0,1])

        sub_GSNs = self.get_sub_GSNs(minsize=minsize, minconf=minconf)

        for sub_gsn,n in zip(sub_GSNs, range(len(sub_GSNs))):
            lines, gsn_map = lineament.plotlinmap(sub_gsn.nodes(), map=map_ax, linewidth=2,\
                                                  color=cmap(float(n)/len(sub_GSNs)))
            if intersections is True:
                edges = sub_gsn.edges(data=True)
                edge_lons = [ edge[2]['lon'] for edge in edges ]
                edge_lats = [ edge[2]['lat'] for edge in edges ]
                gsn_map.scatter(np.degrees(edge_lons), np.degrees(edge_lats), lw=0, color='black', s=50)
                gsn_map.scatter(np.degrees(edge_lons)+360.0, np.degrees(edge_lats), lw=0, color='black', s=50)
                gsn_map.scatter(np.degrees(edge_lons)-360.0, np.degrees(edge_lats), lw=0, color='black', s=50)

        plt.draw()

    #}}}2 end draw_map

    def draw_sort(self, orderby='mean', ax=None, colorby='graph', cmap=pylab.cm.jet_r, title="", weighted=True, normalized=True, minsize=2): #{{{2 TODO: Borken!
        """
        Make a plot showing the retrieved stratigraphic sort.
        
        The features may be ordered by any of the following:
    
             'nid' -- net in degree (default)
            'true' -- the actual ordering, as implied by trueorder
            'mean' -- mean stratigraphic location
              'lb' -- lower bound (earliest possible formation time)
              'ub' -- upper bound (latest possible formation time)
    
        If the true ordering is indicated by the linstacks embedded with the
        GSN, then it will be displayed.

        If any of weighted or normalized are set, they are passed on to the
        GeoSupNet.net_in_degree() function (with orderby='nid')

        minsize determines how many features have to be the connected subgraph
        in order for it to be displayed.  This allows you to avoid showing a
        bunch of sorts that only include 2 or 3 features.

        """
        if ax is None:
            the_fig = plt.figure(figsize=(12,6))
            ax = the_fig.add_axes((0.05,0.1,0.9,0.8))

        sub_GSNs = self.get_sub_GSNs(minsize=minsize, minconf=minconf)
        N_tot = np.sum([len(sub_gsn) for sub_gsn in sub_GSNs ])

        x0 = 0
        for sub_gsn, n_col in zip(sub_GSNs, range(len(sub_GSNs))):
            sub_sort  = sub_gsn.stratigraphic_sort()
            assert (len(sub_sort) == 1)
            sub_sort = sub_sort[0] # should have only one connected component...
            sub_stack = sub_gsn.linstack

            # Figure out what order we're going to plot the features in, and
            # generate the appropriate x, y, and (yerr_lo, yerr_hi) arrays for
            # plotting:
            if orderby == 'mean':
                sort_key = [ 0.5*(sub_sort[lin][0]+sub_sort[lin][1]) for lin in sub_sort.keys() ]
                lins_to_plot = sub_sort.keys()
            elif orderby == 'lb':
                sort_key = [ sub_sort[lin][0] for lin in sub_sort.keys() ]
                lins_to_plot = sub_sort.keys()
            elif orderby == 'ub':
                sort_key = [ sub_sort[lin][1] for lin in sub_sort.keys() ]
                lins_to_plot = sub_sort.keys()
            elif orderby == 'nid':
                nids = self.net_in_degree(normalized=normalized, weighted=weighted, with_labels=True)
                sort_key = [ nids[node] for node in sub_sort.keys() ]
                lins_to_plot = sub_sort.keys()
            elif orderby == 'true':
                assert(len(sub_stack) > 0)
                sort_key = []
                lins_to_plot = []
                for lin in sub_stack:
                    if lin in sub_sort.keys():
                        sort_key.append(sub_stack.index(lin))
                        lins_to_plot.append(lin)
            else: # Uh, we should never get here
                print("Bad orderby string found in draw_stratsort")

            # Create a structured array so we can sort by whatever key we want:
            if len(sub_stack) > 0:
                dtype = [('lin',object),('sort_key',float),('trueorder',int)]
            else:
                dtype = [('lin',object),('sort_key',float)]

            arr_to_sort = np.zeros(len(lins_to_plot), dtype=dtype)
            arr_to_sort['lin'] = lins_to_plot
            arr_to_sort['sort_key'] = sort_key
            if len(sub_stack) > 0:
                arr_to_sort['trueorder'] = [ sub_stack.index(lin) for lin in lins_to_plot ]
            arr_to_sort.sort(order='sort_key')

            symb_size = 810.0/N_tot
            X = np.arange(len(arr_to_sort))
            lower_bounds = np.array([ sub_sort[lin][0] for lin in arr_to_sort['lin'] ])-0.5
            upper_bounds = np.array([ sub_sort[lin][1] for lin in arr_to_sort['lin'] ])+0.5
            ax.vlines(X+x0, lower_bounds, upper_bounds, colors=cmap(float(n_col)/len(sub_GSNs)), linewidth=symb_size)

            if len(sub_stack) > 0:
                ax.plot(X+x0, arr_to_sort['trueorder'], 'k_', markeredgewidth=2, markersize=0.85*symb_size, color='k', lw=0)

            x0 = x0 + len(sub_gsn)

        ax.set_xlim(-1,N_tot)
        ax.set_xlabel('N')
        ax.set_ylim(-1,len(sub_GSNs[0]))
        ax.set_ylabel('N')
        sortnames = {  'ub':'upper stratigraphic bound',\
                       'lb':'lower stratigraphic bound',\
                     'mean':'mean stratigraphic location',\
                     'true':'a priori formation time',\
                      'nid':'net in degree'}

        if title != "":
            title = " of "+title
        ax.set_title("StratSort"+title+" (ordered by "+sortnames[orderby]+")")
        ax.grid(True)

    #}}}2 end draw_sort

    def draw_idist(self, iterations=100, minsize=4, cmap=pylab.cm.jet_r): #{{{2 TODO: Borken!
        """
        Calculate some statistics of the possible information content within the
        GSN, depending on what order the features formed in.

        """
        sub_GSNs = self.get_sub_GSNs(minsize=minsize, minconf=minconf) 
        GIDs = geometric_information_distribution(self, iterations=iterations, minsize=minsize)
        nc_array = np.array([ cmap(float(i)/len(GIDs)) for i in range(len(GIDs)) ])
        the_fig = plt.figure(figsize=(12,6))
        hist_ax = the_fig.add_subplot(111)
        for idist,sub_gsn,nc in zip(GIDs, sub_GSNs, nc_array):
            n, bins, patches = hist_ax.hist(idist, facecolor=nc, normed=True, histtype='stepfilled', lw=2, edgecolor=nc, alpha=0.5)
            patches[0].set_label("%d features" % (len(sub_gsn),) )

        hist_ax.set_xlim((0,1))
        hist_ax.grid(True)
        hist_ax.legend(loc='upper right')
        hist_ax.set_title('%d samples' % (iterations,) )
        hist_ax.set_xlabel('I (proportion of recoverable information)')
        hist_ax.set_ylabel('N (normalized)')

        return(GIDs)
    #}}}2
#}}} end GeoSupNet

################################################################################
# GSN Generation:
################################################################################

# There are three different representations of a GSN:
#
#   1.) The canonical graph representation, based on the nx.MultiDiGraph object
#
#   2.) A set of two shapefiles, one defining the geometries of the lineaments
#       and another defining their crosscutting relationships at those points
#       where they intersect.
#
#   3.) A "linstack", which is an ordered list of lists of Lineament objects,
#       with all the members of each set having indistinguishable stratigraphic
#       locations.  The list is ordered chronologically, such that the first
#       set is the oldest features, and the last set is the youngest.  A
#       linstack thus defines a known a priori order of formation, from which
#       the crosscutting relationships at their intersections can be inferred.
#
#       A linstack is the structure used both to create a synthetic GSN, and to
#       phrase a hypothesized order of formation to the GSN for corroboration
#       or denial.  When creating a synthetic GSN from a linstack, any
#       intersection involving two features from the same stratigraphic set
#       (the ordering of which we have no information about) is assigned a
#       confidence of 0.5.  When posing a hypothesized ordering, 
#
# Additionally, there is a data structure which is derivable from a GSN, which
# describes all of the stratigraphic relationships it implies, and their
# confidences, based on the shortest (highest confidence) path between all
# connected ordered pairs of features in the graph.  This is the GSN's
# path_conf, and it is represented as a dictionary of dictionaries, indexed
# first by source node (older feature) and second by target node (younger
# feature), with the value stored being the confidence of the shortest path
# from the source to the target.

def gsn2gsn(old_gsn, minconf=0.0): #{{{
    """
    Create a new GSN from an old one.  Only include those intersections whose
    confidence is greater than or equal to minconf (useful for pruning out
    intersections deemed too ambiguous for inclusion)

    """
    new_gsn = GeoSupNet()
    old_edges = old_gsn.edges(data=True, keys=True)
    good_old_edges = [ e for e in old_edges if np.exp(-e[3]['weight']) >= minconf ]
    new_gsn.add_edges_from(good_old_edges)
    for lins in old_gsn.linstack:
        new_gsn.linstack.append([ lin for lin in lins if lin in new_gsn ])

    return(new_gsn)
#}}}

def shp2gsn(lin_shp, cross_shp): #{{{
    """
    Given the paths to two shapefiles, one defining a set of mapped lineaments
    (lin_shp) and one describing the superposition relationships at their
    intersections (cross_shp) return a GSN.

    Each feature within lin_shp must have a field in its attribute table named
    'gid', which contains an integer uniquely identifying it within the
    layer.

    Each feature within cross_shp must have three fields in its attribute
    table, named 'top_lin', 'bottom_lin', and 'conf', indicating the gid of the
    features within lin_shp which participate in the intersection, which one
    appears to be on top, and which on the bottom, and the confidence with
    which that assignment is made, as a subjective probability, whose value
    ranges between 0.5 and 1.0, with conf=0.5 indicating that the imaging
    from which the features was digitized did not allow any determination of
    the relative ages of the features at that point.  Probabilities less than
    0.5 are prohibited, because they would in effect be indicating that the
    ordering of features ought to be reversed.

    The confidence is transformed to -log(conf) in the graph weighting to allow
    us to use the sums of weights to infer the products of probabilities.

    """

    # this is the GSN that we're going to construct and return
    the_gsn = GeoSupNet()

    # Create a dictionary of all the lineaments defined in lin_shp, keyed by
    # their gid values:
    lins = lineament.shp2lins(lin_shp)
    # need to be able to look up a lineament object based on its gid:
    lindict = {}
    for lin in lins:
        lindict[lin.gid] = lin
    # check and make sure that all of the gid values were unique:
    assert(len(lindict) == len(lins))

    # Add all of the lineament object to the GSN, just in case any of them don't
    # have any intersections, so they don't get left out:
    the_gsn.add_nodes_from(lins)

    # Iterate over the set of points described in cross_shp, and for each one
    # look up the lineament objects involved in the intersection in the
    # dictionary we just defined, and associate those Lineament objects with the
    # edge data we can read from the points
    cross_pt_data = ogr.Open(cross_shp, update=0)
    cross_pt_layer = cross_pt_data.GetLayer(0)
    cross_pt = cross_pt_layer.GetNextFeature()
    while cross_pt is not None:
        bottom_lin_idx = cross_pt.GetFieldIndex('bottom_lin')
        bottom_lin_id  = cross_pt.GetField(bottom_lin_idx)
        bottom_lin     = lindict[bottom_lin_id]

        top_lin_idx = cross_pt.GetFieldIndex('top_lin')
        top_lin_id  = cross_pt.GetField(top_lin_idx)
        top_lin     = lindict[top_lin_id]

        # The confidence is the probability (0.5 < P < 1.0) that the 
        # lineaments have the observed temporal ordering.  P < 0.5 is
        # prohibited because it would be equivalent to P=(1-P), with
        # the observed temporal ordering reversed.
        conf_idx = cross_pt.GetFieldIndex('conf')
        conf     = cross_pt.GetField(conf_idx)
        # Because most weighted graph algorithms sum the weights of the
        # edges, and what we're really interested in is their product
        # (the probability of a particular path being true) we need to
        # transform the confidence into a weight, logarithmically:
        pt_wt = -np.log(conf)

        # read in the lon/lat location of the feature:
        cross_pt_geom = cross_pt.GetGeometryRef()
        pt_lon, pt_lat, no_z_coord = cross_pt_geom.GetPoint()
        pt_lon = np.mod(np.radians(pt_lon),2.0*np.pi)
        pt_lat = np.radians(pt_lat)

        # add the edge:
        the_gsn.add_edge(bottom_lin, top_lin, key=lonlat2key(pt_lon,pt_lat),\
                         lon=pt_lon, lat=pt_lat, weight=pt_wt)

        # If conf=0.5, no superposition relationship is defined, and we add
        # edges in both directions:
        if np.fabs(conf-0.5) < 1e-6:
            the_gsn.add_edge(top_lin, bottom_lin, key=lonlat2key(pt_lon,pt_lat),\
                             lon=pt_lon, lat=pt_lat, weight=pt_wt)

        # grab the next intersection:
        cross_pt = cross_pt_layer.GetNextFeature()

    return(the_gsn)
#}}}

def before(linstack, n): #{{{
    """
    Return a flat list of all the features in the layers prior to layer n

    """
    assert(n < len(linstack))
    assert(n >= 0)

    before_lins = []
    if n == 0:
        pass
    else:
        before_lins = [ j for i in linstack[:n] for j in i ]

    return(before_lins)
#}}}

def after(linstack, n): #{{{
    """
    Return a flat list of all the features in the layers following layer n.

    """
    assert(n < len(linstack))
    assert(n >= -1)
    after_lins = []
    if n == len(linstack)-1:
        pass
    else:
        after_lins = [ j for i in linstack[n+1:] for j in i ]

    return(after_lins)
#}}}

def linstack2pairs(linstack): #{{{
    """
    Create a list of pairs of lineaments, (A,B) where B comes after A in the
    linstack.

    """
    pairs = set([])
    for n in range(len(linstack)-1):
        for A in linstack[n]:
            for B in [ j for i in linstack[n+1:] for j in i ]:
                pairs.add((A,B))
    return(pairs)

#}}}

def linstack2gsn(linstack, conf_dist=None): #{{{
    """
    A linstack is a list of lists of Lineament objects.  The ordering of the
    list of lists indicates the order of formation, with the first list being
    the oldest features, and the last list being the newest.  Each sub-list
    describes a set of features having indistinguishable stratigraphic
    locations.  This allows one to specify a hypothesized order of formation
    of classes of features, without having to artificially specify the 
    temporal relationships between the features within one particular class.

    If conf_dist is None, any intersection between two features in different
    classes will be given a confidence of 1.0.  Intersections between two
    features in the same class are given a confidence of 0.5 (indicating that
    there is no information contained in the intersection), and edges in both
    directions will be added (as there is no reason to prefer one direction
    over the other when the confidence is 0.5)

    If conf_dist is not None, it is taken to be a pool of confidences to draw
    from, when assigning a confidence to each intersection in the network.
    This allows the creation of ambiguous, imperfect networks for testing, or
    for accurately modeling the dynamics of a given mapped network.
    
    """

    # Initialize the GSN
    the_gsn = GeoSupNet()
    # clear out any empty slots in the linstack:
    linstack = [ lins for lins in linstack if len(lins) > 0 ]
    # Save the linstack for later reference:
    the_gsn.linstack = linstack

    # Add all of the lineaments in the stack to the GSN as nodes.
    [ the_gsn.add_nodes_from(lins) for lins in linstack ]

    # For each layer within the linstack we want to do two things.  First, any
    # intersections between features in the same layer of the linstack should
    # be added with confidences of 0.5.  Then, any intersections between two
    # features in different linstack layers should be added with confidence 1.0
    for n in range(len(linstack)):
        new_edges = []
        # if we haven't been given a distribution of confidences to use, assume
        # that we're building an idealized network, in which we get all of the
        # confidences correct (based on the ordering of the linstack):
        if conf_dist is None:
            new_edges += build_gsn_edges(bottom_lins = linstack[n],\
                                         top_lins    = linstack[n],\
                                         conf_dist   = np.array([0.5,]))

            new_edges += build_gsn_edges(bottom_lins = linstack[n],\
                                          top_lins    = after(linstack, n),\
                                          conf_dist   = np.array([1.0,]))

        # if we've been given a distribution of confidences (conf_dist) then
        # we're trying to recreate the reliability of an existing network, or
        # synthesize an imperfectly understood set of relationships:
        else:
            new_edges += build_gsn_edges(bottom_lins = linstack[n],\
                                         top_lins    = after(linstack, n-1),\
                                         conf_dist   = conf_dist)
            
        the_gsn.add_edges_from(new_edges)

    return(the_gsn)

#}}} end linstack2gsn

def build_gsn_edges(bottom_lins, top_lins, conf_dist): #{{{
    """
    Create and return a list of edge tuples representing intersections in the
    GSN between bottom_lins and top_lins.  Assign those intersections a
    confidence of conf.

    In order to do this all of the intersections between the features (which
    are implied by their geometry), have to be calculated.  This is done with
    the GDAL/OGR library.  The geometries of the lineaments are converted to
    OGRGeometry objects, and the OGRGeometry.intersection() method is used.

    The coordinate system associated with the lineaments is the IAU 2000
    geographic coordinate system for Europa.  If you're using another body, you
    should find the appropriate definition string here:

    http://spatialreference.org/ref/iau2000/

    """

    # See beginning of file for definition of the WKT SRS
    Europa_SRS = osr.SpatialReference(IAU2000_Europa_SRS_WKT)

    edges = []
    for bottom_lin in bottom_lins:
        bottom_ogr = ogr.CreateGeometryFromWkt(bottom_lin.wkt(), Europa_SRS)
        for top_lin in top_lins:
            if top_lin == bottom_lin:
                continue
            
            # We have to iterate over a few different ranges of 2*pi to make
            # sure that we get all the intersections, because, unfortunately,
            # OGR doesn't know the world is round.
            for lonshift in 2*np.pi*np.array([-1,0,1]):
                top_ogr = ogr.CreateGeometryFromWkt(top_lin.lonshift(lonshift).wkt(), Europa_SRS)
                lin_cross = bottom_ogr.Intersection(top_ogr)

                if lin_cross is not None:
                    for point in [ lin_cross.GetPoint(i) for i in range(lin_cross.GetPointCount()) ]:
                        lon = np.mod(point[0],2*np.pi)
                        lat = point[1]
                        # choose a random confidence from the proffered distribution:
                        weight = -np.log(conf_dist[np.random.randint(len(conf_dist))])
                        key = lonlat2key(lon, lat)
                        edges.append( (bottom_lin, top_lin, key, {'lon':lon, 'lat':lat, 'weight':weight}) )

    return(edges)

#}}} end build_gsn_edges

################################################################################
# Linstack generation:
################################################################################

def linstack_random(nlins=100, nbins=10, maxlen=1.25): #{{{
    """
    Create a random population of great circle segments, and transform it into a
    binned linstack for testing, having nbins separate peer classes.  To get a
    linstack in which each feature has a unique stratigraphic position, make
    nbins >> nlins.

    """

    linstack = [ [] for i in range(nbins) ]
    lins = nsrhist.random_gclins(nlins, maxlen=maxlen)
    for lin,n in zip(lins,np.random.randint(nbins,size=nlins)):
        lin.lons = lineament.fixlons(lin.lons)
        linstack[n].append(lin)

    return(linstack)
#}}}

def linstack_nsr(nlins=100, nbins=18, length_scale=1.0): #{{{
    """
    Create a set of synthetic NSR lineaments at a variety of values of
    backrotation, and transform them into a linstack having nbins separate peer
    classes.  To get a linstack in which each feature has a unique peer class,
    make nbins >> nlins.

    """

    # read in what we need to calculate the NSR stresses:
    satfile = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = satstress.Satellite(open(satfile,'r'))
    nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])

    # create some random NSR features
    lins = []
    while len(lins) < nlins:
        lon, lat = lineament.random_lonlatpoints(1)
        length = np.abs(np.random.normal(scale=length_scale))
        newlin = lineament.lingen_nsr(nsr_stresscalc, init_lon=lon[0], init_lat=lat[0], max_length=length)
        lins.append(newlin)

    bs = np.random.uniform(low=-np.pi, high=0.0, size=nlins)

    linstack = [ [] for i in range(nbins) ]
    for lin, b in zip(lins, bs):
        lin.lons = lin.lons - b
        n = int(nbins*b/np.pi)
        linstack[n].append(lin)

    return(linstack)

#}}} end linstack_nsr

def linstack_regular(nlins=10, overpole=False): #{{{
    """
    Create a regular map with a known order of formation for testing the GSN
    algorithms. nlins indicates how many features there are in each of the
    north-south and east-west orientations.

    """

    northsouth_init_lons = np.linspace(-np.pi/3.1, np.pi/3.1, nlins)
    northsouth_fin_lons  = np.linspace(-np.pi/3.1, np.pi/3.1, nlins)
    northsouth_init_lats = np.array([-np.pi/2.1,]*nlins)
    northsouth_fin_lats  = np.array([ np.pi/2.1,]*nlins)

    eastwest_init_lons = np.array([-np.pi/3.0]*nlins)
    eastwest_fin_lons  = np.array([ np.pi/3.0]*nlins)
    eastwest_lats      = np.linspace(-np.pi/3.0, np.pi/3.0, nlins)

    if overpole is True:
        northsouth_fin_lons += np.pi
        northsouth_fin_lats = np.array([ np.pi/2.01,]*nlins)

    init_lons = np.concatenate([northsouth_init_lons, eastwest_init_lons])
    init_lats = np.concatenate([northsouth_init_lats, eastwest_lats])
    fin_lons  = np.concatenate([northsouth_fin_lons,  eastwest_fin_lons])
    fin_lats  = np.concatenate([northsouth_fin_lats,  eastwest_lats])

    outlins = []
    for lon1, lat1, lon2, lat2 in zip(init_lons, init_lats, fin_lons, fin_lats):
        outlins.append(lineament.lingen_greatcircle(lon1, lat1, lon2, lat2))

    pylab.shuffle(outlins)

    return([ [lin,] for lin in outlins ])

#}}} end linstack_regular 

def linstack_file(nlins=0, linfile='output/lins/map_nsrfit', spin=False, tpw=False, nbins=10): #{{{
    """
    Draw nlins lineaments from a saved pool of features, and create a linstack
    from them.

    """
    # read in the features, and update them to make sure we've got whatever
    # interesting data they have, in the new Lineament object form.
    lins = lineament.update_lins(nsrhist.load_lins(linfile))

    if nlins == 0:
        nlins = len(lins)

    # If we're not randomizing, we can't resample with replacement, or we'll
    # get repeats that overlap themselves.
    if spin is False and tpw is False:
        # choose nlins random features from the dataset:
        # this will fail if they've asked for more features than we can get without randomizing.
        assert(nlins <= len(lins))
        pylab.shuffle(lins)
    else:
        lins = nsrhist.make_crazy(nsrhist.linresample_byN(lins), tpw=tpw, spin=spin)

    linstack = [ [] for i in range(nbins) ]
    for lin,n in zip(lins[:nlins],np.random.randint(nbins,size=nlins)):
        linstack[n].append(lin)

    return(linstack)
#}}}

def nsr_linstack(lins, nsr_stresscalc=None, nb=180, nbins=180): #{{{
    """
    Take a list of lineaments, and turn them into a linstack ordered by
    the value of their best fit backrotation.

    """
    if nsr_stresscalc is None:
        print("Initializing satellite")
        satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
        europa = satstress.Satellite(open(satfile,'r'))
        nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])

    print("Calculating NSR fits:")
    linstack = [ [] for i in range(nbins) ]
    # we only need to calculate the fits if they haven't been done already
    for lin,n in zip(lins,range(len(lins))):
        if lin.bs is None or len(lin.bs) != nb:
            if np.mod(n+1,10) == 0:
                print("    %d / %d" % (n+1,len(lins)))
            lin.calc_nsrfits(nb=nb, stresscalc=nsr_stresscalc,\
                             init_doppel_res=0.0, doppel_res=0.1,\
                             num_subsegs=10)

        best_fit, best_b = lin.best_fit()
        if best_fit > 0.0:
            n = -int(nbins*best_b/np.pi)
            linstack[n].append(lin)

    linstack = [ lins for lins in linstack if len(lins) > 0 ]

    return(linstack)

#}}}

def linstack_swap_element(linstack, nswaps=1, by_layer=False): #{{{
    """
    Take two randomly selected elements in the linstack and swap them.  Perform
    this operation nswaps times.

    If by_layer is True, swap two entire layers within the linstack.  If it is
    False, choose two individual features and swap them.

    """
    tmp_linstack = [ [ lin for lin in lins ] for lins in linstack ]
    for swaps in range(nswaps):
        if by_layer is True:
            # pick two random layers and swap them:
            i1, i2 = np.random.randint(len(tmp_linstack), size=2)
            tmp_1 = [ lin for lin in tmp_linstack[i1] ]
            tmp_2 = [ lin for lin in tmp_linstack[i2] ]
            tmp_linstack[i1] = tmp_2
            tmp_linstack[i2] = tmp_1

        else:
            # pick two random features and swap them:
            lay_idx_1, lay_idx_2 = np.random.randint(len(tmp_linstack), size=2)
            lin_idx_1 = np.random.randint(len(tmp_linstack[lay_idx_1]))
            lin_idx_2 = np.random.randint(len(tmp_linstack[lay_idx_2]))
            lin_1 = tmp_linstack[lay_idx_1][lin_idx_1]
            lin_2 = tmp_linstack[lay_idx_2][lin_idx_2]
            tmp_linstack[lay_idx_1][lin_idx_1] = lin_2
            tmp_linstack[lay_idx_2][lin_idx_2] = lin_1

    return(tmp_linstack)
#}}}

################################################################################
# Helper Functions:
################################################################################
def wkt2proj(wkt_str=None): #{{{
    """
    Take a "Well Known Text" (WKT) spatial reference system (SRS) and convert
    it to a set of Proj.4 command line arguments, for consumption by qGIS.

    """
    if wkt_str is None:
        wkt_str = """PROJCS["Europa_Mercator_AUTO",GEOGCS["Europa 2000",DATUM["D_Europa_2000",SPHEROID["Europa_2000_IAU_IAG",1564130.0,488.79062499999998]],PRIMEM["Greenwich",0],UNIT["Decimal_Degree",0.0174532925199433]],PROJECTION["Mercator_1SP"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Central_Meridian",180],PARAMETER["Standard_Parallel_1",0],UNIT["Meter",1]]"""
    coord_ref_sys = osr.SpatialReference(wkt_str)
    return(coord_ref_sys.ExportToProj4())
#}}}

def lonlat2key(lon, lat): #{{{
    """
    Take a longitude and a latitude, in radians, and turn them into a usable
    key for identifying an intersection.  Returns a 2-tuple of integers.

    """
    lon = np.mod(lon, 2*np.pi)
    lon_key = np.int(1e4*lon)
    lat_key = np.int(1e4*lat)

    return(lon_key, lat_key)
#}}}

################################################################################
# Testing Functions:
################################################################################

def test(lintype='nsr', nlins=300, nbins=9, minconf=0.5, maxconf=1.0,\
         clon=0, clat=0, maxdist=np.pi/4.0, maxlen=1.0, minsize=8, iterations=100,\
         draw_map=True, draw_graph=True, draw_idist=False, draw_sort=False): #{{{
   """
   Create a set of synthetic NSR features, binned according to the amount of
   shell rotation at which they formed.  Transform them into a GSN, calculate
   some metrics, and display the data in several ways, to see if things are
   working, generally.
   
   """

   print("Generating %s linstack (N=%d)" % (lintype,nlins) )
   if lintype == 'nsr':
       linstack = linstack_nsr(nlins=nlins, nbins=nbins, length_scale=maxlen)
   else:
       assert(lintype == 'random')
       linstack = linstack_random(nlins=nlins, nbins=nbins, maxlen=maxlen)

   lines, linmap = lineament.plotlinmap([])
   colors = pylab.cm.jet(np.linspace(0,1,len(linstack)))
   for n in range(len(linstack)):
       linstack[n] = lineament.crop_circle(linstack[n], clon=clon, clat=clat, maxdist=maxdist)
       lines, linmap = lineament.plotlinmap(linstack[n], map=linmap, color=colors[n], linewidth=2.0)

   print("Converting map to GSN")
   the_gsn = linstack2gsn(linstack)

   print("Plotting results")
   if draw_map is True:
       the_gsn.draw_map(minsize=minsize, minconf=minconf)
   if draw_graph is True:
       the_gsn.draw_graph(minsize=minsize, minconf=minconf)
   if draw_sort is True:
       test_gsn.draw_sort(orderby=orderby, title=lintype, minsize=minsize)
   if draw_idist is True:
       test_gsn.draw_idist(minsize=minsize, minconf=minconf, iterations=iterations)

   return(the_gsn)

#}}} end test()

def test_agreement(nlins=300, nbins=10, conf_dist=None, layers=True, keep_order=True): #{{{
    """
    Create a linstack having nlins features, organized into nbins layers, with
    intersections conforming to the confidence distribution passed in via
    conf_dist (or having perfect confidences, if conf_dist is None).

    If layers is True, remove one or several layers at a time from the
    linstack, and re-calculate the agreement with the previously generated GSN.
    If layers is False, then remove individual lineaments randomly instead of
    whole layers.

    If keep_order is True, don't re-order the linstack, just remove things.  If
    keep_order is False, then shuffle things around instead, swapping either
    layers or lineaments (depending on whether layers is True), and see how the
    agreement degrades.  Perform iters swaps, and then return.

    """

    print("Generating synthetic linstack and GSN (nlins=%d, nbins=%d)" % (nlins,nbins) )
    the_linstack = linstack_random(nlins=nlins, nbins=nbins)

    # the density and connectivity of the network is important for these
    # results, so we're going to confine it:
    the_linstack = [ lineament.crop_circle(lins, clon=0, clat=0, maxdist=np.pi/4.0) for lins in the_linstack ]

    the_gsn = linstack2gsn(the_linstack, conf_dist=conf_dist)
    # just keep the biggest connected component:
    the_gsn = the_gsn.get_sub_GSNs()[0]
    the_linstack = the_gsn.linstack
    # see what we got:
    the_gsn.draw_map()
    plt.draw()

    results = []
    # this should always start out 1.0...
    results.append(the_gsn.agreement(the_linstack))
    print("A = %g" % (results[-1],))

    if keep_order is True:
        print("Testing agreement with removed elements...")
        # Only makes sense to continue while the linstack actually still
        # defines some relationships.
        while len(the_linstack) >= 2:
            if layers is True:
                del the_linstack[np.random.randint(len(the_linstack))]
            else:
                # randomly remove one feature at a time, and any empty layers:
                layer_idx = np.random.randint(len(the_linstack))
                lin_idx = np.random.randint(len(the_linstack[layer_idx]))
                del the_linstack[layer_idx][lin_idx]
                the_linstack = [ lins for lins in the_linstack if len(lins) > 0 ]

            results.append(the_gsn.agreement(the_linstack))
            print("A = %g" % (results[-1],))

    else: # we're shuffling things
        if layers is True:
            maxiters = len(the_linstack)
            print("Running %d tests with %d shuffled layers." % (maxiters, maxiters))
        else:
            maxiters = np.sum( [ len(lins) for lins in the_linstack ] )
            print("Running %d tests with %d shuffled lins." % (maxiters, maxiters))

        while len(results) < maxiters+1:
            if layers is True:
                the_linstack = linstack_swap_element(the_linstack, by_layer=True, nswaps=1)
            else:
                the_linstack = linstack_swap_element(the_linstack, by_layer=False, nswaps=1)
            
            results.append(the_gsn.agreement(the_linstack))
            print("A = %g" % (results[-1],))

    return(results)
#}}}

def transpose_degradation(linstack, perfect_gsn=None, nsamples=50, nswaps=10, iters=100): #{{{
    results = np.zeros((iters, nsamples))
    if perfect_gsn is None:
        print("converting linstack to GSN")
        perfect_gsn = linstack2gsn(linstack)

    for i in range(iters):
        print(" %d / %d " % (i+1,iters) )
        swapped = linstack
        for j in range(nsamples):
            results[i,j] = perfect_gsn.agreement(swapped)
            swapped = linstack_swap_element(swapped, by_layer=False, nswaps=nswaps)

    results = np.mean(results, axis=0)

    return(results)
#}}}

################################################################################
# Figures:
################################################################################

def makefigs(fast=True, cycle_samples=50, num_iters=500): #{{{
    """
    Generate all the plots required for the Geological Superposition Networks
    portion of my thesis.

    """
    if fast is True:
        cycle_samples = 10
        num_iters  = 100

    e15_gsn = shp2gsn(e15_lin_shp, e15_cross_shp)
    # keep only the biggest sub-GSN:
    e15_gsn = e15_gsn.get_sub_GSNs()[0]

    # Show the features we mapped:
    GSN_Map(e15_gsn, label='E15', min_cross_conf=0.7, K_min=0.5)
    plt.draw()

    # Create a histogram of the mapped intersection confidences:
    Intersection_Confs(e15_gsn, label='E15')
    plt.draw()

    # Show how completeness and number of paths vary with minconf:
    Path_Confs(e15_gsn, label='E15')
    plt.draw()

    # Show how number of cycles varies with minconf, K_min:
    Cycle_Confs(e15_gsn, label='E15', numsamples=cycle_samples)
    plt.draw()

    GSN_vs_NSR(e15_gsn, iters=num_iters, label='E15', linfile='output/lins/E15_nsrfit')
    plt.draw()

    ActHist_Spun(e15_gsn, iters=num_iters, label='E15', linfile='output/lins/E15_nsrfit')
    plt.draw()

#}}} end makefigs

def GSN_Map(the_gsn, label='GSN', min_cross_conf=0.7, K_min=0.0): #{{{
    """
    Plot the features making up the GSN.  Color intersections by confidence,
    and highlight any features participating in cycles in red.

    """

    high_conf_gsn = gsn2gsn(the_gsn, minconf=min_cross_conf)
    cycle_confs, cycles = high_conf_gsn.enumerate_cycles(minconf=K_min)
    path_confs = high_conf_gsn.shortest_path_confidences(minconf=K_min)

    cycle_lins = []
    for lins in cycles:
        cycle_lins += lins

    unique_cycle_lins = []
    for lin in cycle_lins:
        if lin not in unique_cycle_lins:
            unique_cycle_lins.append(lin)

    mapped_lins = the_gsn.nodes()

    map_cycle_fig = plt.figure(figsize=(12,12))
    llcrnrlon, llcrnrlat = (-110.0, -65.0)
    urcrnrlon, urcrnrlat = (-70.0,  -25.0)
    gridspace = 5

    linmap = basemap.Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat)
    linmap.drawmeridians(np.arange(llcrnrlon,urcrnrlon+1,gridspace), labels=[1,0,0,1])
    linmap.drawparallels(np.arange(llcrnrlat,urcrnrlat+1,gridspace), labels=[1,0,0,1])
    e15_map_ax = mapc_cycle_fig.axes[0]

    cycle_lines, e15_map = lineament.plotlinmap(unique_cycle_lins, linewidth=6, color='red', map=linmap)
    dag_lines, e15_map = lineament.plotlinmap(mapped_lins, map=linmap, linewidth=2, color='black')
    e15_map_ax.set_title('%d shortest paths and %d shortest path cycles found in %s map (minconf=%g, $K_{min}=$%g)'\
                          % (len(path_confs), len(cycles), label, min_cross_conf, K_min) )

    e15_map_ax.legend((dag_lines[0], cycle_lines[0]),\
                      ('all mapped features (N=%d)' % (len(mapped_lins),       ),\
                       'features in cycles (N=%d)'  % (len(unique_cycle_lins), )),\
                      loc='upper left')

    map_cycle_fig.savefig('output/figs/E15_Map.pdf')

#}}} end E15_Map

def Intersection_Confs(the_gsn, label='Mapped'): #{{{
    """
    Create a histogram of intersection confidences, showing the overall
    distribution, with only a single confidence value in each bin (0.5-1.0).

    """

    int_confs = np.sort(the_gsn.intersection_confidences())
    # because we're actually interested in the number of intersections, and not
    # the number of edges in this case, we need to chop off half of the 0.5
    # edges (they're added bi-directionally for consistency)
    half_fives = len(int_confs[np.where(int_confs == 0.5)])/2
    int_confs = int_confs[half_fives:]
    the_fig = plt.figure(figsize=(8,5))
    hist_ax = the_fig.add_subplot(1,1,1)

    hist_ax.set_ylabel('N')
    hist_ax.set_xlabel('Confidence')
    hist_ax.set_title('%s Intersection Confidences (N=%d)' % (label,len(int_confs)) )
    hist_ax.grid()

    n, bins, counts = hist_ax.hist(int_confs, bins=6, range=(0.45,1.05),\
                                   facecolor='k', rwidth=0.8)
    hist_ax.set_xticks(np.linspace(0.5,1.0,6))

    the_fig.savefig('output/figs/%s_Intersection_Confs.pdf' % (label,) )

#}}}

def Path_Confs(the_gsn, label='GSN'): #{{{
    """
    Calculate and plot both the number of shortest paths defined in the GSN and
    the GSNs completeness, as a function of the lowest confidence intersections
    we include.

    Also show how the distribution of path confidences varies as the minimum
    intersection confidence changes.

    """

    C = []
    npaths = []
    minconfs = np.linspace(0.5,1.0,6)
    path_conf_list = []
    for minconf in minconfs:
        test_gsn = gsn2gsn(the_gsn, minconf=minconf)
        path_confs = test_gsn.shortest_path_confidences()
        path_conf_list.append(path_confs)
        npaths.append(len(path_confs))
        C.append(test_gsn.completeness()[0])

    plot_fig = plt.figure(figsize=(8,8))
    ax1 = plot_fig.add_subplot(1,1,1)
    ax1.plot(minconfs, C, linewidth=3, color='black', label='C')
    ax1.set_ylabel('C (completeness)')
    ax1.set_xlabel('minimum intersection confidence')
    ax1.set_title('%s completeness and number of paths' % (label,))
    ax1.set_ybound(lower=0)
    ax1.set_xlim(0.5,1.0)
    ax1.grid()

    ax2 = ax1.twinx()
    ax2.plot(minconfs, np.array(npaths)/1000, linewidth=3, color='black',
             linestyle='dashed', label='|R(X)|')
    ax2.set_ylabel('|R(X)| (number of paths defined, in thousands)')
    ax2.set_ybound(lower=0)
    ax2.set_xlim(0.5,1.0)

    lines = ax1.lines + ax2.lines
    labels = [ line.get_label() for line in lines ]
    ax1.legend(lines, labels, 'upper right')

    plot_fig.savefig('output/figs/%s_Completeness.pdf' % (label,))

    hist_fig = plt.figure(figsize=(12,8))
    hist_ax = hist_fig.add_subplot(1,1,1)

    for minconf,conf_hist in zip(minconfs,path_conf_list):
        hist_ax.hist(conf_hist, bins=np.logspace(-1,0,110), range=(0.1,1), linewidth=0, label=str(minconf))

    hist_ax.set_xlim(0.1,1.05)
    hist_ax.set_xscale('log')
    hist_ax.grid()
    hist_ax.legend(loc='upper right')
    xtick_vals = np.linspace(0.1,1,10)
    hist_ax.set_xticks(xtick_vals)
    hist_ax.set_xticklabels([ str(x) for x in xtick_vals ])
    hist_ax.invert_xaxis()
    hist_ax.set_title('Distribution of path confidences in %s' % (label,) )
    hist_ax.set_ylabel('N')
    hist_ax.set_xlabel('Path Confidence')

    hist_fig.savefig('output/figs/%s_Path_Conf_Dist.pdf' % (label,) )
#}}}

def Cycle_Confs(the_gsn, label='GSN', numsamples=50): #{{{
    """
    Show how the number and confidence of cycles in the GSN changes as we vary
    K_min.

    """
    K_mins = np.linspace(0,1,numsamples)

    the_fig = plt.figure(figsize=(8,8))
    ax = the_fig.add_subplot(1,1,1)

    for minconf,lw,color in zip((0.6,0.7,0.8), (15,9,3), ('0.0','0.66','0.0')):
        hiconf_gsn = gsn2gsn(the_gsn, minconf=minconf)
        num_cycles = []
        for K_min in K_mins:
            cycle_confs, cycles = hiconf_gsn.enumerate_cycles(minconf=K_min)
            num_cycles.append(len(cycle_confs))

        ax.plot(K_mins, num_cycles, linewidth=lw, color=color, label=str(minconf))

    ax.grid()
    ax.legend(labelspacing=1, handletextpad=1)
    ax.set_title('Number of cycles vs. minimum cycle confidence in %s' % (label,) )
    ax.set_ylabel('Number of cycles')
    ax.set_xlabel(r'$K_{min}$ (minimum cycle confidence)')

    # don't want to wipe out our good fig on testing runs...
    if numsamples > 20:
        the_fig.savefig('output/figs/%s_Cycle_Confs.pdf' % (label,) )

#}}} end Cycle_Confs

def GSN_vs_NSR(the_gsn, iters=100, nb=180, nbins=180, linfile=None, label='E15'): #{{{
    """
    Calculate the agreement between the lineaments in the GSN, and the ordering
    inferred from their NSR fits.

    """

    # Hopefully, we've already calculated the fits for the features in the GSN,
    # and can just read them in from linfile, but if not, we need to be ready
    # to calculate them:
    
    # read in what we need to calculate the NSR stresses:
    satfile = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = satstress.Satellite(open(satfile,'r'))
    nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])

    if linfile is not None:
        the_gsn.load_fits(linfile)

    A_nsr, A_rand = the_gsn.nsr_compare(nb=nb, nbins=nbins, nsr_stresscalc=nsr_stresscalc, num_iter=iters)
    old_A_rand = pickle.load(open('output/E15_A_rand','r'))
    A_rand = np.concatenate([A_rand,old_A_rand])
    nsrhist.safe_pickle(A_rand, name='E15_A_rand', dir='output')

    the_fig = plt.figure(figsize=(8,8))
    x_min,x_max = (0.3,0.8)
    hist_ax = the_fig.add_subplot(1,1,1)
    hist_ax.hist(A_rand, facecolor='gray', linewidth=0, bins=len(A_rand)/20, range=(x_min,x_max), normed=True)
    hist_ax.grid()
    hist_ax.axvline(x=A_nsr, color='black', linewidth=5)
    mu = np.mean(A_rand)
    sigma = np.std(A_rand)
    x = np.linspace(x_min,x_max,1000)
    hist_ax.plot(x,pylab.normpdf(x,mu,sigma), color='black', linewidth=5, linestyle='dashed')
    hist_ax.set_xlim(x_min,x_max)
    hist_ax.set_xlabel('Agreement')
    hist_ax.set_ylabel('N (normalized)')

    hist_ax.set_title(r'N=%d, $\mu$=%.3g, $\sigma$=%.3g, $A_{nsr}$=%.3g (%.3g $\sigma$)' %
                      (len(A_rand), mu, sigma, A_nsr, np.fabs(A_nsr-mu)/sigma) )

    the_fig.savefig('output/figs/%s_GSN_vs_NSR.pdf' % (label,) )
#}}} end GSN_vs_NSR

def ActHist_Spun(the_gsn, iters=100, nb=180, nbins=180, label='E15', linfile='output/lins/E15_nsrfit'): #{{{
    """
    Activity history of the GSN's region, with some spun lineaments for comparison.

    """

    # read in what we need to calculate the NSR stresses:
    satfile = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = satstress.Satellite(open(satfile,'r'))
    nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])

    the_gsn.load_fits(linfile)

    spun_maps = [ nsrhist.make_crazy(nsrhist.linresample_byN(the_gsn.nodes()),tpw=False) for i in range(iters) ]
    spun_maps.append(the_gsn.nodes())

    # pseudo-maps in gray
    colors = ['gray',]*iters
    # real map in black
    colors.append('black')

    # pseudo-maps nearly transparent
    alphas = [15./iters,]*iters
    # real map opaque
    alphas.append(1.0)

    labels = [None,]*(iters-1)
    labels.append('Spun Subsamples (N=%d)' % (iters,))
    labels.append('Mapped Lineaments')

    acthist_fig = plt.figure(figsize=(12,8))
    nsrhist.activity_history(spun_maps, labels=labels, colors=colors, alphas=alphas, the_fig=acthist_fig, verbose=False)
    acthist_fig.axes[0].set_title('Apparent Lineament Activity History for E15 region')
    acthist_fig.savefig('output/figs/%s_ActHist.pdf' % (label,) )

#}}}

def Type_Intersections():
    import Image
    from matplotlib import rcParams

    k05_1_png = Image.open('input/type_intersections/k05_1.png')
    k05_2_png = Image.open('input/type_intersections/k05_2.png')
    k05_3_png = Image.open('input/type_intersections/k05_3.png')
    k06_1_png = Image.open('input/type_intersections/k06_1.png')
    k06_2_png = Image.open('input/type_intersections/k06_2.png')
    k06_3_png = Image.open('input/type_intersections/k06_3.png')
    k07_1_png = Image.open('input/type_intersections/k07_1.png')
    k07_2_png = Image.open('input/type_intersections/k07_2.png')
    k07_3_png = Image.open('input/type_intersections/k07_3.png')

    the_fig_low_k = plt.figure(figsize=(6,6))
    gax_k05_1 = the_fig_low_k.add_subplot(331)
    gax_k05_1.set_ylabel('k=0.5')
    gax_k05_2 = the_fig_low_k.add_subplot(332)
    gax_k05_2.set_title('Low Confidence Type Intersections')
    gax_k05_3 = the_fig_low_k.add_subplot(333)
    gax_k06_1 = the_fig_low_k.add_subplot(334)
    gax_k06_1.set_ylabel('k=0.6')
    gax_k06_2 = the_fig_low_k.add_subplot(335)
    gax_k06_3 = the_fig_low_k.add_subplot(336)
    gax_k07_1 = the_fig_low_k.add_subplot(337)
    gax_k07_1.set_ylabel('k=0.7')
    gax_k07_2 = the_fig_low_k.add_subplot(338)
    gax_k07_3 = the_fig_low_k.add_subplot(339)

    im1 = gax_k05_1.imshow(k05_1_png, origin='lower')
    im2 = gax_k05_2.imshow(k05_2_png, origin='lower')
    im3 = gax_k05_3.imshow(k05_3_png, origin='lower')
    im4 = gax_k06_1.imshow(k06_1_png, origin='lower')
    im5 = gax_k06_2.imshow(k06_2_png, origin='lower')
    im6 = gax_k06_3.imshow(k06_3_png, origin='lower')
    im7 = gax_k07_1.imshow(k07_1_png, origin='lower')
    im8 = gax_k07_2.imshow(k07_2_png, origin='lower')
    im9 = gax_k07_3.imshow(k07_3_png, origin='lower')

    for gax in the_fig_low_k.axes:
        gax.set_xticks([])
        gax.set_yticks([])

    the_fig_low_k.savefig('output/figs/ExampleIntersections_LowK.pdf')

    k08_1_png = Image.open('input/type_intersections/k08_1.png')
    k08_2_png = Image.open('input/type_intersections/k08_2.png')
    k08_3_png = Image.open('input/type_intersections/k08_3.png')
    k09_1_png = Image.open('input/type_intersections/k09_1.png')
    k09_2_png = Image.open('input/type_intersections/k09_2.png')
    k09_3_png = Image.open('input/type_intersections/k09_3.png')
    k10_1_png = Image.open('input/type_intersections/k10_1.png')
    k10_2_png = Image.open('input/type_intersections/k10_2.png')
    k10_3_png = Image.open('input/type_intersections/k10_3.png')

    the_fig_hi_k = plt.figure(figsize=(6,6))
    gax_k08_1 = the_fig_hi_k.add_subplot(331)
    gax_k08_1.set_ylabel('k=0.8')
    gax_k08_2 = the_fig_hi_k.add_subplot(332)
    gax_k08_2.set_title('High Confidence Type Intersections')
    gax_k08_3 = the_fig_hi_k.add_subplot(333)
    gax_k09_1 = the_fig_hi_k.add_subplot(334)
    gax_k09_1.set_ylabel('k=0.9')
    gax_k09_2 = the_fig_hi_k.add_subplot(335)
    gax_k09_3 = the_fig_hi_k.add_subplot(336)
    gax_k10_1 = the_fig_hi_k.add_subplot(337)
    gax_k10_1.set_ylabel('k=1.0')
    gax_k10_2 = the_fig_hi_k.add_subplot(338)
    gax_k10_3 = the_fig_hi_k.add_subplot(339)

    im1 = gax_k08_1.imshow(k08_1_png, origin='lower')
    im2 = gax_k08_2.imshow(k08_2_png, origin='lower')
    im3 = gax_k08_3.imshow(k08_3_png, origin='lower')
    im4 = gax_k09_1.imshow(k09_1_png, origin='lower')
    im5 = gax_k09_2.imshow(k09_2_png, origin='lower')
    im6 = gax_k09_3.imshow(k09_3_png, origin='lower')
    im7 = gax_k10_1.imshow(k10_1_png, origin='lower')
    im8 = gax_k10_2.imshow(k10_2_png, origin='lower')
    im9 = gax_k10_3.imshow(k10_3_png, origin='lower')

    for gax in the_fig_hi_k.axes:
        gax.set_xticks([])
        gax.set_yticks([])

    the_fig_hi_k.savefig('output/figs/ExampleIntersections_HiK.pdf')


# Possible Figures:
# -----------------
# A.) Significance of correlation:
#     - Histogram of random orderings and their agreement with the GSN
#     - Vertical line on the histogram showing measured agreement
#     - Fit gaussian curve to the histogram:
#       - find mean and standard deviation
#       - based on fit, calculate probability that measured value is chance

# B.) Degradation of correlation with number of lineament transpositions:
#     - Generate GSN from NSR linstack
#     - use test_agreement to perform lineament transpositions
#     - run many tests, showing how agreement degrades with number of transpositions
#     - How variable is the degradation?
#     - How many transpositions does it take on average to get to our measured A=0.70?
#     - Plot the curves with x-axis = # of transpositions, and y-axis = A
#     - Needs to be normalized according to the average path confidence one
#       would expect given the distribution of intersection confidences that we
#       have in the real GSN and unfortunately, this can't currently be done.

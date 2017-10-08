from __future__ import division
from cctbx.array_family import flex
from libtbx.development.timers import Profiler
import math
"""
Initial implementation of clustering by Rodriguez 2014 algorithm
"""
def run_detail(show_plot, save_plot):
    P = Profiler("0. Read data")
    import pickle
    DP = pickle.load(open("mbhresults.pickle","rb"))
    coord_x = DP.get("coord_x")#[:500]
    coord_y = DP.get("coord_y")#[:500]

    print("There are %d points"%(len(coord_x)))
    if show_plot or save_plot:
      import matplotlib
      if not show_plot:
        # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
        matplotlib.use('Agg') # use a non-interactive backend
      from matplotlib import pyplot as plt
      #plt.plot(coord_x.select(selection),coord_y.select(selection),"r.", markersize=2.)
      #plt.plot(coord_x.select(~selection),coord_y.select(~selection),"k.", markersize=3.)
      plt.plot(coord_x,coord_y,"k.", markersize=3.)
      #plt.plot(coord_x.select(selection_1F400),coord_y.select(selection_1F400),"r.", markersize=3.)
      #plt.plot(coord_x.select(selection_2F400),coord_y.select(selection_2F400),"k.", markersize=3.)
      plt.axes().set_aspect("equal")
      if save_plot:
        plt.savefig(plot_name,
                    size_inches=(10,10),
                    dpi=300,
                    bbox_inches='tight')
      if show_plot:
        plt.show()

    print "Now constructing a Dij matrix.  Do it faster with C++ outer product"
    #from IPython import embed; embed()
    P = Profiler("1. compute Dij matrix")
    NN = len(coord_x)
    embedded_vec2 = flex.vec2_double(coord_x,coord_y)
    Dij_C = embedded_vec2.distance_matrix(embedded_vec2)
    print Dij_C.focus()

    Dij = flex.double(flex.grid(NN,NN))
    print Dij.focus()
    for i in xrange(NN):
      for j in xrange(NN):
        Dij[(i,j)] = math.sqrt( (coord_x[i]-coord_x[j])**2 + (coord_y[i]-coord_y[j])**2 )
        assert Dij[(i,j)] == Dij_C[(i,j)]

    P = Profiler("2. calculate rho density")
    print "finished Dij, now calculating rho_i, the density"
    max_Dij = flex.max(Dij)
    d_c = 0.04 # a guess, for now
    rho_i = flex.size_t(NN)
    for i in xrange(NN):
      for j in xrange(NN):
        if Dij[(i,j)] < d_c:  rho_i[i]+=1
    ave_rho = flex.mean(rho_i.as_double())
    print "The average rho_i is %5.2f, or %4.1f%%"%(ave_rho, 100*ave_rho/NN)
    i_max = flex.max_index(rho_i)
    P = Profiler("3.transition")
    print "the index with the highest density is %d"%(i_max)
    delta_i_max = flex.max(flex.double([Dij[i_max,j] for j in xrange(NN)]))
    print "delta_i_max",delta_i_max
    rho_order = flex.sort_permutation(rho_i,reverse=True)
    rho_order_list = list(rho_order)
    delta_i = flex.double(NN,delta_i_max)
    # delta_i is measured by computing the minimum distance between the point i
    # and any other point with higher OR EQUAL density (emphasis mine)
    P = Profiler("4. delta")
    for p in xrange(1,NN): # first iteration through rho order
      i = rho_order[p]
      for q in xrange(p):
        j = rho_order[q]
        #print p,q,i,j
        if rho_i[j] >= rho_i[i] and Dij[(i,j)] < delta_i[i]:
          delta_i[i] = Dij[(i,j)]
    P = Profiler("5. find cluster maxima")
    #---- Now hunting for clusters
    cluster_id = flex.int(NN,-1) # default -1 means no cluster
    delta_order = flex.sort_permutation(delta_i,reverse=True)
    N_CLUST = 10 # maximum of 10 points to be considered as possible clusters
    MAX_PERCENTILE_DELTA = 0.10 # cluster centers have to be in the top 10% percentile delta
    MAX_PERCENTILE_RHO = 0.75 # cluster centers have to be in the top 75% percentile rho
    n_cluster = 0
    max_n_delta = min(N_CLUST, int(MAX_PERCENTILE_DELTA*NN))
    for ic in xrange(max_n_delta):
      # test the density, rho
      item_idx = delta_order[ic]
      if delta_i[item_idx] < 0.25 * delta_i[delta_order[0]]: # too low (another heuristic!)
        continue
      item_rho_order = rho_order_list.index(item_idx)
      if item_rho_order/NN < MAX_PERCENTILE_RHO :
        cluster_id[item_idx] = n_cluster
        print ic,item_idx,item_rho_order,cluster_id[item_idx]
        n_cluster += 1
    print "Found %d clusters"%n_cluster
    for x in xrange(NN):
      if cluster_id[x]>=0:
        print "XC",x,cluster_id[x],rho_i[x],delta_i[x]

    from matplotlib import pyplot as plt
    #plt.plot(rho_i,delta_i,"r.", markersize=3.)
    #for x in xrange(NN):
    #  if cluster_id[x]>=0:
    #    plt.plot([rho_i[x]],[delta_i[x]],"ro")
    #plt.show()

#start here.
# we are not assigning remaining points correctly
    P = Profiler("6. assign all points")
    # one pass to identify cluster id
    # assign each point to its nearest neighbor (Dij) of higher density
    for p in xrange(NN):
      item_idx = rho_order[p]
      if cluster_id[item_idx] == -1: # still unassigned
        trial_Dij = max_Dij
        i_neighbor = None
        for q in xrange(p):
          if Dij[(item_idx,rho_order[q])] < trial_Dij:
            i_neighbor = rho_order[q]
            trial_Dij = Dij[(item_idx,rho_order[q])]
        cluster_id[item_idx] = cluster_id[i_neighbor]
    #for idx in xrange(NN):
    #
    #  plt.plot([coord_x[idx]],[coord_y[idx]],"%so"%(
    #    {-1:'k',0:'b',1:'r',2:'g'}[cluster_id[idx]]
    #))
    #plt.show()
# assign the halos
    P = Profiler("7. assign halos")
    halo = flex.bool(NN,False)
    border = flex.bool(NN,False)
    for i in xrange(0,NN):
      for j in xrange(i+1,NN):
        if Dij[(i,j)] < d_c: # find points being within d_c of those in another cluster
          #print Dij[(i,j)] , d_c, i,j, cluster_id[i], cluster_id[j]
          if cluster_id[i] != cluster_id[j]:
            border[i]=True; border[j]=True
    for ic in range(n_cluster): #loop thru all border regions; find highest density
      print "cluster",ic, "in border",border.count(True)
      this_border = (cluster_id == ic) & (border==True)
      print len(this_border), this_border.count(True)
      if this_border.count(True)>0:
        highest_density = flex.max(rho_i.select(this_border))
        halo_selection = (rho_i < highest_density) & (this_border==True)
        if halo_selection.count(True)>0:
          cluster_id.set_selected(halo_selection,-1)
        core_selection = (cluster_id == ic) & ~halo_selection
        highest_density = flex.max(rho_i.select(core_selection))
        #from IPython import embed; embed()
        too_sparse = core_selection & (rho_i.as_double() < highest_density/10.) # another heuristic
        if too_sparse.count(True)>0:
          cluster_id.set_selected(too_sparse,-1)

    for idx in xrange(NN):
      if cluster_id[idx]==-1:
        plt.plot([coord_x[idx]],[coord_y[idx]],"%s."%(
        {-1:'k',0:'b',1:'r',2:'g'}[cluster_id[idx]]))
      else:
        plt.plot([coord_x[idx]],[coord_y[idx]],"%so"%(
        {-1:'k',0:'b',1:'r',2:'g'}[cluster_id[idx]]))

    plt.show()


if __name__=="__main__":

  run_detail(show_plot=False, save_plot=False)
""" Benchmark, 6672 MBH lattices
                   Python/Flex arrays C++ Dij
         0. Read data: CPU,    0.000s;  0.00s;
1. compute Dij matrix: CPU,  171.810s;  1.92s;
 2.calculate rho dens: CPU,  129.620s;143.39s;
         3.transition: CPU,    0.020s;  0.03s;
             4. delta: CPU,  121.590s;121.04s;
5.find cluster maxima: CPU,    0.570s;  0.59s;
 6. assign all points: CPU,   83.540s; 84.43s;
      7. assign halos: CPU,   63.750s; 87.83s;
TOTAL                : CPU,  570.900s;439.23s;

"""
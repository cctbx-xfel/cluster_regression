from __future__ import division, print_function
from six.moves import range
from cctbx.array_family import flex
from libtbx.development.timers import Profiler
from libtbx import group_args
"""
Initial implementation of clustering by Rodriguez 2014 algorithm
"""
class clustering_manager(group_args):
  def __init__(self, **kwargs):
    group_args.__init__(self, **kwargs)
    # require Dij, d_c
    P = Profiler("2. calculate rho density")
    print("finished Dij, now calculating rho_i, the density")
    from xfel.clustering import Rodriguez_Laio_clustering_2014
    R = Rodriguez_Laio_clustering_2014(distance_matrix = self.Dij, d_c = self.d_c)
    self.rho = rho = R.get_rho()
    ave_rho = flex.mean(rho.as_double())
    NN = self.Dij.focus()[0]
    print("The average rho_i is %5.2f, or %4.1f%%"%(ave_rho, 100*ave_rho/NN))
    i_max = flex.max_index(rho)

    P = Profiler("3.transition")
    print("the index with the highest density is %d"%(i_max))
    delta_i_max = flex.max(flex.double([self.Dij[i_max,j] for j in range(NN)]))
    print("delta_i_max",delta_i_max)
    rho_order = flex.sort_permutation(rho,reverse=True)
    rho_order_list = list(rho_order)

    P = Profiler("4. delta")
    self.delta = delta = R.get_delta(rho_order=rho_order, delta_i_max=delta_i_max)

    P = Profiler("5. find cluster maxima")
    #---- Now hunting for clusters
    cluster_id = flex.int(NN,-1) # default -1 means no cluster
    delta_order = flex.sort_permutation(delta,reverse=True)
    N_CLUST = 10 # maximum of 10 points to be considered as possible clusters
    MAX_PERCENTILE_DELTA = 0.10 # cluster centers have to be in the top 10% percentile delta
    MAX_PERCENTILE_RHO = 0.75 # cluster centers have to be in the top 75% percentile rho
    n_cluster = 0
    max_n_delta = min(N_CLUST, int(MAX_PERCENTILE_DELTA*NN))
    for ic in range(max_n_delta):
      # test the density, rho
      item_idx = delta_order[ic]
      if delta[item_idx] < 0.25 * delta[delta_order[0]]: # too low (another heuristic!)
        continue
      item_rho_order = rho_order_list.index(item_idx)
      if item_rho_order/NN < MAX_PERCENTILE_RHO :
        cluster_id[item_idx] = n_cluster
        print(ic,item_idx,item_rho_order,cluster_id[item_idx])
        n_cluster += 1
    print("Found %d clusters"%n_cluster)
    for x in range(NN):
      if cluster_id[x]>=0:
        print("XC",x,cluster_id[x],rho[x],delta[x])
    self.cluster_id_maxima = cluster_id.deep_copy()

    P = Profiler("6. assign all points")
    R.cluster_assignment(rho_order,cluster_id)

    self.cluster_id_full = cluster_id.deep_copy()

    # assign the halos
    P = Profiler("7. assign halos")
    halo = flex.bool(NN,False)
    border = R.get_border( cluster_id = cluster_id )

    for ic in range(n_cluster): #loop thru all border regions; find highest density
      print("cluster",ic, "in border",border.count(True))
      this_border = (cluster_id == ic) & (border==True)
      print(len(this_border), this_border.count(True))
      if this_border.count(True)>0:
        highest_density = flex.max(rho.select(this_border))
        halo_selection = (rho < highest_density) & (this_border==True)
        if halo_selection.count(True)>0:
          cluster_id.set_selected(halo_selection,-1)
        core_selection = (cluster_id == ic) & ~halo_selection
        highest_density = flex.max(rho.select(core_selection))
        too_sparse = core_selection & (rho.as_double() < highest_density/10.) # another heuristic
        if too_sparse.count(True)>0:
          cluster_id.set_selected(too_sparse,-1)
    self.cluster_id_final = cluster_id.deep_copy()
    print("%d in the excluded halo"%((cluster_id==-1).count(True)))

def run_detail(show_plot, save_plot):
    P = Profiler("0. Read data")
    import pickle
    DP = pickle.load(open("mbhresults.pickle","rb"))
    coord_x = DP.get("coord_x")
    coord_y = DP.get("coord_y")
    del P

    print(("There are %d points"%(len(coord_x))))
    if show_plot or save_plot:
      import matplotlib
      if not show_plot:
        # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
        matplotlib.use('Agg') # use a non-interactive backend
      from matplotlib import pyplot as plt
      plt.plot(coord_x,coord_y,"k.", markersize=3.)
      plt.axes().set_aspect("equal")
      if save_plot:
        plt.savefig(plot_name,
                    size_inches=(10,10),
                    dpi=300,
                    bbox_inches='tight')
      if show_plot:
        plt.show()

    print("Now constructing a Dij matrix.")
    P = Profiler("1. compute Dij matrix")
    NN = len(coord_x)
    embedded_vec2 = flex.vec2_double(coord_x,coord_y)
    Dij = embedded_vec2.distance_matrix(embedded_vec2)
    del P

    d_c = 0.04 # the distance cutoff, such that average item neighbors 1-2% of all items
    CM = clustering_manager(Dij=Dij, d_c=d_c)

    appcolors = ['b', 'r', '#ff7f0e', '#2ca02c',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
    if show_plot:
      #Decision graph
      from matplotlib import pyplot as plt
      plt.plot(CM.rho,CM.delta,"r.", markersize=3.)
      for x in range(NN):
        if CM.cluster_id_maxima[x]>=0:
          plt.plot([CM.rho[x]],[CM.delta[x]],"ro")
      plt.show()

      #No-halo plot
      from matplotlib import pyplot as plt
      colors = [appcolors[i%10] for i in CM.cluster_id_full]

      plt.scatter(coord_x,coord_y,marker='o',color=colors,
          linewidths=0.4, edgecolor='k')
      plt.axes().set_aspect("equal")
      plt.show()

      #Final plot
      halo = (CM.cluster_id_final==-1)
      core = ~halo
      plt.plot(coord_x.select(halo),coord_y.select(halo),"k.")
      colors = [appcolors[i%10] for i in CM.cluster_id_final.select(core)]
      plt.scatter(coord_x.select(core),coord_y.select(core),marker="o",
          color=colors,linewidths=0.4, edgecolor='k')
      plt.axes().set_aspect("equal")
      plt.show()

if __name__=="__main__":

  run_detail(show_plot=True, save_plot=False)
""" Benchmark, 6672 MBH lattices
                   Python/Flex arrays;C++ Dij;C++ rho;C++delt;C clust;C bordr;
         0. Read data: CPU,    0.000s;  0.00s;  0.01s;  0.00s;  0.00s;  0.00s;
1. compute Dij matrix: CPU,  171.810s;  1.92s;  1.47s;  1.11s;  1.66s;  2.02s;
 2.calculate rho dens: CPU,  129.620s;143.39s;  0.27s;  0.27s;  0.19s;  0.28s;
         3.transition: CPU,    0.020s;  0.03s;  0.02s;  0.04s;  0.02s;  0.02s;
             4. delta: CPU,  121.590s;121.04s;114.22s;  0.27s;  0.23s;  0.40s;
5.find cluster maxima: CPU,    0.570s;  0.59s;  0.56s;  0.48s;  0.47s;  0.53s;
 6. assign all points: CPU,   83.540s; 84.43s; 79.74s; 90.77s;  0.36s;  0.72s;
      7. assign halos: CPU,   63.750s; 87.83s; 77.63s; 82.06s; 83.06s;  0.06s;
TOTAL                : CPU,  570.900s;439.23s;273.92s;175.00s; 85.99s;  4.03s;
"""

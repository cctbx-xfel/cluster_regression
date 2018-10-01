from __future__ import division
from six.moves import range
import os
import libtbx.load_env
# must have numpy--Will crash with OpenMP otherwise, when importing singleframe
import numpy as np # import dependency
from xfel.clustering.singleframe import SingleFrame
from cctbx.uctbx.determine_unit_cell import NCDist,NCDist2017

"""Program performs a comparison between two implementations of NCDist:
1) original version from Andrews/Berstein checked in to cctbx_project
   repo by Oliver Zeldin in ~2015.
2) live NCDist2017 version, from current yayahjb/ncdist repo
"""

def generate_unit_cells_from_text(path):
  data = open (path,"r").readlines()
  for line in data:
    tokens = line.strip().split()
    cell_param = tuple([float(t) for t in tokens[0:6]])
    yield cell_param

def run_one(path):
  cells = [ g for g in generate_unit_cells_from_text(path) ]
  g6 = [ SingleFrame.make_g6(u) for u in cells ]

  # for the purpose of this test, cycle through pairs of g6 vectors
  for ix in range(len(g6)-1):
    a = g6[ix]
    b = g6[ix + 1]
    old = NCDist(a,b)
    # workaround allows use of non-thread-safe NCDist, even if openMP is enabled elsewhere in the Python program
    import os,omptbx
    workaround_nt = int(os.environ.get("OMP_NUM_THREADS",1))
    omptbx.omp_set_num_threads(1)
    new = NCDist2017(a,b)
    com = NCDist2017(b,a)
    omptbx.omp_set_num_threads(workaround_nt)
    assert old==new, "Zeldin, AB2017"
    assert new==com, "Pair %d NCDist(a,b) %f != NCDist(b,a) %f"%(ix,new,com)

def run_all():
  textfiles = ["lysozyme1341.txt"]
  for f in textfiles:
    calib_dir = libtbx.env.find_in_repositories("cluster_regression")
    path = os.path.join(calib_dir,"examples",f)
    run_one(path)

if __name__=="__main__":
  for i in range(10):
    print i
    run_all()
  print "OK"

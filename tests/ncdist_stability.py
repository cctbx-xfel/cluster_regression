from __future__ import division
import os
import libtbx.load_env
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
  for ix in xrange(len(g6)-1):
    a = g6[ix]
    b = g6[ix + 1]
    old = NCDist(a,b)
    new = NCDist2017(a,b)
    com = NCDist2017(b,a)
    assert old==new
    assert new==com

def run_all():
  textfiles = ["lysozyme1341.txt"]
  for f in textfiles:
    calib_dir = libtbx.env.find_in_repositories("cluster_regression")
    path = os.path.join(calib_dir,"examples",f)
    run_one(path)

if __name__=="__main__":
  run_all()
  print "OK"

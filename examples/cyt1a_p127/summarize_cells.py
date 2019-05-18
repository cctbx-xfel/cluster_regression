from __future__ import division, print_function
import sys

"""Simply convert integration pickles to a text file with unit cell and
space group symbol"""
if __name__=="__main__":
  args = sys.argv[1:]
  filename = args[0]
  from six.moves import cPickle as pickle
  data = pickle.load(open(filename,"rb"))

  obs = data["observations"][0]

  unit_cell = obs.unit_cell()
  params = unit_cell.parameters()
  print (" ".join([p for p in params]))

  space_group = obs.space_group_info().type().lookup_symbol()
  print("".join(space_group.split()), end=' ')

  import os.path
  print(os.path.basename(filename))


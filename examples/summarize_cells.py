from __future__ import division, print_function
import sys

"""Simply convert integration pickles to a text file with unit cell and
space group symbol"""
from xfel.command_line.print_pickle import generate_data_from_streams
if __name__=="__main__":
  args = sys.argv[1:]

  for data in generate_data_from_streams(args, verbose=False):

    obs = data["observations"][0]

    unit_cell = obs.unit_cell()
    params = unit_cell.parameters()
    print(" ".join([p for p in params]))


    space_group = obs.space_group_info().type().lookup_symbol()
    print("".join(space_group.split()))

from __future__ import division, print_function
import sys
data = open (sys.argv[1],"r").readlines()
a_cell = []
c_cell = []
for line in data:
  tokens = line.strip().split()
  a = float(tokens[0])
  c = float(tokens[2])
  a_cell.append(a)
  c_cell.append(c)

from matplotlib import pyplot as plt
plt.plot(a_cell,c_cell,"r,")
plt.axes().set_aspect("equal")
plt.show()

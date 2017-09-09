from __future__ import division
from libtbx import test_utils
import libtbx.load_env

tst_list = (
  "$D/tests/ncdist_stability.py",
  )

def run_standalones():
  build_dir = libtbx.env.under_build("cluster_regression")
  dist_dir = libtbx.env.dist_path("cluster_regression")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run_standalones()

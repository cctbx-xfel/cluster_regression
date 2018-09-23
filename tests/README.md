Instructions to test Andrews/Bernstein NCDist in the context of CCTBX. rev. 9/23/2018. Purpose:
1) Develop with confidence that nothing has broken (assert that past results are repeated)
2) Compare results and performance of G6, S6, and D7 metrics.
3) Measure OpenMP performance.

Pre-requisites.
1) Linux box with multiple cores (for OpenMP).  Use bash shell.
2) Make sure your conda environment file (~/.condarc) has the following:
channels:
  - cctbx
  - conda-forge 
  - defaults
  - bioconda

report_errors: false
3) Create and change to a working directory ${WORK}

Download three files to configure the environment
  (cctbx bootstrap installer, conda-Python dependency definitions, and miniconda installer [tested with conda 4.5.2]):

wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py --no-check-certificate
wget https://raw.githubusercontent.com/cctbx-xfel/cluster_regression/master/tests/dials_env.txt --no-check-certificate
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh --no-check-certificate
chmod u+x Miniconda2-latest-Linux-x86_64.sh
./Miniconda2-latest-Linux-x86_64.sh
# install in directory ${WORK}/miniconda2
# do not prepend Miniconda2 to .bashrc

source miniconda2/etc/profile.d/conda.sh
conda create -y --name myEnv --file dials_env.txt
conda activate myEnv
which python # verify that you are using the miniconda myEnv python

python -m pip install procrunner
python bootstrap.py hot --builder=dials
python bootstrap.py update --builder=dials
cd ${WORK}/modules
git clone git@github.com:cctbx-xfel/cluster_regression.git
mkdir ${WORK}/build
cd ${WORK}/build
python ../modules/cctbx_project/libtbx/configure.py --enable_openmp_if_possible=True cluster_regression

source build/setpaths.sh
cd build; make; cd -
mkdir test
cd test
rm -rf *
export OMP_NUM_THREADS=8
libtbx.python ${WORK}/modules/cluster_regression/tests/public-test-all.py

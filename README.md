ssa
===

Sparse Surface Adjustment:
SSA is an open-source C++ tool for post optimization of 2D or 3D SLAM solutions. 
SSA iteratively refines robot poses and surface points in one global graph optimization 
system and produces highly accurate 3D point clouds. 
3D sensor observations are not treated as rigid body and might be refined during the 
optimization procedure. This leads to substantially less accumulated noise in the 
resulting model. SSA builds upon g2o as optimization back-end.

More detailed information can be found in the Sparse Surface Adjustment publication:
http://ais.informatik.uni-freiburg.de/publications/papers/ruhnke12icra.pdf

Authors:
Michael Ruhnke; Rainer Kuemmerle; Giorgio Grisetti; Wolfram Burgard;

License:
SSA is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
SSA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

===

How to get the code?
git clone https://github.com/MichaelRuhnke/ssa.git

Dependencies:
- libqglviewer (for opengl visualization)
- Suitesparse  (math library for sparse matrices) 
  SSA offers an option to download and compile openBLAS and suitesparse.
  usually this gives a speedup up between 2 and 10 for the optimization.
- Pointcloud Library (for import / export of data)
- FLANN (also a dependency of pcl)
- Eigen3 (also a dependency of pcl)
- g2o (will be checked out from the developer git repository during the build process)

  On Ubuntu machines you might install the dependencies with:
  sudo apt-get install libqglviewer-qt4-dev 
  sudo apt-get install libsuitesparse-dev
  sudo add-apt-repository ppa:v-launchpad-jochen-sprickerhof-de/pcl
  sudo apt-get update
  sudo apt-get install libpcl-all

Compiling ssa:
mkdir build
cd build
cmake-gui ../
#if qmake-NOTFOUND appears change it into qmake-qt4
make

Run an example:
There is an example data set of a dark mug in example/alufr_black_mug/

First step to convert a dataset into ssa format is: 
pcd_to_ssa

$ pcd_to_ssa: converts a set of PointXYZRGBA pcd files into a ssa 3d graph file
$ usage pcd_to_ssa [options] <ssa3d_file> options:
$ -p [string]	prefix of the pcd file list.
$ -r [double]	resolution of the resulting model (default resolution 0.01m)
$ -n [int] 	count of point clouds
$ example:
$ pcd_to_ssa -p alufr_black_mug_raw -n 60 -r 0.001 alufr_black_mug_raw.ssa3d

this will load alufr_black_mug_raw00000.pcd to alufr_black_mug_raw00059.pcd
and subsample to a resolution of 1mm.

The pcd files should have proper sensor poses and orientations (from your registration
algorithm). The points itself should be in sensor frame and not transformed!!!

> ssa_viewer_3d -ini ssa.ini alufr_black_mug_raw.ssa3d

runs the optimizer gui. An example ssa.ini file is located in the ssa root dir.

If you apply parameters that result in an unconnected graph, the optimizer will crash.
Right now I do not check for the connectivity, but I will provide that feature pretty soon.

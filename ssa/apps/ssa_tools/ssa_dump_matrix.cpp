// Sparse Surface Optimization
// Copyright (C) 2011 M. Ruhnke, R. Kuemmerle, G. Grisetti, W. Burgard
// 
// SSA2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// SSA2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <signal.h>

#include "ssa/core/ssa_graph_3d.h"

#include "ssa/core/allocate_solver.h"
#include "ssa/core/sparse_surface_adjustment.h"

using namespace std;
using namespace ssa;

const char *message[]={
  "ssa_dump_matrix: dumps non zero entries of approximated hessian into ssa-matrix-dump.txt",
  "usage ssa_dump_matrix [options] <ssa_file>",
  "options:",
  "-nn    apply pairwise nearest neighbor data association.",
  "-ns    apply normal shooting data association.",
  "-tree  use tree approximation as data association strategie.",
  0
};

int main(int argc, char **argv)
{

  const char**v=message;
  while (*v){
    cout << *v << endl;
    v++;
  }

  // HACK reset numeric locale, we need C
  setlocale (LC_NUMERIC,"C");
  int c=1;
  const char* logfile=0;
  bool useFullNN = false;
  bool useNormalShooting = false;
  bool useTreeNN = false;

  bool useConfigFile = false;
  const char* configFile=0;

  while (c<argc){
    if (!strcmp(argv[c],"-debug")){
//       debug=true;
      c++;
    } else 
    if (!strcmp(argv[c],"-nn")){
      useFullNN = true;
      c++;
    } else
    if (!strcmp(argv[c],"-ns")){
      useNormalShooting = true;
      c++;
    } else
    if (!strcmp(argv[c],"-tree")){
      c++;
      useTreeNN = true;
    } else
    if (!strcmp(argv[c],"-ini")){
      useConfigFile=true;
      c++;
      configFile = argv[c];
      c++;
    } else 
    if (! logfile){
      logfile=argv[c];
      c++;
      break;
    }
  }

  if(!logfile){
    cerr << "Error: no input ssa file given." << endl;
    exit(0);
  }

  /** initialize ssa */
  SparseSurfaceAdjustmentT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov> ssa;
  ssa.setVerbose(true);
  if(logfile){
    ssa.graph()->load(logfile);
  }

  /** Load configuration from ini file */
  if(useConfigFile){
    ifstream configStream(configFile);
    ssa.params().readParams(configStream);
    configStream.close();
    ssa.params().printParams();
  }

  /** recalc data association */
  if(useFullNN || useNormalShooting || useTreeNN)
    ssa.graph()->dropDataAssociation();

  if(useFullNN)
    NearestNeighbor<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::apply(*ssa.graph(), ssa.params().nearestNeighbor);

  if(useNormalShooting)
    NormalShootingFlann<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::shootNormals(*ssa.graph(), ssa.params().normalShooting);

  if(useTreeNN){
    NearestNeighbor<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::applyStarAssignment(*ssa.graph(), ssa.params().nearestNeighbor);
    //NearestNeighbor<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::pruneToMinSpanTree(*ssa.graph());
  }

  ssa.dumpConnectivityMatrix();
}

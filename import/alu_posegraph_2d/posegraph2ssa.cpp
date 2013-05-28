// Sparse Surface Optimization 2D
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
#include <stdlib.h>
#include <string.h>
#include <istream>
#include <fstream>
#include <sstream>
#include <list>
#include <cstring>
#include <limits>
#include <cmath>

#include "import_posegraph2d.h"
#include "ssa/core/ssa_graph_2d.h"
#include "ssa/sensor_models/laser_sensor_model_2d.h"

using namespace std;
using namespace g2o;
using namespace ssa;

const char *message[]={
  "ssa_from_posegraph: converts gm2dl file into a ssa graph file",
  "usage ssa_from_posegraph [options] <gm2dl_file> <ssa_file>",
  "options:",
  "-window [meter]     radius of the neighorhood search for the normal calculation.",
  0
};


int main(int argc, char **argv)
{

 if (argc<2){
    const char**v=message;
    while (*v){
      cout << *v << endl;
      v++;
    }
    return 0;
  }

  bool debug=false;
  //double maxrange=82;

  const char* logfile=0;
  const char* outputFile=0;
  double windowSize = 0.5; //radius

  int c=1;
  while (c<argc){
    if (!strcmp(argv[c],"-debug")){
      debug=true;
      c++;
    } else
    if (!strcmp(argv[c],"-window")){
      c++;
      windowSize=atof(argv[c]);
      c++;
    } else
//     if (!strcmp(argv[c],"-cutoff")){
//       c++;
//       cutoff=atof(argv[c]);
//       c++;
//     } else
//     if (!strcmp(argv[c],"-maxrange")){
//       c++;
//       maxrange=atof(argv[c]);
//       c++;
//     } else
    if (! logfile){
      logfile=argv[c];
      c++;
    } else
    if (! outputFile){
      outputFile=argv[c];
      c++;
      break;
    }
  }
  cerr << logfile << " >> " << outputFile << std::endl;

  SparseSurfaceAdjustmentGraph2D ssaGraph;
  SSAPoseGraph2D::importPoseGraph2D(logfile, ssaGraph);
  cerr << "constructed graph with " << ssaGraph._optimizer.vertices().size() << " vertices and " << ssaGraph._optimizer.edges().size() << " edges." << endl;

  SparseSurfaceAdjustmentParams ssaParams;
  ssaParams.normalExtractionMaxNeighborDistance = windowSize;
  ssaParams.normalExtractionMinNeighbors = 3;
  ssaParams.normalExtractionMaxNeighbors = 10;
  ssaGraph.fillNeighborCache(ssaParams);
  ssaGraph.calcMeanCov(ssaParams);

  LaserSensorParams params;
  params.maxRange = 80.0;
  params.angularResolution = 0.5;
  params.sensorPrecision = 0.01;
  params.scale = 1e-3;
  LaserSensorModel2D::applySensorModel(ssaGraph, params);
  ssaGraph.save(outputFile);
}


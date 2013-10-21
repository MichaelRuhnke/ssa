// Sparse Surface Optimization
// Copyright (C) 2011 M. Ruhnke, R. Kuemmerle, G. Grisetti, W. Burgard
//
// SSA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SSA is distributed in the hope that it will be useful,
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
#include <ostream>
#include <fstream>
#include <sstream>
#include <list>
#include <cstring>
#include <limits>
#include <cmath>

#include "ssa/core/ssa_graph_2d.h"
#include "ssa/core/sparse_surface_adjustment.h"

using namespace std;
using namespace g2o;
using namespace ssa;

const char *message[]={
  "ssa_to_posegraph: converts gm2dl and a ssa file into a gm2dl file",
  "usage ssa_to_posegraph [options] <gm2dl_file> <ssa_file> <gm2dl_file>",
  "options:",
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

  //double maxrange=82;

  const char* logfile=0;
  const char* ssafile=0;
  const char* outputFile=0;

  int c=1;
  while (c<argc){
    if (! logfile){
      logfile=argv[c];
      c++;
    } else
    if (! ssafile){
      ssafile=argv[c];
      c++;
    } else
    if (! outputFile){
      outputFile=argv[c];
      c++;
      break;
    }
  }

  /** Load ssa graph */
  ssa::SparseSurfaceAdjustmentT<g2o::EdgeSE2, ssa::EdgeSE2PointXYCov, ssa::EdgePointXYCovPointXYCov> ssa;
  ssa.setVerbose(true);

  if (! ssafile){
    cerr << "Missing ssa graph file." << endl;
    exit(-1);
  } else {
    ssa.graph()->load(ssafile);
  }

  cerr << logfile << " and  " << ssafile << " -> " << outputFile << std::endl;
  ofstream outFile(outputFile);
  ifstream is(logfile);
  if (! is ){
    cerr << "error in loading the graph" << endl;
    return 0;
  }

  VertexSE2* previousVertex = 0;
  while(is){
      char buf[40960];
      is.getline(buf,40960);
      istringstream ls(buf);
      string tag;
      ls >> tag;
      if (tag=="VERTEX" || tag=="VERTEX2" || tag=="VERTEX_SE2"){
        int id;
        double x,y,theta;
        ls >> id >> x >> y >> theta;

        VertexSE2* v = ssa.graph()->_verticies_poses[id];
        if(v->id() == id)
          outFile << tag.c_str() << " " << id << " " << v->estimate().translation()(0) << " " << v->estimate().translation()(1) << " " << v->estimate().rotation().angle() << std::endl;
        previousVertex = v;
      } else if(tag=="ROBOTLASER1" && previousVertex)
      {
        int type;
        double angle, fov, res, maxrange, acc;
        int remission_mode;
        int beams, remisioncount;
        std::vector<double> ranges;
        std::vector<double> remissions;
        double x,y,theta;
        double sensorX,sensorY,sensorTheta;
        double laserTv, laserRv, forwardSafetyDist, sideSaftyDist, turnAxis, timestamp, loggerTimestamp;
        string hostname;

        //parsing log line
        ls >> type >> angle >> fov >> res >> maxrange >> acc >> remission_mode;
        ls >> beams;
        ranges.resize(beams);
        for (int i=0; i<beams; i++)
          ls >> ranges[i];

        ls >> remisioncount;
        remissions.resize(remisioncount);
        for (int i = 0; i < remisioncount; i++)
          ls >> remissions[i];

        // special robot laser stuff
        ls >> x >> y >> theta;
        ls >> sensorX >> sensorY >> sensorTheta;
        ls >> laserTv >>  laserRv >>  forwardSafetyDist >> sideSaftyDist >> turnAxis;

        // timestamp + host
        ls >> timestamp;
        ls >> hostname;
        ls >> loggerTimestamp;

        SE2 robotPose = SE2(x,y,theta);
        SE2 laserPose = SE2(sensorX,sensorY,sensorTheta);
        SE2 laserOffSet = robotPose.inverse()*laserPose;

        x = previousVertex->estimate().translation()(0);
        y = previousVertex->estimate().translation()(1);
        theta = previousVertex->estimate().rotation().angle();
        laserPose = previousVertex->estimate() * laserOffSet;
        outFile << fixed << tag.c_str() << " " << type << " " << angle << " " << fov << " " << res << " " << maxrange << " " << acc << " " << remission_mode;
        outFile << " " << beams;
        for (int i=0; i<beams; i++)
          outFile << " " << ranges[i];

        outFile << " " << remisioncount;
        for (int i = 0; i < remisioncount; i++)
          outFile << " " << remissions[i];

        // special robot laser stuff
        outFile << " " << x << " " << y << " " << theta;
        outFile << " " << laserPose.toVector()[0] << " " << laserPose.toVector()[1] << " " << laserPose.toVector()[2];
        outFile << " " << laserTv << " " <<  laserRv << " " <<  forwardSafetyDist << " " << sideSaftyDist << " " << turnAxis;

        // timestamp + host
        outFile << " " << timestamp;
        outFile << " " << hostname;
        outFile << " " << loggerTimestamp << std::endl;
      }
    }
  outFile.close();
}


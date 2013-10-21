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

#include <stdlib.h>
#include <string.h>
#include <istream>
#include <fstream>
#include <sstream>
#include <list>
#include <cstring>
#include <limits>
#include <omp.h>
#include "gridlinetraversal.h"

#include "import/alu_posegraph_2d/import_posegraph2d.h"
#include "ssa/core/sparse_surface_adjustment.h"
#include "ssa/core/ssa_graph_2d.h"
#include "ssa/sensor_models/laser_sensor_model_2d.h"

using namespace std;

const char *message[]={
  "ssa_to_pgm: computes a frequency map out of a ssa file",
  "usage ssa_to_pgm [options] <ssa_file>",
  "options:",
  " -res         <meters>                : the resolution of a grid cell",
  "                                        default [0.1m]",
  " -maxRange    <meters>                : the maximum range of the scanner",
  " -usableRange <meters>                : the maximum usable range of the scanner",
  0
};

bool isInside(Eigen::Vector2i& point, Eigen::Vector2i& size){
  return (point(0) >= 0 && point(0) < size(0) && point(1) >= 0 && point(1) < size(1));
}

Eigen::Vector2i world2map(Eigen::Vector2f& endpoint, Eigen::Vector2f& offset, double& resolution){
  return Eigen::Vector2i(lrint((endpoint(0)-offset(0))/resolution),
		    lrint((endpoint(1)-offset(1))/resolution));
}

int main (int argc, const char ** argv){

  if (argc<2){
    const char**v=message;
    while (*v){
      cout << *v << endl;
      v++;
    }
    return 0;
  }
  double maxRange=80;
  double usableRange=49;
  double resolution=0.1;
  const char* logfile=0;
  int c=1;
  while (c<argc){
    if (!strcmp(argv[c],"-res")){
      c++;
      resolution=atof(argv[c]);
      c++;
    } else
    if (!strcmp(argv[c],"-maxrange")){
      c++;
      maxRange=atof(argv[c]);
      c++;
    } else
    if (!strcmp(argv[c],"-usablerange")){
      c++;
      usableRange=atof(argv[c]);
      c++;
    } else
    if (! logfile){
      logfile=argv[c];
      c++;
      break;
    }
  }

  /** Load ssa graph */
  ssa::SparseSurfaceAdjustmentT<g2o::EdgeSE2, ssa::EdgeSE2PointXYCov, ssa::EdgePointXYCovPointXYCov> ssa;
  ssa.setVerbose(true);

  if (! logfile){
    cerr << "Missing ssa graph file." << endl;
    exit(-1);
  } else {
    ssa.graph()->load(logfile);
  }

  ///compute bbox
  double minX = 1e10;
  double maxX = -1e10;
  double minY = 1e10;
  double maxY = -1e10;
  for (int i=0; i < (int)ssa.graph()->_verticies_points.size(); ++i){
    minX = std::min(ssa.graph()->_verticies_points[i]->estimate()(0), minX);
    maxX = std::max(ssa.graph()->_verticies_points[i]->estimate()(0), maxX);
    minY = std::min(ssa.graph()->_verticies_points[i]->estimate()(1), minY);
    maxY = std::max(ssa.graph()->_verticies_points[i]->estimate()(1), maxY);
  }
  ///setting size and offset
  Eigen::Vector2f size( ceil(maxX - minX) + 1.0, ceil(maxY - minY) + 1.0);
  Eigen::Vector2f offset(minX - 0.5, minY - 0.5);
  Eigen::Vector2i isize((int)round(size(0)/resolution), (int)round(size(1)/resolution));

  ///initialize map
  std::vector< std::pair<int, int> > frequencyMap(isize(0) * isize(1)); ///pair first -> hits, second -> misses
  for (int y=1; y<=(int)frequencyMap.size(); y++)
  {
    frequencyMap[y].first = 0;
    frequencyMap[y].second = 0;
  }

  if (usableRange<0)
    usableRange=maxRange;

  for (g2o::OptimizableGraph::EdgeSet::iterator it=ssa.graph()->_optimizer.edges().begin(); it!=ssa.graph()->_optimizer.edges().end(); it++){
    ssa::EdgeSE2PointXYCov* edge=dynamic_cast<ssa::EdgeSE2PointXYCov*>(*it);
    if(edge){
      g2o::VertexSE2*         vp1=dynamic_cast<g2o::VertexSE2*>(edge->vertices()[0]);
      ssa::VertexPointXYCov*  vp2=dynamic_cast<ssa::VertexPointXYCov*>(edge->vertices()[1]);
      if(vp1 && vp2){
        Eigen::Vector2f rp= Eigen::Vector2f(vp1->estimate().translation()(0), vp1->estimate().translation()(1));
        Eigen::Vector2i start= world2map(rp, offset, resolution);
        double r=edge->measurement().norm();
        if (r>=maxRange)
          continue;

        bool cropped=false;
        if (r>usableRange){
          r=usableRange;
          cropped=true;
        }

        static ssa::GridLineTraversalLine line;
        Eigen::Vector2f bp(vp2->estimate()(0),  vp2->estimate()(1));
        //Eigen::Vector2f bp((vp1->estimate() * edge->measurement())(0), (vp1->estimate() * edge->measurement())(1));
        Eigen::Vector2i end = world2map(bp, offset, resolution);
        ssa::GridLineTraversal::gridLine(start, end, &line);

        for (int i=0; i<line.num_points; i++){
          if( isInside(line.points[i], isize) ){
            frequencyMap[line.points[i](1) * isize(0) + line.points[i](0)].second+=1;
          }
        }
        if (! isInside(end, isize)){
          continue;
        }
        if (! cropped){
	  frequencyMap[end(1) * isize(0) + end(0)].first+=1;
        }
      }
    }
  }

  string ppmName = logfile;
  ppmName = ppmName  + ".pgm";
  ofstream ppmFile(ppmName.c_str());

  ppmFile << "P6" << endl;
  ppmFile << "#resolution " << resolution << endl;
  ppmFile << "#offset "     << offset(0) << " " << offset(1) << endl;
  ppmFile << isize(0) << " " << isize(1) << endl  << 255 << endl;

  int height=isize(1);
  int width=isize(0);
  float max=1.;
  for (int i=1; i<=height; i++){
    for (int x=0; x<width; x++){
      int y=height-i;
      float fo= 0.0f;
      if(frequencyMap[y*width+x].first == 0 && frequencyMap[y*width+x].second == 0){
	fo=-1.;
      } else {
        if(frequencyMap[y*width+x].second > 0)
	  fo= (float) frequencyMap[y*width+x].first / (float) frequencyMap[y*width+x].second;
	if(frequencyMap[y*width+x].first > 0 && frequencyMap[y*width+x].second == 0)
	  fo= 1.0f;
      }

      float occ=fo*max;
      unsigned char c=(unsigned char)(255.-255*occ);
      unsigned char r=c, g=c, b=c;
      if (fo==-1.) { // unknown
	b=(unsigned char)(210);
	g=(unsigned char)(190);
	r=(unsigned char)(190);
      } else if (fo==-2.) {// red
	b=(unsigned char)(64);
	g=(unsigned char)(64);
	r=(unsigned char)(255);
	//cerr << "r";
      } else if (fo==-3.){
	b=(unsigned char)(255);
	g=(unsigned char)(64);
	r=(unsigned char)(64);
	//cerr << "b";
      }
      ppmFile.put(r);
      ppmFile.put(g);
      ppmFile.put(b);
    }
  }
  ppmFile.close();

}

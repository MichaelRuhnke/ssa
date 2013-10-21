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

#include "Eigen/Core"

#include "pcl/io/pcd_io.h"
#include "pcl/point_types.h"
#include "pcl/point_cloud.h"
#include "pcl/common/transforms.h"
#include "pcl/common/eigen.h"

#include "pcl_ssa_hierarchical.h"
using namespace std;
using namespace ssa;

const char *message[]={
  "ssa_to_pcd: converts a ssa3d graph into one pcd point cloud file",
  "usage ssa_convert_to_pcd: [options] <ssa3d_file> <pcd_file>",
  "options:",
  "-p [string]       set prefix and enable dumping of single view point clouds to disk.",
  0
};


int main(int argc, char **argv)
{

 if (argc<3){
    const char**v=message;
    while (*v){
      cout << *v << endl;
      v++;
    }
    return 0;
  }
  const char* ssaFile=0;
  const char* outputFile=0;
  string pcPrefix = "ssa_model";
  bool saveViews = false;
  bool reproject = false;
  int c=1;
  while (c<argc){
    if (!strcmp(argv[c],"-p")){
      c++;
      pcPrefix=argv[c];
      c++;
      saveViews = true;
    }
    if (!strcmp(argv[c],"-r")){
      reproject=true;
      c++;
    }
    if (! ssaFile){
      ssaFile=argv[c];
      c++;
    }
    if (! outputFile){
      outputFile=argv[c];
      c++;
      break;
    }
  }
  cerr << ssaFile << " -> " << outputFile << endl;

  //Create ssa 3d graph
  SparseSurfaceAdjustmentGraph3D ssaGraph;

  //Load SSA3D graph
  ssaGraph.load(ssaFile);
  pcl::PointCloud<pcl::PointXYZRGBA> cloud;
  if(reproject){
    cloud = PCLSSAHierarchicalT<pcl::PointCloud<pcl::PointXYZRGBA> >::landmarksToPointCloud(ssaGraph._edges_observations, reproject);
  } else {
    int i = 0;
    std::map<int, int> keys = ssaGraph.getScanIds();
    for(std::map<int, int>::iterator it = keys.begin(); it != keys.end(); ++it){
      int scanId = it->first;
      pcl::PointCloud<pcl::PointXYZRGBA> tmp = PCLSSAHierarchicalT<pcl::PointCloud<pcl::PointXYZRGBA> >::landmarksToPointCloud(ssaGraph._verticies_observations[scanId]);

      //copy points to global point cloud
      std::copy(tmp.points.begin(), tmp.points.end(), back_inserter(cloud.points));

      if(saveViews){
        //transform points in local coordinate frame (scan)
        g2o::VertexSE3* v = dynamic_cast<g2o::VertexSE3*>(ssaGraph._optimizer.vertex(scanId));
        Eigen::Affine3f transformation = v->estimate().cast<float>();
        pcl::transformPointCloud(tmp, tmp, transformation.inverse());

        float x, y, z, roll, pitch, yaw;
        pcl::getTranslationAndEulerAngles(transformation, x, y, z, roll, pitch, yaw);

        tmp.sensor_origin_(0) = x;
        tmp.sensor_origin_(1) = y;
        tmp.sensor_origin_(2) = z;
        tmp.sensor_origin_(3) = 1;
        tmp.sensor_orientation_ = transformation.rotation();

        char filename[128];
        sprintf(filename, "%s%05d.pcd", pcPrefix.c_str(), i);
        cerr << "saving file:"  << filename << endl;
        pcl::io::savePCDFileBinary(filename, tmp);
        i++;
      }
    }
  }
  pcl::io::savePCDFileBinary(outputFile , cloud);
}


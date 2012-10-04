#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <istream>
#include <fstream>
#include <sstream>
#include <list>
#include <cstring>
#include <limits>
// #include <cmath>

#include <QApplication>

#include "Eigen/Core"

#include "pcl/io/pcd_io.h"
#include "pcl/io/ply_io.h"
#include "pcl/point_types.h"
#include "pcl/point_cloud.h"
#include "pcl/common/transforms.h"

#include "pcl/filters/filter.h"
#include "pcl/filters/voxel_grid.h"

#include "g2o/stuff/timeutil.h"
#include "pcl_ssa_hierarchical.h"

using namespace std;
using namespace ssa;

const char *message[]={
  "ply_to_ssa_rgb: converts a set of pcd files into a ssa 3d graph file",
  "usage ply_to_ssa_rga [options] <ssa3d_file>",
  "options:",
  "-e [int]            take every e.th scan.",
  "-n [int]            number of scans.",
  "-p [string]         prefix of the ply file list.",
  "-lms                use lms sensor model instead of rgbd sensor model.",
  "-r [double]         resolution of the resulting model (default resolution 0.01m)",
  "example:",
  "ply_to_ssa_rgb -p alufr_black_mug_raw -n 60 -r 0.001 alufr_black_mug_raw.ssa3d",
  0
};


int main(int argc, char **argv)
{

 //Debuging stuff
 //OptimizableGraph::printRegisteredTypes(cout, true);
  string pcPrefix  = "pcloud";
 if (argc<2){
    const char**v=message;
    while (*v){
      cout << *v << endl;
      v++;
    }
    return 0; 
  }

  const char* outputFile=0;

  int numOfScans = 100;
  int every = 1;
  double voxelSize = 0.01;
  bool useLMS = false;
  int c=1;
  while (c<argc){
    if (!strcmp(argv[c],"-e")){
      c++;
      every=atoi(argv[c]);
      c++;
    } else 
    if (!strcmp(argv[c],"-n")){
      c++;
      numOfScans=atoi(argv[c]);
      c++;
    } else 
    if (!strcmp(argv[c],"-p")){
      c++;
      pcPrefix=argv[c];
      c++;
    } else 
    if (!strcmp(argv[c],"-r")){
      c++;
      voxelSize=atof(argv[c]);
      c++;
    } else
    if (!strcmp(argv[c],"-lms")){
      useLMS = true;
      c++;
    } else
    if (! outputFile){
      outputFile=argv[c];
      c++;
      break;
    }
  }
  //Create ssa 3d graph
  SparseSurfaceAdjustmentGraph3D ssaGraph;

  vector< pcl::PointCloud<pcl::PointXYZRGB>::Ptr > clouds;
  cerr << "Building ssa graph with resolution " << voxelSize << "m."<< endl;

  PCLSSAHierarchicalT<pcl::PointCloud<pcl::PointXYZRGB> >  pclToSSA;
  //pclToSSA.setInput(&clouds);
  pclToSSA.setSensor(PCLSSAHierarchicalT<pcl::PointCloud<pcl::PointXYZRGB> >::KINECT);
  pclToSSA.setLevels(3);
  pclToSSA.setResolution(0, 0.001);
  pclToSSA.setResolution(1, 0.005);
  pclToSSA.setResolution(2, 0.01);
  pclToSSA.setResolution(3, 0.02);
  pclToSSA.setResolution(4, 0.05);
  pclToSSA.params().normalExtractionMaxNeighborDistance = 0.2;
  pclToSSA.params().normalExtractionMinNeighbors = 9;
  //pclToSSA.constructGraph(ssaGraph);

  for(int i = 0; i < numOfScans; i=i+every)
  {
    double timing = get_time();
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    char cloudFileName[2048];
    //sprintf(cloudFileName, "%s%03d.pcd", pcPrefix.c_str(), i);
    sprintf(cloudFileName, "%s%d.ply", pcPrefix.c_str(), i);
    pcl::io::loadPLYFile(cloudFileName, *cloud);
    cerr << "loading from disk took " << (get_time() - timing) * 1000 << "ms." << endl;
    pclToSSA.addCloudToGraph(ssaGraph, cloud);
  }

  cout << "Writing ssa3d file to " << outputFile << endl;
  ssaGraph.save(outputFile);
}


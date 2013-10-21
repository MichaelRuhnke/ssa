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
  "pcd_to_ssa: converts a set of pcd files into a ssa 3d graph file",
  "usage pcd_to_ssa [options] <ssa3d_file>",
  "options:",
  "-e [int]            take every e.th scan.",
  "-n [int]            number of scans.",
  "-p [string]         prefix of the pcd file list.",
  "-lms                use lms sensor model instead of kinect.",
  "-r [double]         resolution of the resulting model (default resolution 0.01m)",
  "example:",
  "pcd_to_ssa -p alufr_black_mug_raw -n 60 -r 0.001 alufr_black_mug_raw.ssa3d",
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

  vector< pcl::PointCloud<pcl::PointXYZRGBA>::Ptr > clouds;
  cerr << "Building ssa graph with resolution " << voxelSize << "m."<< endl;

  PCLSSAHierarchicalT<pcl::PointCloud<pcl::PointXYZRGBA> >  pclToSSA;
  //pclToSSA.setInput(&clouds);
  if(useLMS){
    pclToSSA.setSensor(PCLSSAHierarchicalT<pcl::PointCloud<pcl::PointXYZRGBA> >::LMS);
    pclToSSA.setLevels(4);
    pclToSSA.setResolution(0, 0.001);
    pclToSSA.setResolution(1, 0.01);
    pclToSSA.setResolution(2, 0.1);
    pclToSSA.setResolution(3, 0.2);
    pclToSSA.params().normalExtractionMaxNeighborDistance = 1.0;
  } else {
    pclToSSA.setSensor(PCLSSAHierarchicalT<pcl::PointCloud<pcl::PointXYZRGBA> >::KINECT);
    pclToSSA.setLevels(6);
    pclToSSA.setResolution(0, 0.001);
    pclToSSA.setResolution(1, 0.005);
    pclToSSA.setResolution(2, 0.01);
    pclToSSA.setResolution(3, 0.02);
    pclToSSA.setResolution(4, 0.1);
    pclToSSA.setResolution(5, 0.4);
  pclToSSA.params().normalExtractionMaxNeighborDistance = 0.5;
  }
  pclToSSA.params().normalExtractionMinNeighbors = 9;

  for(int i = 0; i < numOfScans; i=i+every)
  {
    double timing = get_time();
    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGBA>);
    char cloudFileName[2048];
    sprintf(cloudFileName, "%s%05d.pcd", pcPrefix.c_str(), i);
    pcl::io::loadPCDFile(cloudFileName , *cloud);
    cerr << "loading from disk took " << (get_time() - timing) * 1000 << "ms." << endl;
    pclToSSA.addCloudToGraph(ssaGraph, cloud);
  }

  cout << "Writing ssa3d file to " << outputFile << endl;
  ssaGraph.save(outputFile);
}


/* \author Michael Ruhnke, Bastian Steder*/

#include <iostream>
#include <fstream>
using namespace std;

#include <boost/filesystem.hpp>

#include "pcl/range_image/range_image_planar.h"
#include "pcl/common/common_headers.h"
#include "pcl/io/pcd_io.h"
#include "pcl/console/parse.h"
#include <pcl/common/transformation_from_correspondences.h>
#include <pcl/keypoints/sift_keypoint.h>
#include <pcl/features/range_image_border_extractor.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/passthrough.h>
#include <pcl/features/normal_3d_omp.h>
#include "pcl/common/transforms.h"

using namespace pcl;

void
printUsage (const char* progName)
{
  cout << "\n\nUsage: "<<progName<<" [options]\n\n"
       << "Options:\n"
       << "-------------------------------------------\n"
       << "-p               point cloud prefix\n"
       << "-m <int>         maximum number of scans that should be processed\n"
       << "-e <int>         process only every e-th scan\n"
       << "-h              this help\n"
       << "\n\n";
}

int main (int argc, char** argv)
{
  string pcPrefix  = "pcloud";
  string ocPrefix  = "ncloud";
  // --------------------------------------
  // -----Parse Command Line Arguments-----
  // --------------------------------------
  if (pcl::console::find_argument (argc, argv, "-h") >= 0)
  {
    printUsage (argv[0]);
    return 0;
  }

  pcl::console::parse (argc, argv, "-p", pcPrefix);

  int  maxStep = 1;
  pcl::console::parse (argc, argv, "-m", maxStep);

  int  skipScansPerStep = 1;
  pcl::console::parse (argc, argv, "-e", skipScansPerStep);


  for(int i = 0; i < maxStep; i=i+skipScansPerStep)
  {
    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGBA>);
    char cloudFileName[2048];
    sprintf(cloudFileName, "%s%05d.pcd", pcPrefix.c_str(), i);
    pcl::io::loadPCDFile(cloudFileName , *cloud);
    if(cloud->size() == 0)
      continue;
    std::cerr << ".";

    pcl::NormalEstimationOMP<pcl::PointXYZRGBA, pcl::Normal> ne;
    pcl::search::KdTree<pcl::PointXYZRGBA>::Ptr tree1 (new pcl::search::KdTree<pcl::PointXYZRGBA>);
    tree1->setInputCloud (cloud);
    ne.setInputCloud (cloud);
    ne.setSearchMethod (tree1);
    ne.setKSearch (200);
    ne.setViewPoint(0.0f, 0.0f, 0.0f);
    pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);
    ne.compute (*normals);


    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud_remapped(new pcl::PointCloud<pcl::PointXYZRGBNormal>(640,480));
    cloud_remapped->sensor_orientation_ = cloud->sensor_orientation_;
    cloud_remapped->sensor_origin_ = cloud->sensor_origin_;
    cloud_remapped->width = cloud->width;
    cloud_remapped->height = cloud->height;
    for (size_t j = 0; j < cloud->size (); j++)
    {
      const pcl::PointXYZRGBA& point = cloud->at(j);
      pcl::PointXYZRGBNormal& pt_remapped = cloud_remapped->at(j);
      pt_remapped.x = point.x;
      pt_remapped.y = point.y;
      pt_remapped.z = point.z;
      pt_remapped.r = point.r;
      pt_remapped.g = point.g;
      pt_remapped.b = point.b;
      pt_remapped.getNormalVector3fMap() = normals->at(j).getNormalVector3fMap();
    }

    sprintf(cloudFileName, "%s%05d.pcd", ocPrefix.c_str(), i);
    if(cloud_remapped->size() > 0)
      pcl::io::savePCDFileBinary (cloudFileName,  *cloud_remapped);
  }
}


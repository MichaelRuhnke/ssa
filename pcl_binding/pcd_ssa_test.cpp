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
#include <pcl/features/normal_3d_omp.h>
#include <pcl/common/time.h>

#include <pcl/search/organized.h>

#include "g2o/stuff/timeutil.h"
#include "g2o/core/solver.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/optimization_algorithm_levenberg.h"

#include "ssa/core/allocate_solver.h"
#include "ssa/types_3d/edge_se3_point_xyz_normal.h"
#include "pcl_ssa_hierarchical.h"
#include "ssa/data_association/data_association.h"

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

float maxRange = 4.0f;
float maxAngle = g2o::deg2rad(30);

struct Correspondence{
  size_t scanId;
  float  x;
  float  y;
  float  z;
};

typedef std::tr1::unordered_map<size_t, std::vector< Correspondence > >  DataAssociation;


void computeKeyFrameSkeleton(std::vector< pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr >& clouds, std::vector< Eigen::Affine3f >& poses, pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& skeleton,  std::vector<int>& keyframes, float& resolution){
  skeleton->points.clear();
  //pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr globalCloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  std::tr1::unordered_map< int , std::tr1::unordered_map< int,  std::tr1::unordered_map< int, std::vector< std::pair< int, int > > > > > virtualGrid;
//   std::vector<Eigen::Vector3i> indices;
  for(size_t i = 0; i < keyframes.size(); ++i){
    int& cloudId = keyframes[i];
    Eigen::Affine3f& pose = poses[cloudId];
    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& cloud = clouds[i];
    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud_transformed(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
    pcl::transformPointCloudWithNormals(*cloud, *cloud_transformed, pose);
    for(size_t j = 0; j < cloud_transformed->size(); ++j){
      pcl::PointXYZRGBNormal& point = cloud_transformed->points[j];
      if(!std::isnan(point.z) && cloud->points[j].z < maxRange){
        Eigen::Vector3i index;
        index(0) = lrint(point.x / resolution);
        index(1) = lrint(point.y / resolution);
        index(2) = lrint(point.z / resolution);
        if(virtualGrid[index(0)][index(1)][index(2)].size() == 0){
          virtualGrid[index(0)][index(1)][index(2)].push_back(std::pair<int, int>(cloudId, j));
          //indices.push_back(index);
          skeleton->points.push_back(point);
        }
      }
    }
  }
//   for(size_t j = 0; j < indices.size(); ++j){
//     Eigen::Vector3i& index = indices[j];
//     std::pair< int, int >& other = virtualGrid[index(0)][index(1)][index(2)][0];
//     pcl::PointXYZRGBNormal point = clouds[other.first]->points[other.second];
//     point.curvature = point.z;
//     Eigen::Affine3f& pose = poses[other.first];
//     point.getVector3fMap() = pose * point.getVector3fMap();
//     skeleton->points.push_back(point);
//   }
}

void solveDA(std::vector< pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr >& clouds, std::vector< Eigen::Affine3f >& poses, pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& skeleton, SparseSurfaceAdjustmentGraphT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZNormal, ssa::EdgePointXYZNormalPointXYZNormal>& ssaGraph, std::vector<VertexSE3*>& vertices, std::vector<VertexPointXYZNormal*>& pointVertices, float maxDistance){
  float angleNinty = g2o::deg2rad(90);
//   double timing = pcl::getTime();
  ssa::RGBDSensorModel  sensorModel;

  ssaGraph.dropSensorEdges();
  std::cerr << "check point vertices \t  ";
  ///check if vertices had been cleared
  for(size_t i = 0; i < pointVertices.size(); ++i){
    if(pointVertices[i] != 0){
      ssa::VertexPointXYZNormal* v =  pointVertices[i];
//       std::cerr << PVAR(i) << " " << pointVertices[i] << std::endl;
      pointVertices[i] = 0;
      ssaGraph.removeVertex(v);
    }
  }
  std::cerr << "[done]";

  ///solve data association
  for(size_t i = 0; i < clouds.size(); ++i){
    std::cerr << "processing da for cloud " << i << std::endl;
    Eigen::Affine3f& pose = poses[i];
    g2o::VertexSE3* vertex = vertices[i];

    pcl::search::KdTree<pcl::PointXYZRGBNormal>::Ptr kdTree (new pcl::search::KdTree<pcl::PointXYZRGBNormal>);
    kdTree->setInputCloud (clouds[i]);

    for(size_t j = 0; j < skeleton->size(); ++j){
      const pcl::PointXYZRGBNormal& point = skeleton->points[j];
      pcl::PointXYZRGBNormal  tmpPoint = point;
      tmpPoint.getVector3fMap() = pose.inverse() * skeleton->points[j].getVector3fMap(); ///local frame
      std::vector< int > k_indices;
      std::vector< float > k_sqr_distances;
      //kdTree->nearestKSearch(tmpPoint, 1, k_indices, k_sqr_distances);
      kdTree->radiusSearch(tmpPoint, maxDistance, k_indices, k_sqr_distances);
      if(k_indices.size() > 0){
        int pid = k_indices[0];
        ///check if skeleton point was already assigned
        const pcl::PointXYZRGBNormal& p = clouds[i]->points[pid];

        ///check direction of assignment accordign to normal
        Eigen::Vector3f diff = tmpPoint.getVector3fMap() - p.getVector3fMap();
        float angle = fabs(acos(diff.normalized().dot(p.getNormalVector3fMap().normalized())));
        if(angle > angleNinty)
          angle = (angleNinty + angleNinty) - angle;
//         std::cerr << PVAR(diff.transpose()) << " " << PVAR(g2o::rad2deg(angle)) << std::endl;
        if(angle > maxAngle)
          continue;

        Eigen::Vector3d measurement;
        measurement(0) = p.x;
        measurement(1) = p.y;
        measurement(2) = p.z;

        int pointToAssign = j;

        if(pointVertices[pointToAssign] == 0){
          ssa::VertexPointXYZNormal* v = new ssa::VertexPointXYZNormal();
          v->normal()(0) = p.normal_x;
          v->normal()(1) = p.normal_y;
          v->normal()(2) = p.normal_z;
          ///move point to world frame
          v->setEstimate(vertex->estimate() * measurement);
          v->setParentVertex(vertices[i]);
          v->cr = p.r;
          v->cg = p.g;
          v->cb = p.b;
          pointVertices[j] = v;
          ssaGraph.addVertex(v);
//           std::cerr << "v";
        }

        ssa::EdgeSE3PointXYZNormal* e = new ssa::EdgeSE3PointXYZNormal;
        e->vertices()[0] = vertices[i];
        e->vertices()[1] = pointVertices[pointToAssign];
        e->setMeasurement(measurement);
        e->setLevel(0);
        e->information().setIdentity();
        sensorModel.getInformationMatrix(measurement, e->information());
        ssaGraph.addEdge(e);
//         std::cerr << "e";
      }
    }
  }
}

int width = 640;
int height = 480;
float focal_length_x = 525.;
float focal_length_y = 525.;
float principal_point_x = 319.5;
float principal_point_y = 239.5;

///fr1 dataset
//float focal_length_x = 517.3f;
//float focal_length_y = 516.5f;
//float principal_point_x = 318.6;
//float principal_point_y = 255.3;

bool getProjectedPointIndex(pcl::PointXYZRGBNormal& point, Vector2i& index){
  if(pcl_isnan(point.z) || point.z == 0.0f)
    return false;
  index(0) = (point.x / point.z) * focal_length_x + principal_point_x;
  index(1) = (point.y / point.z) * focal_length_y + principal_point_y;
  return (index(0) >= 0 && index(0) < width && index(1) >= 0 && index(1) < height);
}

int searchProjectedNeighborhood(pcl::PointXYZRGBNormal& point, pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& cloud, Vector2i& index, float& maxDistance, int& width, int& height, float& distance){
  int searchRange = 2;
  float minDistance = 100 * maxDistance * maxDistance;
  //float minRGBDistance = 1000000;
//   std::tr1::unordered_map<int, bool> checked;
  bool foundMinimum = false;

  Eigen::Vector2i offset(0, 0);
  Eigen::Vector2i last_offset(0, 0);
  while(!foundMinimum){
    last_offset = offset;
    for(int dy = -searchRange; dy <= searchRange; ++dy){
      for(int dx = -searchRange; dx <= searchRange; ++dx){
        Eigen::Vector2i cell((index(0)+offset(0)+dx), (index(1)+offset(1)+dy));
        int refIndex = cell(1) * width + cell(0);
        if( cell(0) < 0 || cell(0) >= width || cell(1) < 0 || cell(1) >= height) //checked[refIndex] ||
          continue;

        const pcl::PointXYZRGBNormal& candidate = cloud->points[refIndex];
        if(isnan(candidate.z)){
//           checked[refIndex] = true;
          continue;
        }
        float d = (point.getVector3fMap() - candidate.getVector3fMap()). squaredNorm();
        //float rgbDistance = (point.getRGBVector3i() - cloud->points[refIndex].getRGBVector3i()).cast<float>().norm();
        if(d < minDistance){ //&& rgbDistance  < 1.1 * minRGBDistance
          minDistance = d;
          //minRGBDistance = rgbDistance;
          offset(0) = last_offset(0) + dx;
          offset(1) = last_offset(1) + dy;
        }
//         checked[refIndex] = true;
      }
    }

    if(last_offset == offset) ///break if no point in range found
      foundMinimum = true;
  }
  int result = -1;
  index += offset;
  distance = sqrt(minDistance);
  if(distance < maxDistance)
    result = index(1) * width + index(0);
  return result;
}

void solveDAProjective(std::vector< pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr >& clouds, std::vector< Eigen::Affine3f >& poses, pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& skeleton, SparseSurfaceAdjustmentGraphT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZNormal, ssa::EdgePointXYZNormalPointXYZNormal>& ssaGraph, std::vector<VertexSE3*>& vertices, std::vector<VertexPointXYZNormal*>& pointVertices, float maxDistance, int threads){
  float angleNinty = g2o::deg2rad(90);
//   double timing = pcl::getTime();
  ssa::RGBDSensorModel  sensorModel;

  ssaGraph.dropDataAssociation();
  ssaGraph.dropSensorEdges();
  std::cerr << "check point vertices \t  ";
  ///check if vertices had been cleared
  ssaGraph._verticies_points.clear();
  for(size_t i = 0; i < pointVertices.size(); ++i){
    if(pointVertices[i] != 0){
      ssa::VertexPointXYZNormal* v =  pointVertices[i];
      pointVertices[i] = 0;
      ssaGraph.optimizer()->removeVertex(v);
    }
  }
  std::cerr << "[done]" << std::endl;


  std::cerr << "processing da: "<< std::endl;
  ///solve data association
  for(size_t i = 0; i < clouds.size(); ++i){
    int width = clouds[i]->width;
    int height = clouds[i]->height;
    std::cerr <<  i << " ";
    Eigen::Affine3f& pose = poses[i];
    g2o::VertexSE3* vertex = vertices[i];
    Eigen::Affine3f inversePose = pose.inverse();
    int count = 0;

    std::vector< std::pair<int, float> > assignments(clouds[i]->size());
    for(size_t j = 0; j < clouds[i]->size(); ++j){
      assignments[j].first = 0;
      assignments[j].second = 1e10;
    }

    #pragma omp parallel for default(shared) schedule(dynamic, 2) num_threads(threads)
    for(size_t j = 0; j < skeleton->size(); ++j){
      pcl::PointXYZRGBNormal  tmpPoint = skeleton->points[j];
      tmpPoint.getVector3fMap() = inversePose * tmpPoint.getVector3fMap(); ///local frame

      Eigen::Vector2i proj(0,0);
      if(getProjectedPointIndex(tmpPoint, proj)){
//         int pid = searchProjectedNeighborhood(tmpPoint, clouds[i], proj, maxDistance);

        int pid = proj(1) * width + proj(0);
        if(pid < 0 || proj(0) < 2 || proj(0) > width - 2 || proj(1) < 2 || proj(1) > height - 2){
          continue;
        }

        ///check if skeleton point was already assigned
        const pcl::PointXYZRGBNormal& p = clouds[i]->points[pid];
        if(isnan(p.z)){
          continue;
        }

        ///check direction of assignment accordign to normal
        Eigen::Vector3f diff = tmpPoint.getVector3fMap() - p.getVector3fMap();
        float distance = diff.norm();
//         std::cerr << PVAR(diff.norm()) << std::endl;
        if(distance > 0.1f) //5.0 * sensorModel.depthUncertainty(p.z))
           continue;

        float angle = fabs(acos(diff.normalized().dot(p.getNormalVector3fMap().normalized())));
        if(angle > angleNinty)
          angle = (angleNinty + angleNinty) - angle;
//         std::cerr << PVAR(diff.norm()) << " " << PVAR(g2o::rad2deg(angle)) << std::endl;
//         if(angle > maxAngle)
//           continue;
        #pragma omp critical
        if(distance < assignments[pid].second)
        {
          assignments[pid].second = distance;
          assignments[pid].first = j;
        }
      }
    }

    for(size_t j = 0; j < assignments.size(); ++j){
      std::pair<int, float>& a = assignments[j];
      if(a.first == 0)
        continue;

      const pcl::PointXYZRGBNormal& p = clouds[i]->points[j];

      Eigen::Vector3d measurement;
      measurement(0) = p.x;
      measurement(1) = p.y;
      measurement(2) = p.z;
      int pointToAssign = a.first;

      if(pointVertices[pointToAssign] == 0){
        ssa::VertexPointXYZNormal* v = new ssa::VertexPointXYZNormal();
        v->normal()(0) = p.normal_x;
        v->normal()(1) = p.normal_y;
        v->normal()(2) = p.normal_z;
        ///move point to world frame
        v->setEstimate(vertex->estimate() * measurement);
        v->setParentVertex(vertex);
        v->cr = p.r;
        v->cg = p.g;
        v->cb = p.b;
        pointVertices[pointToAssign] = v;
        ssaGraph.addVertex(v);
      }


      ssa::EdgeSE3PointXYZNormal* e = new ssa::EdgeSE3PointXYZNormal;
      e->vertices()[0] = vertex;
      e->vertices()[1] = pointVertices[pointToAssign];
      e->setMeasurement(measurement);
      e->setLevel(0);
      e->information().setIdentity();
      sensorModel.getInformationMatrix(measurement, e->information());
      #pragma omp critical
      {
        ssaGraph.addEdge(e);
        count++;
      }
    }
    std::cerr << "assigned " << count  << " of " << clouds[i]->size() << " points to skeleton " << std::endl;
  }
  std::cerr << std::endl << "[done]"<< std::endl;
}


void solveDAProjectiveFull(std::vector< pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr >& clouds, std::vector< Eigen::Affine3f >& poses, SparseSurfaceAdjustmentGraphT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZNormal, ssa::EdgePointXYZNormalPointXYZNormal>& ssaGraph, std::vector<VertexSE3*>& vertices, float maxDistance, int threads){
  ssa::RGBDSensorModel  sensorModel;
  bool gicp = false;
  ssaGraph.dropDataAssociation();
  ssaGraph.dropSensorEdges();
  std::cerr << "check point vertices \t  ";
  ///check if vertices had been cleared

  for(size_t i = 0; i < ssaGraph._verticies_points.size(); ++i){
    if(ssaGraph._verticies_points[i] != 0){
      ssa::VertexPointXYZNormal* v =  ssaGraph._verticies_points[i];
      ssaGraph._verticies_points[i] = 0;
      ssaGraph.optimizer()->removeVertex(v);
    }
  }
  ssaGraph._verticies_points.clear();

  std::cerr << "[done]" << std::endl;
  std::tr1::unordered_map<int, std::tr1::unordered_map<int, int > > pointAssigned;

  std::cerr << "processing da: "<< std::endl;
  ///solve data association
  for(size_t i = 0; i < clouds.size(); ++i){
    int width = clouds[i]->width;
    int height = clouds[i]->height;
    std::cerr <<  i << " ";
    Eigen::Affine3f& pose = poses[i];
    g2o::VertexSE3* vertex = vertices[i];
    Eigen::Affine3f inversePose = pose.inverse();
    int count = 0;

    //
    #pragma omp parallel for default(shared) schedule(dynamic, 2) num_threads(threads)
    for(size_t j = 0; j < clouds[i]->size(); ++j){
      bool plop;
      #pragma omp critical
        plop = pointAssigned[i][j];

      if(plop)
        continue;

      pcl::PointXYZRGBNormal  referencePoint = clouds[i]->points[j];
      if(isnan(referencePoint.z) || referencePoint.z > maxRange)
        continue;

      std::vector< std::pair<int, int> > assignments;
/*      int y = floor(j / width);
      int x = j % width;    */
      for(size_t k = 0; k < clouds.size(); ++k){
        if(k == i)
          continue;
        referencePoint.getVector3fMap() = poses[k].inverse() * pose * clouds[i]->points[j].getVector3fMap(); ///transform to frame k

        Eigen::Vector2i proj(0,0);
        if(getProjectedPointIndex(referencePoint, proj)){
          float distance = 0.0f;;
          int pid = searchProjectedNeighborhood(referencePoint, clouds[k], proj, maxDistance, width, height, distance);
          //int pid = proj(1) * width + proj(0);

          if(pid < 0 || proj(0) < 2 || proj(0) > width - 2 || proj(1) < 2 || proj(1) > height - 2){
            continue;
          }

          ///check if point is valid
          pcl::PointXYZRGBNormal corPointK = clouds[k]->points[pid];
          if(isnan(corPointK.z)){
            continue;
          }

          corPointK.getVector3fMap() = inversePose * poses[k] * corPointK.getVector3fMap(); ///transform into frame i
          if(getProjectedPointIndex(corPointK, proj)){
            float tmpDistance = 0.0f;
            int bid = searchProjectedNeighborhood(corPointK, clouds[i], proj, maxDistance, width, height, tmpDistance);
            if(bid < 0 || proj(0) < 2 || proj(0) > width - 2 || proj(1) < 2 || proj(1) > height - 2 || fabs(distance - tmpDistance) > 0.001){
              continue;
            }
          }


          float angle = fabs(acos((poses[k].linear() * corPointK.getNormalVector3fMap()).normalized().dot((poses[i].linear() * referencePoint.getNormalVector3fMap()).normalized())));
          if(angle > maxAngle)
            continue;


          ///difference in global frame
          Eigen::Vector3f diff =  (poses[i] * clouds[i]->points[j].getVector3fMap()) - (poses[k] * clouds[k]->points[pid].getVector3fMap());
          if(diff.norm() > 10.0 * sensorModel.depthUncertainty(std::max(clouds[i]->points[j].z, clouds[k]->points[pid].z))) //5.0 * sensorModel.depthUncertainty(p.z))
             continue;
            assignments.push_back(std::pair<int, int>(k , pid));
        }
      }
//       std::cerr << assignments.size() << " ";
      if(assignments.size() == 0)
        continue;

      if(!gicp){
        const pcl::PointXYZRGBNormal& p = clouds[i]->points[j];
        Eigen::Vector3d measurement;
        measurement(0) = p.x;
        measurement(1) = p.y;
        measurement(2) = p.z;

        ssa::VertexPointXYZNormal* v = new ssa::VertexPointXYZNormal();
        v->normal()(0) = p.normal_x;
        v->normal()(1) = p.normal_y;
        v->normal()(2) = p.normal_z;
        ///move point to world frame
        v->setEstimate(vertex->estimate() * measurement);
        v->setParentVertex(vertex);
        v->cr = p.r;
        v->cg = p.g;
        v->cb = p.b;
        #pragma omp critical
          ssaGraph.addVertex(v);

        ssa::EdgeSE3PointXYZNormal* e = new ssa::EdgeSE3PointXYZNormal;
        e->vertices()[0] = vertex;
        e->vertices()[1] = v;
        e->setMeasurement(measurement);
        e->setLevel(0);
        e->information().setIdentity();
        sensorModel.getInformationMatrix(measurement, e->information());
        #pragma omp critical
        {
          ssaGraph.addEdge(e);
          count++;
          pointAssigned[i][j] = true;
        }

        for(size_t j = 0; j < assignments.size(); ++j){
          std::pair<int, int>& a = assignments[j];
          bool blib;
          #pragma omp critical
            blib = pointAssigned[a.first][a.second];
          if(blib)
            continue;
          const pcl::PointXYZRGBNormal& op = clouds[a.first]->points[a.second];

          Eigen::Vector3d m;
          m(0) = op.x;
          m(1) = op.y;
          m(2) = op.z;

          ssa::EdgeSE3PointXYZNormal* e = new ssa::EdgeSE3PointXYZNormal;
          e->vertices()[0] = vertices[a.first];
          e->vertices()[1] = v;
          e->setMeasurement(m);
          e->setLevel(0);
          sensorModel.getInformationMatrix(m, e->information());
          #pragma omp critical
          {
            ssaGraph.addEdge(e);
            count++;
            pointAssigned[a.first][a.second] = true;
          }
        }
      } else {

         const pcl::PointXYZRGBNormal& p1 = clouds[i]->points[j];
         pointAssigned[i][j] = true;

         for(size_t l = 0; l < assignments.size(); ++l){
           std::pair<int, int>& a = assignments[l];
//             bool blib;
//             #pragma omp critical
//               blib = pointAssigned[a.first][a.second];
//             if(blib)
//               continue;

           const pcl::PointXYZRGBNormal& p2 = clouds[a.first]->points[a.second];

           g2o::Edge_V_V_GICP * e           // new edge with correct cohort for caching
              = new g2o::Edge_V_V_GICP();

          e->vertices()[0]            // first viewpoint
            = dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex);

          e->vertices()[1]            // second viewpoint
            = dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertices[a.first]);

          g2o::EdgeGICP meas;
          meas.pos0 = p1.getVector3fMap().cast<double>();
          meas.normal0 = Eigen::Vector3d(1.0, 0.0, 0.0); //p1.getNormalVector3fMap().cast<double>();
          meas.pos1 = p2.getVector3fMap().cast<double>();
          meas.normal1 = Eigen::Vector3d(1.0, 0.0, 0.0); //p2.getNormalVector3fMap().cast<double>();

          e->setMeasurement(meas);
          meas = e->measurement();

          // use this for point-plane
          e->information() = meas.prec0(0.01);
          e->information().setIdentity();
          #pragma omp critical
          {
//             std::cerr << e->information() << std::endl << std::endl;
            count++;
            ssaGraph.addEdge(e);
            pointAssigned[a.first][a.second] = true;
          }
        }
      }
    }
    std::cerr << "assigned " << count  << " of " << clouds[i]->size() << " points to skeleton " << std::endl;
  }
  std::cerr << std::endl << "[done]"<< std::endl;
}

bool isValidCorrespondence(pcl::PointXYZRGBNormal& p1, pcl::PointXYZRGBNormal& p2, Eigen::Affine3f& trans1 , Eigen::Affine3f& trans2, float& maxRange, float& maxAngle, float maxSqrDistance){
  float angle = fabs(acos((trans2.linear() * p2.getNormalVector3fMap().normalized()).dot((trans1.linear() * p1.getNormalVector3fMap().normalized()))));
  float distance = ((trans2 * p2.getVector3fMap()) - (trans1 * p1.getVector3fMap())).squaredNorm();
  return (p2.z < maxRange && angle < maxAngle && distance < maxSqrDistance);
}

void solveProjectiveDA(
                      pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& cloudA,
                      size_t& idA,
                      pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& cloudB,
                      size_t& idB,
                      std::vector< Eigen::Affine3f >& poses,
                      std::vector<size_t>& validIndices,
                      DataAssociation&  dataAssociation,
                      std::vector<bool>& pointAssigned,
                      std::vector<pcl::PointXYZRGBNormal>& refPoints,
                      size_t& maxCloudSize,
                      int threads){

    ssa::RGBDSensorModel  sensorModel;
    Affine3f transform = poses[idB].inverse() * poses[idA];
    size_t referenceId = idA * maxCloudSize; ///Assuming all clouds have same dimensions!!!!
    //#pragma omp parallel for default(shared) schedule(dynamic, 2) num_threads(threads)
    for(size_t j = 0; j < cloudA->size(); ++j){
      size_t index = referenceId + j;
      if(pointAssigned[index])
        continue;

      pcl::PointXYZRGBNormal  referencePoint = cloudA->points[j];
      if(isnan(referencePoint.z) || referencePoint.z > maxRange)
        continue;

      referencePoint.getVector3fMap() = transform * cloudA->points[j].getVector3fMap(); ///transform to frame k

      Eigen::Vector2i proj(0,0);
      if(getProjectedPointIndex(referencePoint, proj)){
        int pid = proj(1) * width + proj(0);

        if(pid < 0 || proj(0) < 2 || proj(0) > width - 2 || proj(1) < 2 || proj(1) > height - 2){
          continue;
        }

        ///check if point is valid
        pcl::PointXYZRGBNormal corPointK = cloudB->points[pid];
        if(isnan(corPointK.z)){
          continue;
        }
        size_t corresId = idB * maxCloudSize + pid;
        if(pointAssigned[corresId]) ///avoid double information
          continue;

        float maxErrorDistance = 8.0 * sensorModel.depthUncertainty(std::max(cloudA->points[j].z, cloudB->points[pid].z));
        maxErrorDistance *= maxErrorDistance;
        if(isValidCorrespondence(cloudA->points[j], corPointK, poses[idA], poses[idB], maxRange, maxAngle, maxErrorDistance)){
          #pragma omp critical
          {
            if(dataAssociation[index].size() == 0){
              validIndices.push_back(index);
              refPoints.push_back(cloudA->points[j]);
            }
            if(dataAssociation[index].capacity() == dataAssociation[index].size() + 1)
              dataAssociation[index].reserve(dataAssociation[index].size() + 2);
            Correspondence cor;
            cor.scanId = idB;
            cor.x = corPointK.x;
            cor.y = corPointK.y;
            cor.z = corPointK.z;
            dataAssociation[index].push_back(cor);
            pointAssigned[corresId] = true;
          }
        }
      }
    }
}

void loadPointCloudBlock(std::vector< std::string >& cloudNames,
                         std::vector< pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr >& block,
                         std::vector< size_t >& blockIds,
                         size_t& startIndex, size_t& maxBlockSize){
   block.clear();
   blockIds.clear();
   for(size_t i = startIndex; i < cloudNames.size(); ++i){
    if(blockIds.size() > maxBlockSize){
      continue;
    }
    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
    pcl::io::loadPCDFile(cloudNames[i].c_str() , *cloud);
    block.push_back(cloud);
    blockIds.push_back(i);
  }
}

void solveDAOutOfCore(std::vector< std::string >& cloudNames,
                      std::vector< Eigen::Affine3f >& poses,
                      std::vector<size_t>& validIndices,
                      DataAssociation&  dataAssociation,
                      std::vector<pcl::PointXYZRGBNormal>& refPoints,
                      std::vector<bool>& pointAssigned,
                      size_t& solveBlockX,
                      size_t& maxBlockSize, size_t& maxCloudSize, int threads){
  std::cerr << "Solving full data association problem for every scan pair... this might take a while..." << std::endl;
  size_t blockCount = std::ceil(cloudNames.size() / maxBlockSize);

  std::vector< pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr > blockA;
  std::vector< size_t > blockAIds;
  std::vector< pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr > blockB;
  std::vector< size_t > blockBIds;

  size_t x = solveBlockX;
  size_t blockStartIndex = x * maxBlockSize;
  loadPointCloudBlock(cloudNames, blockA, blockAIds, blockStartIndex, maxBlockSize);
  ///solve da inside blockA
  for(size_t i = 0; i < blockA.size(); ++i){
    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& cloud = blockA[i];
    for(size_t j = i+1; j < blockA.size(); ++j){
      pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& cloud2 = blockA[j];
      solveProjectiveDA(cloud, blockAIds[i], cloud2, blockAIds[j], poses, validIndices, dataAssociation, pointAssigned, refPoints, maxCloudSize, threads);
    }
    std::cerr << "." << validIndices.size() << std::endl;
  }
  for(size_t y = x+1; y < blockCount; ++y){
    blockStartIndex = y * maxBlockSize;
    loadPointCloudBlock(cloudNames, blockB, blockBIds, blockStartIndex, maxBlockSize);
    for(size_t i = 0; i < blockA.size(); ++i){
      pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& cloud = blockA[i];
      for(size_t j = 0; j < blockB.size(); ++j){
        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& cloud2 = blockB[j];
        solveProjectiveDA(cloud, blockAIds[i], cloud2, blockBIds[j], poses, validIndices, dataAssociation, pointAssigned, refPoints, maxCloudSize, threads);
      }
      std::cerr << ".";
    }
    std::cerr << " " << validIndices.size() << std::endl;
  }
  std::cerr << std::endl;
}

void buildGraphForNextNEdges(DataAssociation&  dataAssociation,
                             std::vector<size_t>& validIndices,
                             std::vector<pcl::PointXYZRGBNormal>& refPoints,
                             std::vector< Eigen::Affine3f >& poses,
                             SparseSurfaceAdjustmentGraphT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZNormal, ssa::EdgePointXYZNormalPointXYZNormal>& ssaGraph,
                             std::vector<VertexSE3*>& vertices,
                             size_t cloudSize, int maxEdges){
  ssa::RGBDSensorModel  sensorModel;
  ///cleanup graph
  ssaGraph.dropDataAssociation();
  ssaGraph.dropSensorEdges();

  for(size_t i = 0; i < ssaGraph._verticies_points.size(); ++i){
    if(ssaGraph._verticies_points[i] != 0){
      ssa::VertexPointXYZNormal* v =  ssaGraph._verticies_points[i];
      ssaGraph._verticies_points[i] = 0;
      ssaGraph.optimizer()->removeVertex(v);
    }
  }
  ssaGraph._verticies_points.clear();

  int count = 0;
  int countAssignedIndices = 0;
  std::cerr << "Collecting up to " << maxEdges << " correspondences and building graph...  \t";
  for(size_t j = 0; j < validIndices.size(); ++j){
     std::vector< Correspondence >& assignments = dataAssociation[validIndices[j]];
     if(assignments.size() == 0)
       continue;
     size_t scanId = floor(validIndices[j] / cloudSize);

     if(scanId > vertices.size())
       std::cerr << "Error in data association result, vertex for " << scanId << " does not exist. " << vertices.size() << std::endl;
     g2o::VertexSE3*& vertex = vertices[scanId];

     const pcl::PointXYZRGBNormal& p = refPoints[j];
      Eigen::Vector3d measurement;
      measurement(0) = p.x;
      measurement(1) = p.y;
      measurement(2) = p.z;

      ssa::VertexPointXYZNormal* v = new ssa::VertexPointXYZNormal();
      v->normal()(0) = p.normal_x;
      v->normal()(1) = p.normal_y;
      v->normal()(2) = p.normal_z;
      ///move point to world frame
      v->setEstimate(vertex->estimate() * measurement);
      v->setParentVertex(vertex);
      v->cr = p.r;
      v->cg = p.g;
      v->cb = p.b;
      #pragma omp critical
        ssaGraph.addVertex(v);

      ssa::EdgeSE3PointXYZNormal* e = new ssa::EdgeSE3PointXYZNormal;
      e->vertices()[0] = vertex;
      e->vertices()[1] = v;
      e->setMeasurement(measurement);
      e->setLevel(0);
      e->information().setIdentity();
      sensorModel.getInformationMatrix(measurement, e->information());
      #pragma omp critical
      {
        ssaGraph.addEdge(e);
        count++;
      }
      for(size_t k = 0; k < assignments.size(); ++k){
        Correspondence& a = assignments[k];

        Eigen::Vector3d m;
        m(0) = a.x;
        m(1) = a.y;
        m(2) = a.z;

        ssa::EdgeSE3PointXYZNormal* e = new ssa::EdgeSE3PointXYZNormal;
        e->vertices()[0] = vertices[a.scanId];
        e->vertices()[1] = v;
        e->setMeasurement(m);
        e->setLevel(0);
        sensorModel.getInformationMatrix(m, e->information());
        #pragma omp critical
        {
          ssaGraph.addEdge(e);
          count++;
        }
      }
      std::vector< Correspondence > tmp;
      std::swap(assignments, tmp);
      validIndices[j] = 0;
      countAssignedIndices++;
      if(count > maxEdges)
        break;
    }
  std::vector< size_t > indices;
  std::vector< pcl::PointXYZRGBNormal > points;
  indices.reserve(validIndices.size() - countAssignedIndices);
  points.reserve(validIndices.size() - countAssignedIndices);
  for(size_t j = 0; j < validIndices.size(); ++j){
    if(validIndices[j] > 0){
      indices.push_back(validIndices[j]);
      points.push_back(refPoints[j]);
    }
  }
  std::swap(indices, validIndices);
  std::swap(points, refPoints);

  std::cerr << "[done]"<< std::endl;
}



size_t getKeyFrames(std::vector<std::string>& cloudNames, std::vector< pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr >& keyClouds, std::vector< Eigen::Affine3f >& poses, std::vector<int>& keyFrames, float maxDistance, float overlap){
  std::cerr << "selecting keyFrames: "<< std::endl;
  float maxSqrDistance = maxDistance * maxDistance;
  poses.resize(cloudNames.size());

  PCLSSAHierarchicalT<pcl::PointCloud<pcl::PointXYZRGBNormal> >  pclToSSA;
  size_t maxCloudSize = 0;
  ///solve data association
  for(size_t i = 0; i < cloudNames.size(); ++i){
    std::cerr << i << " ";
    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud2(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
    pcl::io::loadPCDFile(cloudNames[i].c_str() , *cloud);
    poses[i] = pclToSSA.getPose(cloud);

    if(cloud->size() == 0)
      continue;
    maxCloudSize = std::max(maxCloudSize, cloud->size());
    keyClouds.push_back(cloud);
    keyFrames.push_back(i);
    int width = cloud->width;
    int height = cloud->height;

    Eigen::Affine3f& pose = poses[i];
//     Eigen::Affine3f inversePose = pose.inverse();
    for(size_t k = i+1; k < cloudNames.size(); ++k){
      pcl::io::loadPCDFile(cloudNames[k].c_str() , *cloud2);
      poses[k] = pclToSSA.getPose(cloud2);
      if(cloud2->size() == 0)
        continue;

      int validPoints = 0;
      int assignedPoints = 0;
      #pragma omp parallel for default(shared) schedule(dynamic, 1) num_threads(4) reduction(+: validPoints) reduction(+: assignedPoints)
      for(size_t j = 0; j < cloud->size(); ++j){
        pcl::PointXYZRGBNormal  referencePoint = cloud->points[j];
        if(isnan(referencePoint.z) || referencePoint.z > maxRange)
          continue;

        validPoints++;
        referencePoint.getVector3fMap() = poses[k].inverse() * pose * cloud->points[j].getVector3fMap(); ///transform to frame k

        Eigen::Vector2i proj(0,0);
        if(getProjectedPointIndex(referencePoint, proj)){
          float dist = 0.0f;
          int pid = searchProjectedNeighborhood(referencePoint, cloud2, proj, maxDistance, width, height, dist);
          //int pid = proj(1) * width + proj(0);
          if(pid < 0){
            continue;
          }

          ///check if point is valid
          pcl::PointXYZRGBNormal corPointK = cloud2->points[pid];
          if(isValidCorrespondence(cloud->points[j], corPointK, poses[i], poses[k], maxRange, maxAngle, maxSqrDistance)){

//             corPointK.getVector3fMap() = inversePose * poses[k] * corPointK.getVector3fMap(); ///transform into frame i
//             if(getProjectedPointIndex(corPointK, proj)){
//               float otherDist = 0.0f;
//               int bid = searchProjectedNeighborhood(corPointK, cloud2, proj, maxDistance, width, height, otherDist);
//               if(bid < 0 || proj(0) < 2 || proj(0) > width - 2 || proj(1) < 2 || proj(1) > height - 2 || fabs(otherDist - dist) > 0.001){
//                 continue;
//               }
//             }

            assignedPoints++;
          }
        }
      }
      float coverage = (float) assignedPoints / (float) validPoints;
      //std::cerr << k << " " << coverage << " ";
      if(coverage < overlap){
        if(i >= k-2 || k == 1){
          i = k-1;
        } else {
          i = k-2;
        }
        break;
      }
      if(k == cloudNames.size() - 1){
        i = cloudNames.size() - 2;
      }
    }
  }
  return maxCloudSize;
}

void checkPointVertices(SparseSurfaceAdjustmentGraphT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZNormal, ssa::EdgePointXYZNormalPointXYZNormal>& ssaGraph){
  std::cerr << "check point vertices \t  ";
  ssaGraph.dropDataAssociation();
  ssaGraph.dropSensorEdges();
  ///check if vertices had been cleared
  for(size_t i = 0; i < ssaGraph._verticies_points.size(); ++i){
    if(ssaGraph._verticies_points[i] != 0){
      ssa::VertexPointXYZNormal* v =  ssaGraph._verticies_points[i];
      ssaGraph._verticies_points[i] = 0;
      ssaGraph.optimizer()->removeVertex(v);
    }
  }
  ssaGraph._verticies_points.clear();
  std::cerr << "[done]" << std::endl;
}

void alignNonKeyFrames(std::vector<std::string>& cloudNames, std::vector< pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr >& keyClouds, std::vector< Eigen::Affine3f >& poses, std::vector<int>& keyFrames, std::vector<int>& nonKeyFrames, SparseSurfaceAdjustmentGraphT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZNormal, ssa::EdgePointXYZNormalPointXYZNormal>& ssaGraph, std::vector<VertexSE3*>& vertices, float maxDistance, int threads){
  float maxSqrDistance = maxDistance * maxDistance;
  threads = 1;
  std::vector< pcl::search::KdTree<pcl::PointXYZRGBNormal>::Ptr > kdTrees;
  for(size_t i = 0; i < keyClouds.size(); ++i){
    pcl::search::KdTree<pcl::PointXYZRGBNormal>::Ptr kdTree (new pcl::search::KdTree<pcl::PointXYZRGBNormal>);
    kdTree->setInputCloud (keyClouds[i]);
    kdTrees.push_back(kdTree);
  }

  for(size_t i = 0; i < nonKeyFrames.size(); ++i){
    int& id = nonKeyFrames[i];

    std::cerr << "aligning scan " << id << " with keyframe...";
    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
    pcl::io::loadPCDFile(cloudNames[id].c_str() , *cloud);

    int startIndex = 0;
    for(size_t k = 0; k < keyFrames.size() - 1; ++k){
      if(keyFrames[k] < id && id < keyFrames[k+1]){
        startIndex = k;
        break;
      }
    }
    std::cerr << PVAR(startIndex) << " " << keyFrames.size() << std::endl;
    for(int iterations = 0; iterations < 1; ++iterations){
      ///update scan pose with relative transformation to new keyFrame pose
      Eigen::Affine3f pose = vertices[id]->estimate().cast<float>();
      checkPointVertices(ssaGraph);

      std::vector< std::vector< Eigen::Vector4i > > dataAssociation(threads);
      for(int k = startIndex; k < startIndex + 2; ++k){
        //std::cerr << k << " ";
        int keyFrameId = keyFrames[k];
        //std::cerr << keyFrameId << " " << poses[keyFrameId].translation().transpose() << std::endl;

        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& cloud2 = keyClouds[k];
        std::vector<bool> assigned(cloud2->size());
        for(size_t j = 0; j < assigned.size(); ++j)
          assigned[j] = false;
        //#pragma omp parallel for default(shared) schedule(dynamic, 1) num_threads(threads)
        for(size_t j = 0; j < cloud->size(); j+=2){
          pcl::PointXYZRGBNormal  referencePoint = cloud->points[j];
          if(isnan(referencePoint.z) || referencePoint.z > maxRange)
            continue;

          int threadId = omp_get_thread_num();
          std::vector< Eigen::Vector4i >& resultPerThread = dataAssociation[threadId];
          referencePoint.getVector3fMap() = poses[keyFrameId].inverse() * pose * cloud->points[j].getVector3fMap(); ///transform to frame k

          std::vector< int > k_indices;
          std::vector< float > k_sqr_distances;
          kdTrees[k]->nearestKSearch(referencePoint, 1, k_indices, k_sqr_distances);
          //kdTrees[k]->radiusSearch(referencePoint, maxDistance, k_indices, k_sqr_distances);
          if(k_indices.size() > 0 && k_sqr_distances[0] < maxSqrDistance){
            pcl::PointXYZRGBNormal corPointK = cloud2->points[k_indices[0]];
            if(isValidCorrespondence(cloud->points[j], corPointK, poses[id], poses[keyFrameId], maxRange, maxAngle, maxSqrDistance)){
//               corPointK.getVector3fMap() = pose.inverse() * poses[keyFrameId] * corPointK.getVector3fMap();
//               int x = k_indices[0] % cloud->width;
//               int y = (k_indices[0]-x) / cloud->width;
//               std::cerr << " " << k_indices[0] << " " << x << " " << y << " " << y * cloud->width + x << std::endl;
//               Eigen::Vector2i proj(x,y);
//               float otherDist = 0.0f;
              //int bid = searchProjectedNeighborhood(corPointK, cloud2, proj, maxDistance, width, height, otherDist);
              //if((bid == (int) j || fabs(k_sqr_distances[0] - (otherDist * otherDist)) < 0.001) &&
              if(!assigned[k_indices[0]]){
                resultPerThread.push_back(Eigen::Vector4i(id, j, k, k_indices[0]));
                assigned[k_indices[0]] = true;
              }
            }
          }
        }
      }

      int count = 0;
      g2o::OptimizableGraph::EdgeSet eset;
      for(size_t k = 0; k < dataAssociation.size(); ++k){
        for(size_t j = 0; j < dataAssociation[k].size(); ++j){
          Eigen::Vector4i& ref = dataAssociation[k][j];
          const pcl::PointXYZRGBNormal& p1 = cloud->points[ref(1)];
          const pcl::PointXYZRGBNormal& p2 = keyClouds[ref(2)]->points[ref(3)];
          g2o::Edge_V_V_GICP* e = new g2o::Edge_V_V_GICP();
          e->vertices()[0] = dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertices[ref(0)]);
          e->vertices()[1] = dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertices[keyFrames[ref(2)]]);
          //std::cerr << ref(0) << " <-> " << keyFrames[ref(2)] << " " << ref(1) << " <-> " << ref(3) << "     \t     ";
          g2o::EdgeGICP meas;
          meas.pos0 = p1.getVector3fMap().cast<double>();
          meas.normal0 = p1.getNormalVector3fMap().cast<double>();
          meas.pos1 = p2.getVector3fMap().cast<double>();
          meas.normal1 = p2.getNormalVector3fMap().cast<double>();

          e->setMeasurement(meas);
          meas = e->measurement();

          //double uncertainty = 1.0 / ((sensorModel.depthUncertainty(std::max(p1.z, p2.z))) / maxRangeUncertainty);
          /// do plane to plane with additionally using the measurement uncertainty
          //e->information() = meas.prec0(0.01); //.cov0(uncertainty); //
          e->information() = meas.cov0(1000.0); //
          count++;
          ssaGraph.addEdge(e);
          eset.insert(e);
        }
      }
      std::cerr << "assigned " << count  << " points " << std::endl;
      if(count > 0){
        ssaGraph.optimizer()->initializeOptimization(eset);
        ssaGraph.optimizer()->optimize(1);
      }
    }
    vertices[id]->setFixed(true);
    poses[id] = vertices[id]->estimate().cast<float>();
  }
}

void addGICPEdges(SparseSurfaceAdjustmentGraphT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZNormal, ssa::EdgePointXYZNormalPointXYZNormal>& ssaGraph, std::vector< std::vector< Eigen::Vector4i > >& dataAssociation, std::vector<VertexSE3*>& vertices, std::vector< pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr >& keyClouds, std::vector<int>& keyframes){
   //double maxRangeUncertainty = sensorModel.depthUncertainty(maxRange);
  int count = 0;
  for(size_t i = 0; i < dataAssociation.size(); ++i){
    for(size_t j = 0; j < dataAssociation[i].size(); ++j){
      Eigen::Vector4i& ref = dataAssociation[i][j];
      const pcl::PointXYZRGBNormal& p1 = keyClouds[ref(0)]->points[ref(1)];
      const pcl::PointXYZRGBNormal& p2 = keyClouds[ref(2)]->points[ref(3)];
      g2o::Edge_V_V_GICP* e = new g2o::Edge_V_V_GICP();
      e->vertices()[0] = dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertices[keyframes[ref(0)]]);
      e->vertices()[1] = dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertices[keyframes[ref(2)]]);

      g2o::EdgeGICP meas;
      meas.pos0 = p1.getVector3fMap().cast<double>();
      meas.normal0 = p1.getNormalVector3fMap().cast<double>();
      meas.pos1 = p2.getVector3fMap().cast<double>();
      meas.normal1 = p2.getNormalVector3fMap().cast<double>();

      e->setMeasurement(meas);
      meas = e->measurement();

      //double uncertainty = 1.0 / ((sensorModel.depthUncertainty(std::max(p1.z, p2.z))) / maxRangeUncertainty);
      /// do plane to plane with additionally using the measurement uncertainty
      //e->information() = meas.prec0(0.01); //.cov0(uncertainty); //
      e->information() = meas.cov0(1000.0); //
      count++;
      ssaGraph.addEdge(e);
    }
  }
  std::cerr << "assigned " << count  << " points " << std::endl;
}



void solveDAProjectiveKeyFrames(std::vector< pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr >& keyClouds, std::vector< Eigen::Affine3f >& poses, SparseSurfaceAdjustmentGraphT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZNormal, ssa::EdgePointXYZNormalPointXYZNormal>& ssaGraph, std::vector<VertexSE3*>& vertices, std::vector<int>& keyframes, float maxDistance, int threads){

  int search = 1;
  float maxSqrDistance = maxDistance * maxDistance;
  checkPointVertices(ssaGraph);

  std::cerr << "processing da: "<< std::endl;
  std::vector< std::vector< Eigen::Vector4i > > dataAssociation(threads);

  int width = keyClouds[0]->width;
  int height = keyClouds[0]->height;

  //pcl::search::OrganizedNeighbor<PointXYZRGB>
//   std::vector< pcl::search::OrganizedNeighbor<pcl::PointXYZRGBNormal>::Ptr > kdTrees;
//   for(size_t i = 0; i < keyframes.size(); ++i){
//     pcl::search::OrganizedNeighbor<pcl::PointXYZRGBNormal>::Ptr kdTree (new pcl::search::OrganizedNeighbor<pcl::PointXYZRGBNormal>);
//     kdTree->setInputCloud (keyClouds[i]);
//     kdTrees.push_back(kdTree);
//   }

  std::vector< pcl::search::KdTree<pcl::PointXYZRGBNormal>::Ptr > kdTrees;
  for(size_t i = 0; i < keyframes.size(); ++i){
    pcl::search::KdTree<pcl::PointXYZRGBNormal>::Ptr kdTree (new pcl::search::KdTree<pcl::PointXYZRGBNormal>);
    kdTree->setInputCloud (keyClouds[i]);
    kdTrees.push_back(kdTree);
  }

  ///solve data association
  #pragma omp parallel for default(shared) schedule(dynamic, 1) num_threads(threads)
  for(size_t i = 0; i < keyframes.size(); ++i){
    std::cerr <<  i << " ";
    int& id = keyframes[i];
    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& cloud = keyClouds[i];

    int threadId = omp_get_thread_num();
    std::vector< Eigen::Vector4i >& resultPerThread = dataAssociation[threadId];

    Eigen::Affine3f& pose = poses[id];
    Eigen::Affine3f inversePose = pose.inverse();

    for(size_t j = 0; j < cloud->size(); j+=3){
      pcl::PointXYZRGBNormal  referencePoint = cloud->points[j];
      if(isnan(referencePoint.z) || referencePoint.z > maxRange)
        continue;

      size_t maxOtherFrames = std::min(i+2, keyframes.size());
      if(i % 5 == 0 && j % 10 == 0) ///do loop closures only on very 5.th frame for every 2.th point
        maxOtherFrames = keyframes.size();

      for(size_t k = i+1; k < maxOtherFrames; k+=1){
        int& otherId = keyframes[k];
        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& otherCloud = keyClouds[k];
        Eigen::Affine3f& corFramePose = poses[otherId];

        referencePoint.getVector3fMap() = corFramePose.inverse() * pose * cloud->points[j].getVector3fMap(); ///transform to frame k

        if(search == 1){
          std::vector< int > k_indices;
          std::vector< float > k_sqr_distances;
          kdTrees[k]->nearestKSearch(referencePoint, 1, k_indices, k_sqr_distances);
          if(k_sqr_distances.size() > 0 && k_sqr_distances[0] < maxSqrDistance){
            pcl::PointXYZRGBNormal corPointK = otherCloud->points[k_indices[0]];
            if(isValidCorrespondence(cloud->points[j], corPointK, pose, corFramePose, maxRange, maxAngle, maxSqrDistance)){
              resultPerThread.push_back(Eigen::Vector4i(i, j, k, k_indices[0]));
            }
          }
        } else {
          Eigen::Vector2i proj(0,0);
          if(getProjectedPointIndex(referencePoint, proj)){
            float distance = 0.0f;
            int pid = searchProjectedNeighborhood(referencePoint, otherCloud, proj, maxDistance, width, height, distance);
            //int pid = proj(1) * width + proj(0);

            if(pid < 0){
              continue;
            }

            ///check if point is valid
            pcl::PointXYZRGBNormal corPointK = otherCloud->points[pid];
            if(isValidCorrespondence(cloud->points[j], otherCloud->points[pid], pose, corFramePose, maxRange, maxAngle, maxSqrDistance)){

              corPointK.getVector3fMap() = inversePose * corFramePose * corPointK.getVector3fMap(); ///transform into frame i
              if(getProjectedPointIndex(corPointK, proj)){
                float tmpDistance = 0.0f;
                int bid = searchProjectedNeighborhood(corPointK, cloud, proj, maxDistance, width, height, tmpDistance);
                if(bid < 0 || fabs(tmpDistance - distance) > 0.001){
                  continue;
                }
              }

              resultPerThread.push_back(Eigen::Vector4i(i, j, k, pid));
            }
          }
        }
      }
    }
  }
  addGICPEdges(ssaGraph, dataAssociation, vertices, keyClouds, keyframes);
  std::cerr << std::endl << "[done]"<< std::endl;
}




void solveSurfaceDA(pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr& skeleton, SparseSurfaceAdjustmentGraphT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZNormal, ssa::EdgePointXYZNormalPointXYZNormal>& ssaGraph, std::vector<VertexPointXYZNormal*>& pointVertices, float maxDistance){
  ssa::RGBDSensorModel  sensorModel;

  pcl::search::KdTree<pcl::PointXYZRGBNormal>::Ptr kdTree (new pcl::search::KdTree<pcl::PointXYZRGBNormal>);
  kdTree->setInputCloud (skeleton);

  std::tr1::unordered_map< int, std::tr1::unordered_map< int, bool > > assignments;
  ssaGraph.dropDataAssociation();
  int count = 0;
  for(size_t j = 0; j < skeleton->size(); ++j){
    const pcl::PointXYZRGBNormal& point = skeleton->points[j];
//     double maxdistance =  2.0 * sensorModel.depthUncertainty(point.curvature);

    std::vector< int > k_indices;
    std::vector< float > k_sqr_distances;
    kdTree->radiusSearch(point, maxDistance, k_indices, k_sqr_distances);

    for(size_t i = 0; i < k_indices.size(); ++i){
      int& index = k_indices[i];
      if(!assignments[j][index] && (int) j != index){
        if(pointVertices[j] != 0 && pointVertices[index] != 0){
          ssa::EdgePointXYZNormalPointXYZNormal* edge = ssaGraph.createEdge(pointVertices[j], pointVertices[index]);
          assignments[index][j] = true;
          ssaGraph.addEdge(edge);
          count++;
        }
      }
    }
  }
  std::cerr << "added " << count << " surface edges" << std::endl;
}

void usage(){
  std::cerr << "Usage: pcd_ssa_test [options]" << std::endl;
  std::cerr << "-a <int>:           apply global pose optimization for <int> iterations (defaul 20)" << std::endl;
  std::cerr << "-p <string>:        prefix for pointcloud names (default pcloud, starting with pcloud00000.pcd)"  << std::endl;
  std::cerr << "-n <int>:           number of scans (default 1)"  << std::endl;
  std::cerr << "-e <int>:           take only every <int> scan (default 1)" << std::endl;
  std::cerr << "-d <float>:         maximum distance for valid correspondences (default 0.2)" << std::endl;
  std::cerr << "-maxrange <float>:  maximum range for sensor (default 4.0)" << std::endl;
  std::cerr << "-t <int>:           maximum number of threads (default 4)" << std::endl;
}

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

  int numOfScans = 100;
  int threads = 4;
  size_t maxBlockSize = 20;
  int every = 1;
  int gicpIterations = 20;
  double voxelSize = 0.01;
  float maxCorrespondenceDistance = 0.2;
  bool preAlign = false;
  int c=1;
  while (c<argc){
    if (!strcmp(argv[c],"-a")){
      c++;
      preAlign = true;
      gicpIterations=atoi(argv[c]);
      c++;
    } else
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
    if (!strcmp(argv[c],"-d")){
      c++;
      maxCorrespondenceDistance=atof(argv[c]);
      c++;
    } else
    if (!strcmp(argv[c],"-maxrange")){
      c++;
      maxRange=atof(argv[c]);
      c++;
    } else
    if (!strcmp(argv[c],"-r")){
      c++;
      voxelSize=atof(argv[c]);
      c++;
    } else
    if (!strcmp(argv[c],"-t")){
      c++;
      threads=atoi(argv[c]);
      c++;
    } else
    if (!strcmp(argv[c],"-h")){
      usage();
      c++;
    }
  }
  if(argc == 1)
    usage();

  ///Create ssa 3d graph
  SparseSurfaceAdjustmentGraphT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZNormal, ssa::EdgePointXYZNormalPointXYZNormal> ssaGraph;
  g2o::BlockSolverX::LinearSolverType* linearSolver = AllocateLinearSolver<g2o::BlockSolverX>(1); //CHOLMOD
  g2o::BlockSolver< BlockSolverTraits<-1, -1> >* blockSolver = new g2o::BlockSolver< BlockSolverTraits<-1, -1> >(linearSolver);
  g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(blockSolver);
  ssaGraph.optimizer()->setVerbose(true);
  ssaGraph.optimizer()->setAlgorithm(solver);

  PCLSSAHierarchicalT<pcl::PointCloud<pcl::PointXYZRGBNormal> >  pclToSSA;
  std::vector<std::string> cloudNames;
  std::vector< Eigen::Affine3f > poses;
  std::vector<VertexSE3*> vertices;
  std::vector<int> keyframes;

  pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr skeleton(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  std::vector< pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr > keyClouds;

  float resolution = voxelSize;
  for(int i = 0; i < numOfScans; i=i+every)
  {
    char cloudFileName[2048];
    sprintf(cloudFileName, "%s%05d.pcd", pcPrefix.c_str(), i);
    cloudNames.push_back(cloudFileName);
  }

  ///extract keyframes for alignment
  size_t maxCloudSize = getKeyFrames(cloudNames, keyClouds, poses, keyframes, 0.05, 0.75);

  std::vector<int> nonKeyframes;
  for(size_t i = 0; i < keyframes.size() - 1; ++i){
    for(int j = keyframes[i]+1; j < keyframes[i+1]; ++j){
      nonKeyframes.push_back(j);
    }
  }

  ///create full posegraph to "move" also non keyframes
  VertexSE3* lastVertex = 0;
  for(size_t i = 0; i < poses.size(); ++i){
    ///create vertex for graph
    Eigen::Affine3f& pose = poses[i];
    VertexSE3* vertex = pclToSSA.createVertexSE3(pose);
    if(i == 0)
      vertex->setFixed(true);
    ssaGraph.addVertex(vertex);
    vertices.push_back(vertex);

    if(lastVertex){
      EdgeSE3* edge = new EdgeSE3;
      edge->vertices()[0] = lastVertex;
      edge->vertices()[1] = vertex;
      edge->setMeasurement(lastVertex->estimate().inverse() * vertex->estimate());
      edge->information().setIdentity();
      edge->information() *= 0.1;
      ssaGraph.addEdge(edge);
    }
    lastVertex = vertex;
  }

  computeKeyFrameSkeleton(keyClouds, poses, skeleton, keyframes, resolution);
  pcl::io::savePCDFileBinary ("downsampled-input.pcd",  *skeleton);

  ///phase 1 pre aligning
  float stepChange = 0.05 * maxCorrespondenceDistance;
  if(preAlign)
  for(int iterations = 0; iterations < gicpIterations; ++iterations){
    std::cerr << "Iteration: " << iterations << std::endl;
    float maxDist = maxCorrespondenceDistance - (iterations * stepChange);
    maxDist = std::max(maxDist, 0.005f);
    solveDAProjectiveKeyFrames(keyClouds, poses, ssaGraph, vertices, keyframes, maxDist, threads);

    ssaGraph.optimizer()->initializeOptimization();
    ssaGraph.optimizer()->optimize(5);

    ///transfer poses back
    for(size_t j = 0; j < poses.size(); ++j){
      //std::cerr << j << " " << poses[j].translation().transpose() << " <-> ";
      poses[j] = vertices[j]->estimate().cast<float>();
      //std::cerr << poses[j].translation().transpose() << std::endl;
    }

    computeKeyFrameSkeleton(keyClouds, poses, skeleton, keyframes, resolution);
    stringstream ss;
    ss << "skeleton" << iterations << ".pcd";
    pcl::io::savePCDFileBinary (ss.str().c_str(),  *skeleton);

  }

  ///fix pose vertices
  for(size_t i = 0; i < keyframes.size(); ++i){
    vertices[keyframes[i]]->setFixed(true);
  }

  checkPointVertices(ssaGraph);
  ssaGraph.optimizer()->initializeOptimization();
  ssaGraph.optimizer()->optimize(5);

  if(preAlign){
    ///align non keyframes with model
    //alignNonKeyFrames(cloudNames, keyClouds, poses, keyframes, nonKeyframes, ssaGraph, vertices, maxCorrespondenceDistance, threads);
  }

  ofstream poseFile("poses.txt");
  for(size_t i = 0; i < poses.size(); ++i){
    Eigen::Vector3f t = poses[i].translation();
    Eigen::Quaternionf q = (Eigen::Quaternionf) poses[i].linear();
    poseFile << t(0) << " " << t(1) << " " << t(2) << " " << q.w() << " " << q.x() << " "  << q.y() << " " << q.z() << std::endl;
  }
  poseFile.close();

  keyClouds.clear();

  pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr result(new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  ///full data association
  DataAssociation  dataAssociation;
  std::vector<size_t> validIndices;
  std::vector<pcl::PointXYZRGBNormal> refPoints;

  std::vector<bool> pointAssigned(maxCloudSize * cloudNames.size());
  for(size_t t = 0; t < pointAssigned.size(); ++t)
    pointAssigned[t] = false;

  size_t blockCount = std::ceil(cloudNames.size() / maxBlockSize);
  for(size_t x = 0; x < blockCount; ++x){
    solveDAOutOfCore(cloudNames, poses, validIndices, dataAssociation, refPoints, pointAssigned, x, maxBlockSize, maxCloudSize, threads);
    ///phase 2
    for(int iterations = 0; iterations < 1000000; ++iterations){
      buildGraphForNextNEdges(dataAssociation, validIndices, refPoints, poses, ssaGraph, vertices, maxCloudSize, 5000000);

      if(ssaGraph._verticies_points.size() == 0)
        break;

      ssaGraph.optimizer()->initializeOptimization();
      ssaGraph.optimizer()->optimize(5);

      for(size_t j = 0; j < ssaGraph._verticies_points.size(); ++j){
        ssa::VertexPointXYZNormal* v = ssaGraph._verticies_points[j];
        if(v == 0)
          continue;
        pcl::PointXYZRGBNormal point;
        point.getVector3fMap() =  Eigen::Vector3f(v->estimate()(0), v->estimate()(1), v->estimate()(2));
        point.r = v->cr;
        point.g = v->cg;
        point.b = v->cb;
        point.getNormalVector3fMap() = Eigen::Vector3f(v->globalNormal()(0), v->globalNormal()(1), v->globalNormal()(2));
        result->points.push_back(point);
      }
      if(result->size() > 0){
        stringstream resultss;
        resultss << "result" << iterations << ".pcd";
        pcl::io::savePCDFileBinary (resultss.str().c_str(),  *result);
      }
    }
  }


  if(result->size() > 0){
    stringstream resultss;
    resultss << "result" << ".pcd";
    pcl::io::savePCDFileBinary (resultss.str().c_str(),  *result);
    std::cerr << "saving resulting model in " << resultss.str() << std::endl;
  }
}


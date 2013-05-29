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

#ifndef __SSA_PARAMS__
#define __SSA_PARAMS__

#include <iostream>
#include <fstream>

namespace ssa {

  struct BasicParams{
    int version;

    inline BasicParams(){
      version = 0;
    }

    inline void printParams(){
      std::cerr << PVAR(version) << std::endl;
    }
  };

  struct CorrespondenceRejectionParams {
    double maxAngleDifference; /** maximum angle between normals for a valid correspondence in degree */
    double maxColorChannelDiff; /** percentage of allowed color change for a valid correspondence*/
    double maxPoseDistance; /** maximum distance between sensor poses to be considered for assignments*/

    inline CorrespondenceRejectionParams(){
      maxAngleDifference = 10;
      maxColorChannelDiff = 1.0;
      maxPoseDistance = 1.0;
    }

    inline void printParams(){
      std::cerr << PVAR(maxAngleDifference) << std::endl;
      std::cerr << PVAR(maxColorChannelDiff) << std::endl;
      std::cerr << PVAR(maxPoseDistance) << std::endl;
    }
  };

  struct NearestNeighborParams : public CorrespondenceRejectionParams{
    double maxSearchDistance; /** in meter */
    int    numOfNeighbors; /** number of neighbors that should get connected*/
    int    increment; /** search NN for every n-th point */
    bool   onlyIncremental; /** search for NN only for incremental scans */
    int    maxCorrespondencesPerPoint; /** select only the n best correspondences */

    inline NearestNeighborParams(){
      maxSearchDistance = 0.001;
      numOfNeighbors = 1;
      increment = 1;
      onlyIncremental = false;
      maxAngleDifference = 10;
      maxColorChannelDiff = 1.0;
      maxCorrespondencesPerPoint = 5;
    }

    inline void printParams(){
      std::cerr << PVAR(maxSearchDistance) << std::endl;
      std::cerr << PVAR(numOfNeighbors) << std::endl;
      std::cerr << PVAR(increment) << std::endl;
      std::cerr << PVAR(onlyIncremental) << std::endl;
      std::cerr << PVAR(maxAngleDifference) << std::endl;
      std::cerr << PVAR(maxColorChannelDiff) << std::endl;
      std::cerr << PVAR(maxCorrespondencesPerPoint) << std::endl;
    }

  };


  struct NormalShootingParams : public CorrespondenceRejectionParams{
    double maxSearchDistance; /** in meter */
    double stepSize; /** size of single ns step */
    int    steps;
    int    increment;
    bool   onlyIncremental; /** search for NN only for incremental scans */

    inline NormalShootingParams(){
      maxSearchDistance = 0.005;
      stepSize = 0.001;
      steps = 5;
      increment = 1;
      onlyIncremental = false;
      maxAngleDifference = 10;
      maxColorChannelDiff = 1.0;
    }

    inline void printParams(){
      std::cerr << PVAR(maxSearchDistance) << std::endl;
      std::cerr << PVAR(stepSize) << std::endl;
      std::cerr << PVAR(steps) << std::endl;
      std::cerr << PVAR(increment) << std::endl;
      std::cerr << PVAR(onlyIncremental) << std::endl;
      std::cerr << PVAR(maxAngleDifference) << std::endl;
      std::cerr << PVAR(maxColorChannelDiff) << std::endl;
    }
  };

   struct SparseSurfaceAdjustmentParams : public BasicParams{

    enum Sensor{KINECT, LMS};

    int     g2oIterations;
    int     ssaIterations;

    double  normalExtractionMaxNeighborDistance;
    int     normalExtractionMaxNeighbors;
    int     normalExtractionMinNeighbors;

    Sensor  sensor;

    int     outlierRejectionMinConnectedNeighbors;
    bool    optimizeColors;
    int     maxThreads;
    double  targetResolution;

    /** Normal shooting stuff */
    NormalShootingParams    normalShooting;

    /** Nearest neighbot stuff */
    NearestNeighborParams   nearestNeighbor;


    SparseSurfaceAdjustmentParams(){
      g2oIterations = 6;
      ssaIterations = 10;
      normalExtractionMaxNeighborDistance = 0.01;
      normalExtractionMaxNeighbors = 32;
      normalExtractionMinNeighbors = 16;
      sensor = KINECT;
      outlierRejectionMinConnectedNeighbors = 1;
      optimizeColors = false;
      maxThreads = 1;
      targetResolution = 0.002;
    }

    inline void set2DDefaultParams(){
      g2oIterations = 5;
      ssaIterations = 5;
      normalExtractionMaxNeighborDistance = 0.5;
      normalExtractionMaxNeighbors = 8;
      normalExtractionMinNeighbors = 3;
      sensor = LMS;
      outlierRejectionMinConnectedNeighbors = 0;
      optimizeColors = false;
      maxThreads = 4;
      targetResolution = 0.01;
      normalShooting.maxSearchDistance = 0.0;
      normalShooting.stepSize = 0.01;
      normalShooting.steps = 10;
      normalShooting.increment = 1;
      normalShooting.onlyIncremental = 0;
      normalShooting.maxAngleDifference = 60;
      normalShooting.maxColorChannelDiff = 1.0;
      nearestNeighbor.maxSearchDistance = 0.1;
      nearestNeighbor.numOfNeighbors = 1;
      nearestNeighbor.increment = 1;
      nearestNeighbor.onlyIncremental = 0;
      nearestNeighbor.maxAngleDifference = 20;
      nearestNeighbor.maxColorChannelDiff = 1.0;
      nearestNeighbor.maxCorrespondencesPerPoint = 4;
    }


    inline void printParams(){
      std::cerr << PVAR(version) << std::endl;
      std::cerr << PVAR(g2oIterations) << std::endl;
      std::cerr << PVAR(ssaIterations) << std::endl;
      std::cerr << PVAR(normalExtractionMaxNeighborDistance) << std::endl;
      std::cerr << PVAR(normalExtractionMaxNeighbors) << std::endl;
      std::cerr << PVAR(normalExtractionMinNeighbors) << std::endl;
      std::cerr << PVAR(sensor) << std::endl;
      std::cerr << PVAR(outlierRejectionMinConnectedNeighbors) << std::endl;
      std::cerr << PVAR(optimizeColors) << std::endl;
      std::cerr << PVAR(maxThreads) << std::endl;
      std::cerr << PVAR(targetResolution) << std::endl;

      normalShooting.printParams();
      nearestNeighbor.printParams();
    }

    inline void readParams(std::istream& in){
      std::string field;
      while (in.good()){
        in >> field;
        if(field == "version")
          in >> version;
        if(field == "g2oIterations")
          in >> g2oIterations;
        if(field == "ssaIterations")
          in >> ssaIterations;
        if(field == "normalExtractionMaxNeighborDistance")
          in >> normalExtractionMaxNeighborDistance;
        if(field == "normalExtractionMaxNeighbors")
          in >> normalExtractionMaxNeighbors;
        if(field == "normalExtractionMinNeighbors")
          in >> normalExtractionMinNeighbors;
        if(field == "sensor"){
          int i;
          in >> i;
          sensor = (Sensor) i;
        }
        if(field == "outlierRejectionMinConnectedNeighbors")
          in >> outlierRejectionMinConnectedNeighbors;
        if(field == "optimizeColors")
          in >> optimizeColors;
        if(field == "maxThreads")
          in >> maxThreads;
        if(field == "targetResolution")
          in >> targetResolution;

        /** NormalShootingParams */
        if(field == "normalShooting.maxSearchDistance")
          in >> normalShooting.maxSearchDistance;
        if(field == "normalShooting.stepSize")
          in >> normalShooting.stepSize;
        if(field == "normalShooting.steps")
          in >> normalShooting.steps;
        if(field == "normalShooting.increment")
          in >> normalShooting.increment;
        if(field == "normalShooting.onlyIncremental")
          in >> normalShooting.onlyIncremental;
        if(field == "normalShooting.maxAngleDifference")
          in >> normalShooting.maxAngleDifference;
        if(field == "normalShooting.maxColorChannelDiff")
          in >> normalShooting.maxColorChannelDiff;
        if(field == "normalShooting.maxPoseDistance")
          in >> normalShooting.maxPoseDistance;

        /** NearestNeighborsParams */
        if(field == "nearestNeighbor.maxSearchDistance")
          in >> nearestNeighbor.maxSearchDistance;
        if(field == "nearestNeighbor.numOfNeighbors")
          in >> nearestNeighbor.numOfNeighbors;
        if(field == "nearestNeighbor.increment")
          in >> nearestNeighbor.increment;
        if(field == "nearestNeighbor.onlyIncremental")
          in >> nearestNeighbor.onlyIncremental;
        if(field == "nearestNeighbor.maxAngleDifference")
          in >> nearestNeighbor.maxAngleDifference;
        if(field == "nearestNeighbor.maxColorChannelDiff")
          in >> nearestNeighbor.maxColorChannelDiff;
        if(field == "nearestNeighbor.maxPoseDistance")
          in >> nearestNeighbor.maxPoseDistance;
	if(field == "nearestNeighbor.maxCorrespondencesPerPoint")
          in >> nearestNeighbor.maxCorrespondencesPerPoint;
      }
    }

    inline void writeParams(std::ostream& out){
      out << "version " << version << std::endl;
      out << "g2oIterations " << g2oIterations << std::endl;
      out << "ssaIterations " << ssaIterations << std::endl;
      out << "normalExtractionMaxNeighborDistance " << normalExtractionMaxNeighborDistance << std::endl;
      out << "normalExtractionMaxNeighbors " << normalExtractionMaxNeighbors << std::endl;
      out << "normalExtractionMinNeighbors " << normalExtractionMinNeighbors << std::endl;
      out << "sensor " << (int) sensor << std::endl;
      out << "outlierRejectionMinConnectedNeighbors " << outlierRejectionMinConnectedNeighbors << std::endl;
      out << "optimizeColors " << optimizeColors << std::endl;
      out << "maxThreads " << maxThreads << std::endl;
      out << "targetResolution " << targetResolution << std::endl;

      /** NormalShootingParams */
      out << "normalShooting.maxSearchDistance " << normalShooting.maxSearchDistance << std::endl;
      out << "normalShooting.stepSize " << normalShooting.stepSize << std::endl;
      out << "normalShooting.steps " << normalShooting.steps << std::endl;
      out << "normalShooting.onlyIncremental " << normalShooting.onlyIncremental << std::endl;
      out << "normalShooting.increment " << normalShooting.increment << std::endl;
      out << "normalShooting.maxAngleDifference " << normalShooting.maxAngleDifference << std::endl;
      out << "normalShooting.maxColorChannelDiff " << normalShooting.maxColorChannelDiff << std::endl;
      out << "normalShooting.maxPoseDistance " << normalShooting.maxPoseDistance << std::endl;

      out << "nearestNeighbor.maxSearchDistance " << nearestNeighbor.maxSearchDistance << std::endl;
      out << "nearestNeighbor.numOfNeighbors " << nearestNeighbor.numOfNeighbors << std::endl;
      out << "nearestNeighbor.increment " << nearestNeighbor.increment << std::endl;
      out << "nearestNeighbor.onlyIncremental " << nearestNeighbor.onlyIncremental << std::endl;
      out << "nearestNeighbor.maxAngleDifference " << nearestNeighbor.maxAngleDifference << std::endl;
      out << "nearestNeighbor.maxColorChannelDiff " << nearestNeighbor.maxColorChannelDiff << std::endl;
      out << "nearestNeighbor.maxPoseDistance " << nearestNeighbor.maxPoseDistance << std::endl;
      out << "nearestNeighbor.maxCorrespondencesPerPoint " << nearestNeighbor.maxCorrespondencesPerPoint << std::endl;
    }

   };

} //end namespace

#endif
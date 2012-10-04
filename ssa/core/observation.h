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

#ifndef __SPARSE_SURFACE_ADJUSTMENT_OBSERVATION__
#define __SPARSE_SURFACE_ADJUSTMENT_OBSERVATION__

#include <vector>
#include <deque>
#include <Eigen/Core>

#include "ssa/kdtree_flann/kdtree_flann.h"
#include "parameter.h"

namespace ssa {

  template <typename VertexType>
  class Observation : public std::vector<VertexType* >
  {

    static const int Dimension = VertexType::Dimension;
    typedef Eigen::Matrix<double, Dimension, 1> PointVector;
    typedef Eigen::Matrix<double, Dimension, Dimension> PointMatrix;

    typedef KDTreeFlannT<VertexType>     PointTree;

    public:
    inline Observation(){
      computedNeighbors = false;
    }
    void removeVertex(VertexType*& v);

    /** calculate mean */
    static PointVector getMean(std::vector<VertexType* >& neighbors);
    /** calculate covariance */
    static PointMatrix getCovariance(std::vector<VertexType* >& neighbors, PointVector& mean);
    
    /** calculate mean faster */
    static void getMean(std::vector< VertexType* >& neighbors, size_t size, PointVector& mean);
    /** calculate cov faster */
    static void getCovariance(std::vector<VertexType* >& neighbors, size_t size, PointVector& mean, PointMatrix& cov);    

    /** calculates mean and covariance threaded*/
    static void calcMeanCovThreaded(std::vector<VertexType* >& observation, SparseSurfaceAdjustmentParams& params);
    /** calculates mean and covariance threaded and caches the neighborhood information */    
    void calcMeanCovThreadedAndCached(std::vector<VertexType* >& observation, SparseSurfaceAdjustmentParams& params);

    //!returns a kdtree with all observation vertices of this observation
    PointTree getKDTree();

    void calcMeanCov(double& distance);
    void calcMeanCovThreaded(SparseSurfaceAdjustmentParams& params);
    void calcMeanCovThreadedAndCached(SparseSurfaceAdjustmentParams& params);
    void clearCache();
  private:
    void allocNeighborCache(size_t& pointCount, int& maxNeighbors);
    std::vector< std::vector<VertexType*> > neighbors_;
    std::vector< size_t > neighbors_count_;
    //std::map<int, std::deque< VertexType* > > neighbors_;
    bool computedNeighbors;
  };

} //end namespace

#include "observation.hpp"

#endif
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

#ifndef __SSA_LOCAL_NEIGHBORHOOD__
#define __SSA_LOCAL_NEIGHBORHOOD__

#include <vector>
#include <limits>
#include <Eigen/Core>

#include "parameter.h"

namespace ssa {
  /** \brief  this class stores the local neighborhood information for all available surface vertices
  */
  template <typename VertexType>
  class LocalNeighborhoodT
  {

    static const int Dimension = VertexType::Dimension;
    typedef Eigen::Matrix<double, Dimension, 1> PointVector;
    typedef Eigen::Matrix<double, Dimension, Dimension> PointMatrix;

    public:
    inline LocalNeighborhoodT(){
    }

    void addVertex(VertexType*& v, std::vector<VertexType* >& neighborhood);

    void removeVertex(VertexType*& v);
  
    /** clear neighbor cache */
    inline void clear(){ vertex_neighborhood_.clear(); }

    /** create neighborhood cache for full scan */
    void createFromScan(std::vector<VertexType* >& scan, SparseSurfaceAdjustmentParams& params);

    /** calculate mean */
    PointVector getMean(int id);

    /** calculate covariance */
    PointMatrix getCovariance(int id, SparseSurfaceAdjustmentParams& params);

    /** calculate covariance  */
    void getCovariance(int id, PointMatrix& cov, SparseSurfaceAdjustmentParams& params);

    /** access to the neighbors of the vertex with vertexId */
    std::vector< VertexType* >& neighbors(int vertexId);

    VertexType* getNearestNeighbor(VertexType* vertex, int neighborhoodId, double& sqrDistance);

  private:
    std::map<int, std::vector<VertexType* > >    vertex_neighborhood_;
  };

  #include "local_neighborhood.hpp"

} //end namespace

#endif
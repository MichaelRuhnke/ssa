// Sparse Surface Optimization 2D
// Copyright (C) 2011 M. Ruhnke, R. Kuemmerle, G. Grisetti, W. Burgard
// 
// SSA2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// SSA2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __SSA_NEAREST_NEIGHBORT__
#define __SSA_NEAREST_NEIGHBORT__

#include <vector>
#include "ssa/kdtree_flann/kdtree_flann.h"
#include "nearest_neighbor_brute_force.h"
#include <tr1/unordered_map> 

namespace ssa {

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  class NearestNeighborKdtreeT : public NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3> {

    typedef EdgeType1 SLAMEdgeType;
    typedef EdgeType2 SensorEdgeType;
    typedef EdgeType3 DataAssociationEdgeType;

    typedef typename SensorEdgeType::VertexXjType   PointVertex;
    static const int Dj = SensorEdgeType::Dimension;
    typedef Eigen::Matrix<double, Dj, 1> PointVector;
    typedef Eigen::Matrix<double, Dj, Dj> PointMatrix;

    typedef KDTreeFlannT<PointVertex>     PointTree;

    public:

    typedef typename NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::CorrespondenceList  CorrespondenceList;
    typedef typename NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::ScanPairVector      ScanPairVector;


    NearestNeighborKdtreeT();
    ~NearestNeighborKdtreeT();

    using NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::apply;

    /** \brief creates data association edges for all available scan pairs  */
    void apply(ScanPairVector scanPairs, CorrespondenceList& resultingCorrespondences, SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level);

    /** \brief creates an edge between reference and correspondence */
    using NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::assign;

  };
  #include "nearest_neighbor_kdtree.hpp"
}
#endif

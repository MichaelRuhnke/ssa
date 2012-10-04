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

#ifndef __SSA_CORRESPONDENCE_REJECTION__
#define __SSA_CORRESPONDENCE_REJECTION__

#include <vector>
#include "ssa/core/ssa_graph.h"
#include "correspondence.h"

namespace ssa {

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  class CorrespondenceRejectionT{

    typedef EdgeType1 SLAMEdgeType;
    typedef EdgeType2 SensorEdgeType;
    typedef EdgeType3 DataAssociationEdgeType;

    typedef typename SensorEdgeType::VertexXiType   PoseVertex;
    static const int Di = SensorEdgeType::Dimension;
    typedef Eigen::Matrix<double, Di, Di>  PoseMatrix;

    typedef typename SensorEdgeType::VertexXjType   PointVertex;
    static const int Dj = SensorEdgeType::Dimension;
    typedef Eigen::Matrix<double, Dj, 1> PointVector;
    typedef Eigen::Matrix<double, Dj, Dj> PointMatrix;

    typedef KDTreeFlannT<PointVertex>     PointTree;

    public:
    CorrespondenceRejectionT();
    ~CorrespondenceRejectionT();

    static bool isValid(PointVertex*& reference, PointVertex*& correspondence, double& maxAngleDifference);
    static bool isValid(PointVertex*& reference, PointVertex*& correspondence, double& maxAngleDifference, double& maxColorChannel);
    static bool isValid(CorrespondenceT<PointVertex, EdgeType3>& correspondence, double& maxAngleDifference, double& maxColorChannel);

  };
}

#include "correspondence_rejection.hpp"

#endif

// Sparse Surface Adjustment 
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

#ifndef __SSA_CORRESPONDENCE__
#define __SSA_CORRESPONDENCE__

#include <vector>
#include "ssa/core/ssa_graph.h"

namespace ssa {

  template <typename PointVertex, typename DataAssociationEdgeType>
  class CorrespondenceT{
    
    public:     
    CorrespondenceT();
    CorrespondenceT(const CorrespondenceT<PointVertex,DataAssociationEdgeType>& other);
    CorrespondenceT(PointVertex* query, PointVertex* correspondence, double distance);

    PointVertex* query_point;
    PointVertex* cor_point;
    double sqrDistance;
    bool changed;
    DataAssociationEdgeType* edge;
  };
  
  
  template <typename PointVertex, typename DataAssociationEdgeType>
  CorrespondenceT<PointVertex, DataAssociationEdgeType>::CorrespondenceT(): query_point(0), cor_point(0), sqrDistance(0.0), changed(true), edge(0)
  { 
  };
  
  template <typename PointVertex, typename DataAssociationEdgeType>
  CorrespondenceT<PointVertex, DataAssociationEdgeType>::CorrespondenceT(PointVertex* query, PointVertex* correspondence, double distance) : query_point(query), cor_point(correspondence), sqrDistance(distance), changed(true), edge(0)
  { 
    if(query && correspondence){
      double dist = (query->estimate() - correspondence->estimate()).squaredNorm();
      if(fabs(dist-distance) > 1e-7)
        std::cerr << "Provided distance does not match true point pair distance! " << fabs(dist-distance) << std::endl;
    }
  };
  

  template <typename PointVertex, typename DataAssociationEdgeType>
  CorrespondenceT<PointVertex, DataAssociationEdgeType>::CorrespondenceT(const CorrespondenceT<PointVertex, DataAssociationEdgeType>& other): query_point(0), cor_point(0), sqrDistance(0.0), changed(true), edge(0)
  { 
    query_point = other.query_point;
    cor_point = other.cor_point;
    sqrDistance = other.sqrDistance;
  };
}



#endif

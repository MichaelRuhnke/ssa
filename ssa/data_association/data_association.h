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

#ifndef __SSA_DATA_ASSOCIATION__
#define __SSA_DATA_ASSOCIATION__

#include <vector>
#include <unordered_map> 
/** graph and param struct */
#include "ssa/core/parameter.h"
#include "ssa/core/ssa_graph.h"
/** correspondence and search strategies */
#include "correspondence.h"
#include "nearest_neighbor_brute_force.h"
#include "nearest_neighbor_kdtree.h"

namespace ssa {

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  class DataAssociationT{

    typedef EdgeType1 SLAMEdgeType;
    typedef EdgeType2 SensorEdgeType;
    typedef EdgeType3 DataAssociationEdgeType;

    typedef typename SensorEdgeType::VertexXjType   PointVertex;


    typedef typename NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::Correspondence     Correspondence;
    typedef typename NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::CorrespondenceList CorrespondenceList;
    typedef typename NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::OverlapMap         OverlapMap;
  
    public:     
    DataAssociationT();

    /** \brief creates data association edges for all available scan pairs  */
    void apply(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level);

    void applyBasedOnOverlap(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level);

    OverlapMap overlap(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level);

    CorrespondenceList correspondences;

    enum  Strategy{BRUTEFORCE, KDTREE, LOCAL, SPARSE_KDTREE_LOCAL};
    /** getter / setter for search strategy, default search strategy is KDTREE */
    inline void setStrategy(Strategy strategy){strategy_ = strategy;}
    inline Strategy& strategy(){return strategy_ ;}

    private:
      /** \brief search strategy for corresponding vertices of different scans, default KDTREE */
      Strategy    strategy_; 
      OverlapMap  overlap_;
  };

  #include "data_association.hpp"
}



#endif

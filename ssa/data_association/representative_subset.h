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

#ifndef __SSA_REPRESENTATIVE_SUBSET__
#define __SSA_REPRESENTATIVE_SUBSET__

#include <vector>
#include <tr1/unordered_map>
#include "ssa/core/ssa_graph.h"


namespace ssa {

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  class RepresentativeSubsetT{
    
    typedef EdgeType2 SensorEdgeType;    
    typedef typename SensorEdgeType::VertexXjType   PointVertex;
    
    public:     
    RepresentativeSubsetT();

    void createSubset(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>* graph, SparseSurfaceAdjustmentParams& params);
    
    std::vector<PointVertex* > subset;
  };
  
  
  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  RepresentativeSubsetT<EdgeType1, EdgeType2, EdgeType3>::RepresentativeSubsetT()
  { 
  };
  
  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void RepresentativeSubsetT<EdgeType1, EdgeType2, EdgeType3>::createSubset(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>* graph, SparseSurfaceAdjustmentParams& params)
  { 
    subset.clear();
    std::tr1::unordered_map<int, std::tr1::unordered_map<int, std::tr1::unordered_map<int, std::pair<double, PointVertex* > > > > hashGrid;
    std::vector<Eigen::Vector3i> validIndices;
    int numOfVertices = (int) graph->_verticies_points.size();
    for(int i = 0; i < numOfVertices; ++i){
      PointVertex* point = graph->_verticies_points[i];
      Vector3i index;
      index(0) = lrint(point->estimate()(0) / params.targetResolution);
      index(1) = lrint(point->estimate()(1) / params.targetResolution);
      index(2) = lrint(point->estimate()(2) / params.targetResolution);
      double range = 0.0;
      for(g2o::OptimizableGraph::EdgeSet::iterator it=point->edges().begin(); it!=point->edges().end(); it++){
	EdgeType2* e = dynamic_cast<EdgeType2* >(*it);
        if(e){
	  range = e->measurement().norm();
	}
      }
      
      if(range < hashGrid[index(0)][index(1)][index(2)].first || (hashGrid[index(0)][index(1)][index(2)].first == 0 && range > 0.0)){
	if(hashGrid[index(0)][index(1)][index(2)].second != 0){
	  hashGrid[index(0)][index(1)][index(2)].second->covariance().setIdentity();
	} else {
	 validIndices.push_back(index);
	}
        hashGrid[index(0)][index(1)][index(2)].second = point;		
      } else {	
	point->covariance().setIdentity();
      }
    }
    
    numOfVertices = (int) validIndices.size();
    for(int i = 0; i < numOfVertices; ++i){
      Vector3i& index = validIndices[i];
      subset.push_back(hashGrid[index(0)][index(1)][index(2)].second);
    }
    std::cerr << "created subset of size " << subset.size() << std::endl;
  };
  

}



#endif

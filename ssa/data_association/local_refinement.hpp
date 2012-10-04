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


namespace ssa{

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  LocalRefinement<EdgeType1, EdgeType2, EdgeType3>::LocalRefinement()
  {
  };

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  LocalRefinement<EdgeType1, EdgeType2, EdgeType3>::~LocalRefinement()
  {
  };

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void 
  LocalRefinement<EdgeType1, EdgeType2, EdgeType3>::apply(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level = 0)
  {
    #pragma omp parallel for shared(graph) num_threads(params.maxThreads)
    for(int j = 0; j <  (int)graph._edges_data_association.size(); ++j)
    {
      EdgeType3*& edge = graph._edges_data_association[j];
      if(edge->level() >= level)
        refine(edge, graph.neighborCache().neighbors(edge->vertices()[1]->id()));
    }
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void 
  LocalRefinement<EdgeType1, EdgeType2, EdgeType3>::refine(EdgeType3* correspondence, std::vector< PointVertex* >& neighbors)
  {
    PointVertex* from = static_cast<PointVertex* >(correspondence->vertices()[0]);
    PointVertex* to = static_cast<PointVertex* >(correspondence->vertices()[1]);
    double distance = (from->estimate() - to->estimate()).squaredNorm();
    PointVertex* newCorrespondence = 0;
    ///search in local neighborhood for a closer point
    for(size_t j = 0; j < neighbors.size(); ++j){
      PointVertex* candidate =  neighbors[j];
      double candidate_distance = (from->estimate() - candidate->estimate()).squaredNorm();
      if(distance > candidate_distance && candidate->id() != to->id())
      {
        distance = candidate_distance;
        newCorrespondence = candidate;
      }
    }

    ///Update if neccessary
    if(newCorrespondence)
    {
      g2o::HyperGraph::EdgeSet::iterator it = find(to->edges().begin(), to->edges().end(), correspondence);
      if(it != to->edges().end()){
        #pragma omp critical
          to->edges().erase(it);
        correspondence->vertices()[1] = newCorrespondence;
        #pragma omp critical
          newCorrespondence->edges().insert(correspondence);
      } else {
        std::cerr << "Could not find correspondence edge in edge list of target vertex (to)! Graph might be corrupted!" << std::endl;
      }
    }
  }
}
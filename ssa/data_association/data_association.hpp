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


template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
DataAssociationT<EdgeType1, EdgeType2, EdgeType3>::DataAssociationT(): strategy_(DataAssociationT<EdgeType1, EdgeType2, EdgeType3>::KDTREE)
{ 

};


template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
void 
DataAssociationT<EdgeType1, EdgeType2, EdgeType3>::apply(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level = 0)
{
  switch (strategy_) {
    case BRUTEFORCE:
      {
        NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3> daBruteForce;
        daBruteForce.apply(correspondences, graph, params, level);
        break;
      }
    case LOCAL:
      std::cerr << std::endl  << "LOCAL: not yet implemented" << std::endl;
      break;
    case SPARSE_KDTREE_LOCAL:
      {
        std::cerr << std::endl  << "SPARSE_KDTREE_LOCAL: not yet implemented" << std::endl;
//         NearestNeighborSparseKdtreeT<EdgeType1, EdgeType2, EdgeType3> daSparseKdTree;
//         daSparseKdTree.apply(correspondences, graph, params, graph.getMaxLevel());
        break;
      }
    default:
    {
      NearestNeighborKdtreeT<EdgeType1, EdgeType2, EdgeType3> daKdTree;
      daKdTree.apply(correspondences, graph, params, level);
    }
  }

  double timing = get_time();
  std::tr1::unordered_map<int, std::tr1::unordered_map<int, std::tr1::unordered_map<int,  Correspondence* > > > assignments;
  for(size_t i = 0; i < correspondences.size(); ++i)
  {
    Correspondence* correspondence = correspondences[i];
    if(correspondence->query_point && correspondence->cor_point)
      assignments[correspondence->query_point->parentVertex()->id()][correspondence->cor_point->parentVertex()->id()][correspondence->query_point->id()] = correspondence;
  }
  std::cerr << "assignment map computation took " << (get_time() - timing) << "ms. " << std::endl;

  #ifdef __SPARSE_SURFACE_ADJUSTMENT_GRAPH_3D__
  if(level < 1){
  #endif
    timing = get_time();
    //#pragma omp parallel for shared(graph, scanA, correspondences, params) firstprivate(k_indices, k_squared_distances) num_threads(2)   
    for(int i = 0; i < (int) correspondences.size();  ++i)
    {
      Correspondence* correspondence = correspondences[i];
      if(correspondence->query_point && correspondence->cor_point)
      {
        Correspondence* back_correspondence = assignments[correspondence->cor_point->parentVertex()->id()][correspondence->query_point->parentVertex()->id()][correspondence->cor_point->id()];
        if(
          back_correspondence &&
          correspondence->query_point->id() == back_correspondence->cor_point->id() && 
          correspondence->cor_point->id() == back_correspondence->query_point->id() && 
          CorrespondenceRejectionT<EdgeType1, EdgeType2, EdgeType3>::isValid(*correspondences[i], params.nearestNeighbor.maxAngleDifference, params.nearestNeighbor.maxColorChannelDiff)
          )
        {
          EdgeType3* edge = graph.createEdge(correspondences[i]->query_point, correspondences[i]->cor_point);
          if(edge != 0)
          {
            //#pragma omp critical      
            {
              graph.addEdge(edge);
            }
          }
        }
      }
    }
    std::cerr << "edge creation took " << (get_time() - timing) << "ms. " << std::endl;
    
    #ifdef __SPARSE_SURFACE_ADJUSTMENT_GRAPH_3D__
  } else {
//      std::cerr << "linking level " << level << " with gicp edges..." << std::endl;
    for(int i = 0; i < (int) correspondences.size();  ++i)
    {
      Correspondence* correspondence = correspondences[i];
      if(correspondence->query_point && correspondence->cor_point)
      {
        Correspondence* back_correspondence = assignments[correspondence->cor_point->parentVertex()->id()][correspondence->query_point->parentVertex()->id()][correspondence->cor_point->id()];
        if(
          correspondence->query_point->id() == back_correspondence->cor_point->id() && 
          correspondence->cor_point->id() == back_correspondence->query_point->id() && 
          CorrespondenceRejectionT<EdgeType1, EdgeType2, EdgeType3>::isValid(*correspondences[i], params.nearestNeighbor.maxAngleDifference, params.nearestNeighbor.maxColorChannelDiff)
          )
        {

          g2o::Edge_V_V_GICP * e           // new edge with correct cohort for caching
              = new g2o::Edge_V_V_GICP(); 
      
          e->vertices()[0]            // first viewpoint
            = dynamic_cast<g2o::OptimizableGraph::Vertex*>(correspondences[i]->query_point->parentVertex());
      
          e->vertices()[1]            // second viewpoint
            = dynamic_cast<g2o::OptimizableGraph::Vertex*>(correspondences[i]->cor_point->parentVertex());
      
          g2o::EdgeGICP meas;
          meas.pos0 = correspondences[i]->query_point->parentVertex()->estimate().inverse() * correspondences[i]->query_point->estimate();
          meas.normal0 = correspondences[i]->query_point->normal();
          meas.pos1 =  correspondences[i]->cor_point->parentVertex()->estimate().inverse() * correspondences[i]->cor_point->estimate();
          meas.normal1 = correspondences[i]->cor_point->normal();

          e->setMeasurement(meas);
          meas = e->measurement();

          // use this for point-plane
          e->information() = meas.prec0(0.01);
      
          // use this for point-point 
          //e->information().setIdentity();

          graph.addEdge(e);
        }
      }
    }    
  }
  #endif
}

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
void 
DataAssociationT<EdgeType1, EdgeType2, EdgeType3>::applyBasedOnOverlap(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level = 0)
{
  int maxLevel = graph.getMaxLevel();
  OverlapMap overlapMap = overlap(graph, params, maxLevel);

  switch (strategy_){
    case BRUTEFORCE:
      {
        NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3> bruteForce;
        bruteForce.apply(overlapMap, correspondences, graph, params, level);
        break;
      }
    case LOCAL:
      std::cerr << std::endl  << "LOCAL: not yet implemented" << std::endl;
      break;
    case SPARSE_KDTREE_LOCAL:
      {
        std::cerr << std::endl  << "SPARSE_KDTREE_LOCAL: not yet implemented" << std::endl;
        break;
      }
    default:
    {
      NearestNeighborKdtreeT<EdgeType1, EdgeType2, EdgeType3> daKdTree;
      daKdTree.apply(overlapMap, correspondences, graph, params, level);
    }
  }

  std::tr1::unordered_map<int, std::tr1::unordered_map<int, std::tr1::unordered_map<int,  Correspondence* > > > assignments;
  for(size_t i = 0; i < correspondences.size(); ++i)
  {
    Correspondence* correspondence = correspondences[i];
    if(correspondence->query_point && correspondence->cor_point)
      assignments[correspondence->query_point->parentVertex()->id()][correspondence->cor_point->parentVertex()->id()][correspondence->query_point->id()] = correspondence;
  }

  for(int i = 0; i < (int) correspondences.size();  ++i)
  {
    Correspondence* correspondence = correspondences[i];
    if(correspondence->query_point && correspondence->cor_point)
    {
      Correspondence* back_correspondence = assignments[correspondence->cor_point->parentVertex()->id()][correspondence->query_point->parentVertex()->id()][correspondence->cor_point->id()];
      if(
        correspondence->query_point->id() == back_correspondence->cor_point->id() && 
        correspondence->cor_point->id() == back_correspondence->query_point->id() && 
        CorrespondenceRejectionT<EdgeType1, EdgeType2, EdgeType3>::isValid(*correspondences[i], params.nearestNeighbor.maxAngleDifference, params.nearestNeighbor.maxColorChannelDiff)
        )
      {
        EdgeType3* edge = graph.createEdge(correspondences[i]->query_point, correspondences[i]->cor_point);
        if(edge != 0)
        {
	  graph.addEdge(edge);
        }
      }
    }
  }

}

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
typename DataAssociationT<EdgeType1, EdgeType2, EdgeType3>::OverlapMap
DataAssociationT<EdgeType1, EdgeType2, EdgeType3>::overlap(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level = 0)
{
  correspondences.clear();
  NearestNeighborKdtreeT<EdgeType1, EdgeType2, EdgeType3> daKdTree;
  daKdTree.apply(correspondences, graph, params, level);

  std::tr1::unordered_map<int, std::tr1::unordered_map<int, int > > assignmentCount;
  OverlapMap overlap;
  ///prebuild getPointVertices cache
  std::vector<int>  indices = graph.getPoseIds();
  if(indices.size() > 0)
    graph.getPointVertices(indices[0], level);

  for(int j = 0; j < (int) indices.size(); ++j)
  {
    int& id1 = indices[j];
    for(int i = 0; i < (int) indices.size(); ++i)
    {
      int& id2 = indices[i];
       assignmentCount[id1][id2] = 0;
    }
  }

  for(size_t i = 0; i < correspondences.size(); ++i)
  {
    Correspondence* correspondence = correspondences[i];
    if(correspondence->query_point && correspondence->cor_point)
      assignmentCount[correspondence->query_point->parentVertex()->id()][correspondence->cor_point->parentVertex()->id()]++;
  }

  for(int j = 0; j < (int) indices.size(); ++j)
  {
    int& id1 = indices[j];
    std::vector<PointVertex* > vertices_id1 = graph.getPointVertices(id1, level);
    for(int i = 0; i < (int) indices.size(); ++i)
    {
      int& id2 = indices[i];
      overlap[id1][id2] = (float) assignmentCount[id1][id2] / (float) vertices_id1.size();
    }
  }
  return overlap;
}
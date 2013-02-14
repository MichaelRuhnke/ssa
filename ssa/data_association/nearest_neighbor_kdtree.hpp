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

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
NearestNeighborKdtreeT<EdgeType1, EdgeType2, EdgeType3>::NearestNeighborKdtreeT()
{
};

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
NearestNeighborKdtreeT<EdgeType1, EdgeType2, EdgeType3>::~NearestNeighborKdtreeT()
{
};

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
void 
NearestNeighborKdtreeT<EdgeType1, EdgeType2, EdgeType3>::apply(ScanPairVector scanPairs, CorrespondenceList& resultingCorrespondences, SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level = 0)
{
  double timing = get_time();
  /** construct correspondence lists for every thread, for parallel filling */
  std::vector< CorrespondenceList > correspondencesPerThread;
  correspondencesPerThread.resize(params.maxThreads);

  double sqrDistance = params.nearestNeighbor.maxSearchDistance * params.nearestNeighbor.maxSearchDistance;

  ///flann kdtree search result structs
  std::vector<int> k_indices(1);
  std::vector<float> k_squared_distances(1);

  /// pre calculate kdTrees
  int numKdTreeThreads = 4; ///more threads makes no sense from a performance point of view ;) 
  if(params.maxThreads < numKdTreeThreads)
    numKdTreeThreads = params.maxThreads;
  
  std::vector<int>  indices = graph.getPoseIds();
  std::tr1::unordered_map<int, PointTree> kdTrees;
  #pragma omp parallel for schedule(dynamic) shared(graph, indices, params, level, kdTrees) num_threads(numKdTreeThreads) 
  for(int j = 0; j < (int) indices.size(); ++j)
  {
    int& id1 = indices[j];
    std::vector<PointVertex* > vertices_id1 = graph.getPointVertices(id1, level);
    #pragma omp critical
      kdTrees[id1].copyData(vertices_id1, true);
    kdTrees[id1].createKDTree();
  }

  std::cerr << "NearestNeighborKdtreeT pre calculation took " << (get_time() - timing) * 1000.0 << "ms with " << numKdTreeThreads << " threads." << std::endl;
  timing = get_time();
  #pragma omp parallel for schedule(dynamic) shared(graph, params, level, sqrDistance, kdTrees) num_threads(params.maxThreads) firstprivate(k_indices, k_squared_distances)
  for(int j = 0; j < (int) scanPairs.size(); ++j)
  {
    int threadId = omp_get_thread_num();

    int& id1 = scanPairs[j].first;
    int& id2 = scanPairs[j].second;
    PointTree& kDTree = kdTrees[id1];

    std::vector<PointVertex* >& vertices_id1 = graph.getPointVertices(id1, level);
    std::vector<PointVertex* >& vertices_id2 = graph.getPointVertices(id2, level);

    for(size_t k = 0; k < vertices_id2.size(); k=k+params.nearestNeighbor.increment)
    {
      /// FLANN nearest k search
      kDTree.nearestKSearch(vertices_id2[k], 1, k_indices, k_squared_distances);

      if(k_squared_distances[0] < sqrDistance){
        assign(correspondencesPerThread[threadId], vertices_id2[k], vertices_id1[k_indices[0]], k_squared_distances[0]);
      } 
    }
  }
//   std::cerr << "NearestNeighborKdtreeT calculation took " << (get_time() - timing) * 1000.0 << "ms." << std::endl;
//   timing = get_time();
  ///merge results into single list
  for(size_t i = 0; i <  correspondencesPerThread.size(); ++i)
  {
    copy (correspondencesPerThread[i].begin(),correspondencesPerThread[i].end(),back_inserter(resultingCorrespondences));
  }
//   std::cerr << "NearestNeighborKdtreeT result merging took " << (get_time() - timing) * 1000.0 << "ms." << std::endl;
}

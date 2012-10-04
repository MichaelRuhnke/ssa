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
NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::NearestNeighborBruteForceT()
{
};

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::~NearestNeighborBruteForceT()
{
};

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
void 
NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::apply(CorrespondenceList& resultingCorrespondences, SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level = 0)
{
  /** get valid scan pairs */
  std::vector< std::pair <int,int> > scanPairs;
  if(params.nearestNeighbor.onlyIncremental == 1){
    scanPairs = getIncrementalScanPairs(graph);
  }  else {
    scanPairs = getScanPairs(graph);
  } 
  std::cerr << PVAR(scanPairs.size()) << std::endl;
  /** pre build PointVertex Observation Cache*/
  if(scanPairs.size() > 0)
    graph.getPointVertices(scanPairs[0].first, level);
  apply(scanPairs, resultingCorrespondences, graph, params, level);
}

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
void 
NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::apply(OverlapMap overlap, CorrespondenceList& resultingCorrespondences, SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level = 0)
{
  /** get valid scan pairs */
  std::vector< std::pair <int,int> > scanPairs = getScanPairs(graph, overlap, 0.6);
  /** pre build PointVertex Observation Cache*/
  if(scanPairs.size() > 0)
    graph.getPointVertices(scanPairs[0].first, level);
  apply(scanPairs, resultingCorrespondences, graph, params, level);
}

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
void 
NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::apply(ScanPairVector scanPairs, CorrespondenceList& resultingCorrespondences, SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level = 0)
{
  /** construct correspondence lists for every thread, for parallel filling */
  std::vector< CorrespondenceList > correspondencesPerThread;
  correspondencesPerThread.resize(params.maxThreads);

  double sqrDistance = params.nearestNeighbor.maxSearchDistance * params.nearestNeighbor.maxSearchDistance;
  /** search for correspondences */
  #pragma omp parallel for schedule(dynamic, 1) shared(graph, params, level, correspondencesPerThread, sqrDistance) num_threads(params.maxThreads)
  for(int j = 0; j <  (int) scanPairs.size(); ++j)
  {
    int threadId = omp_get_thread_num();
    applyForScanPair(scanPairs[j].first, scanPairs[j].second, correspondencesPerThread[threadId], graph, sqrDistance, level);
  } 

  ///merge results into single list
  for(size_t i = 0; i <  correspondencesPerThread.size(); ++i)
  {
    copy (correspondencesPerThread[i].begin(),correspondencesPerThread[i].end(),back_inserter(resultingCorrespondences));
  }
}

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
void 
NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::applyForScanPair(int id1, int id2, CorrespondenceList& resultingCorrespondences, SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, double sqrDistance, int level = 0)
{
  std::vector<PointVertex* >& vertices_id1 = graph.getPointVertices(id1, level);
  std::vector<PointVertex* >& vertices_id2 = graph.getPointVertices(id2, level);
  /** we assign the points from vertices_id2 to their nearest neighbors in vertices_id1 */
  for(int l = 0; l <  (int) vertices_id2.size(); ++l)
  {
    PointVertex*& vertex = vertices_id2[l];
    PointVertex* cor_vertex = 0;
    double minDistance = sqrDistance;
    for(int m = 0; m <  (int) vertices_id1.size(); ++m)
    {
      PointVertex*& candidate = vertices_id1[m];
      double distance = (candidate->estimate() - vertex->estimate()).squaredNorm();
      if(distance < minDistance)
      {
        minDistance = distance;
        cor_vertex = candidate;
      }
    }
    if(cor_vertex){
      assign(resultingCorrespondences, vertex, cor_vertex, minDistance);
    }
  }
}

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
typename NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::ScanPairVector 
NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::getScanPairs(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph)
{
  std::vector< std::pair <int,int> > scanPairs;
  std::vector<int>  indices = graph.getPoseIds();
  for(int j = 0; j <  (int) indices.size(); ++j)
  {
    int& id1 = indices[j];
    for(size_t k = 0; k <  indices.size(); ++k)
    {
      int& id2 = indices[k];
      if(id1 == id2)
        continue;
      std::pair <int,int> scanPair(id1, id2);
      scanPairs.push_back(scanPair);
    }
  }
  return scanPairs;
}

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
typename NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::ScanPairVector 
NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::getIncrementalScanPairs(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph)
{
  std::vector< std::pair <int,int> > scanPairs;
  std::vector<int>  indices = graph.getPoseIds();
  for(int j = 0; j <  (int) (indices.size() - 1); ++j)
  {
    std::pair <int,int> scanPair(indices[j], indices[j+1]);
    scanPairs.push_back(scanPair);
    std::pair <int,int> scanPairBack(indices[j+1], indices[j]);
    scanPairs.push_back(scanPairBack);
  }
  return scanPairs;
}

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
typename NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::ScanPairVector 
NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::getScanPairs(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, OverlapMap overlap, float minOverlap)
{
  std::vector< std::pair <int,int> > scanPairs;
  std::vector<int>  indices = graph.getPoseIds();
  for(int j = 0; j <  (int) indices.size(); ++j)
  {
    int& id1 = indices[j];
    for(size_t k = 0; k <  indices.size(); ++k)
    {
      int& id2 = indices[k];
      if(id1 == id2)
        continue;
      if(overlap[id1][id2] > minOverlap){
        std::pair <int,int> scanPair(id1, id2);
        scanPairs.push_back(scanPair);
      }
    }
  }
  return scanPairs;
}

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
void 
NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::assign(CorrespondenceList& resultingCorrespondences, PointVertex*& ref_vertex, PointVertex*& cor_vertex, float distance)
{
  Correspondence* correspondence = new Correspondence(ref_vertex, cor_vertex, distance); 
  // #pragma omp critical
  resultingCorrespondences.push_back(correspondence);
}

template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
void 
NearestNeighborBruteForceT<EdgeType1, EdgeType2, EdgeType3>::assign(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, PointVertex*& reference, PointVertex*& correspondence)
{
  EdgeType3* edge = graph.createEdge(reference, correspondence);
  if(edge != 0)
  {
    #pragma omp critical      
    {
      graph.addEdge(edge);
    }
  }
}

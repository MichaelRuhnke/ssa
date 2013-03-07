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


template <typename VertexType>
void LocalNeighborhoodT<VertexType>::addVertex(VertexType*& v, std::vector<VertexType* >& neighborhood)
{
  ///copy neighborhood information
  copy (neighborhood.begin(),neighborhood.end(),back_inserter(vertex_neighborhood_[v->id()]));
}

template <typename VertexType>
void LocalNeighborhoodT<VertexType>::removeVertex(VertexType*& v){
  vertex_neighborhood_[v->id()].clear();
}

template <typename VertexType>
void LocalNeighborhoodT<VertexType>::createFromScan(std::vector<VertexType* >& scan, SparseSurfaceAdjustmentParams& params){
  ///create kdtree of local scan
  KDTreeFlannT<VertexType>  kDTree;
  kDTree.copyData(scan, true);
  kDTree.createKDTree();

  ///search for neighbors within distance
  double sqrDistance = params.normalExtractionMaxNeighborDistance * params.normalExtractionMaxNeighborDistance;
  std::vector<int> k_indices;
  std::vector<float> k_squared_distances;

  for(unsigned int j = 0; j < scan.size(); ++j){
    VertexType*& point = scan[j];
    ///FLANN Neighbor search
    kDTree.nearestKSearch(point, params.normalExtractionMaxNeighbors, k_indices, k_squared_distances);
    for(size_t i  = 0; i < k_indices.size(); ++i){
      if(k_squared_distances[i] < sqrDistance){
        #pragma omp critical
          vertex_neighborhood_[point->id()].push_back(scan[k_indices[i]]);
      }
    }
  }
} 

template <typename VertexType>
typename LocalNeighborhoodT<VertexType>::PointVector LocalNeighborhoodT<VertexType>::getMean(int id)
{
  std::vector< VertexType* >& neighbors = vertex_neighborhood_[id];
  PointVector mean = PointVector::Zero();
  for(unsigned int i = 0; i < neighbors.size(); ++i){
      mean += neighbors[i]->estimate();
  }
  mean = mean / (double) neighbors.size();
  return mean;
}
  
template <typename VertexType>
typename LocalNeighborhoodT<VertexType>::PointMatrix LocalNeighborhoodT<VertexType>::getCovariance(int id, SparseSurfaceAdjustmentParams& params){
  PointMatrix cov = PointMatrix::Identity(); 
  getCovariance(id, cov, params);
  return cov;
}
  
template <typename VertexType>
void LocalNeighborhoodT<VertexType>::getCovariance(int id, typename LocalNeighborhoodT<VertexType>::PointMatrix& cov, SparseSurfaceAdjustmentParams& params){

  std::vector< VertexType* >& neighbors = vertex_neighborhood_[id];
  /// can not calculate proper covariance so skip the next steps
  if((int) neighbors.size() < params.normalExtractionMinNeighbors)
  {
    cov = PointMatrix::Identity(); 
    return;
  }
  PointVector mean = getMean(id);
  cov = PointMatrix::Zero(); 

  for(size_t j = 0; j < neighbors.size(); ++j){
    cov +=  (neighbors[j]->estimate() - mean) *  (neighbors[j]->estimate() - mean).transpose();
  }
  double factor = 1.0 / ((double) neighbors.size() - 1.0);
  cov =  factor * cov;

  //add small value to diag
  for(int i = 0; i < Dimension; ++i){
    cov(i,i) += 1e-7;
  }
}

template <typename VertexType>
std::vector< VertexType* >& LocalNeighborhoodT<VertexType>::neighbors(int vertexId){
  return vertex_neighborhood_[vertexId];
}

template <typename VertexType>
VertexType* LocalNeighborhoodT<VertexType>::getNearestNeighbor(VertexType* vertex, int neighborhoodId, double& sqrDistance){
  double maxSqrDistance = 5.0 * sqrDistance;
  std::vector< VertexType* >& neighbors = vertex_neighborhood_[neighborhoodId];
  VertexType* newCorrespondence = 0;
  for(size_t j = 0; j < neighbors.size(); ++j){
    VertexType* candidate =  neighbors[j];
    double candidate_distance = (vertex->estimate() - candidate->estimate()).squaredNorm();
    if(candidate_distance < maxSqrDistance && fabs(candidate_distance - maxSqrDistance) > std::numeric_limits<double>::epsilon())
    {
      maxSqrDistance = candidate_distance;
      newCorrespondence = candidate;
    }
  }

  VertexType* startCorrespondence = newCorrespondence;
  int refineSteps = 3;
  for(int i = 0; i < refineSteps; ++i){
    if(newCorrespondence){
      neighbors = vertex_neighborhood_[newCorrespondence->id()];
      for(size_t j = 0; j < neighbors.size(); ++j){
        VertexType* candidate =  neighbors[j];
        double candidate_distance = (vertex->estimate() - candidate->estimate()).squaredNorm();
        if(candidate_distance < maxSqrDistance && fabs(candidate_distance - maxSqrDistance) > std::numeric_limits<double>::epsilon())
        {
          maxSqrDistance = candidate_distance;
          newCorrespondence = candidate;
        }
      }
      if(newCorrespondence == startCorrespondence)
        break;
    }
  }
  if(maxSqrDistance > sqrDistance)
    newCorrespondence = 0;
  sqrDistance = maxSqrDistance;
  return newCorrespondence;
}
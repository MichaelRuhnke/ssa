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

namespace ssa {

  template <typename VertexType>
  void Observation<VertexType>::removeVertex(VertexType*& v){
    typename Observation<VertexType>::iterator it = std::find(this->begin(), this->end(), v);
    if(it != this->end()){
      this->erase(it);
    } else {
      std::cerr << "ERROR: vertex " << v->id() << " not found in observation."  << std::endl;
      sleep(1);
    }
  }

  template <typename VertexType>
  typename Observation<VertexType>::PointVector Observation<VertexType>::getMean(std::vector< VertexType* >& neighbors){
    PointVector mean = PointVector::Zero();
    for(unsigned int i = 0; i < neighbors.size(); ++i){
        mean += neighbors[i]->estimate();
    }
    mean = mean / (double) neighbors.size();
    return mean;
  }
  
  template <typename VertexType>
  void Observation<VertexType>::getMean(std::vector< VertexType* >& neighbors, size_t size,  PointVector& mean){
    mean = PointVector::Zero();
    for(size_t i = 0; i < size; ++i){
        mean += neighbors[i]->estimate();
    }
    mean = mean / (double) neighbors.size();
  }
  
  template <typename VertexType>
  typename Observation<VertexType>::PointMatrix Observation<VertexType>::getCovariance(std::vector< VertexType* >& neighbors, PointVector& mean){
    PointMatrix cov = PointMatrix::Zero(); 
    for(size_t j = 0; j < neighbors.size(); ++j){
      cov +=  (neighbors[j]->estimate() - mean) *  (neighbors[j]->estimate() - mean).transpose();
    }
    double factor = 1.0 / ((double) neighbors.size() - 1.0);
    cov =  factor * cov;

    //add small value to diag
    for(int i = 0; i < Dimension; ++i){
      cov(i,i) += 1e-7;
    }

    return cov;
  }
  
  template <typename VertexType>
  void Observation<VertexType>::getCovariance(std::vector< VertexType* >& neighbors, size_t size, PointVector& mean, typename Observation<VertexType>::PointMatrix& cov){
    cov = PointMatrix::Zero(); 
    for(size_t j = 0; j < size; ++j){
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
  void Observation<VertexType>::calcMeanCovThreaded(std::vector<VertexType* >& observation, SparseSurfaceAdjustmentParams& params){

    PointTree kDTree;
    kDTree.copyData(observation, true);
    kDTree.createKDTree();

    double sqrDistance = params.normalExtractionMaxNeighborDistance * params.normalExtractionMaxNeighborDistance;
    #pragma omp parallel for schedule(dynamic, 20) shared(observation, sqrDistance, kDTree) 
    for(unsigned int j = 0; j < observation.size(); ++j){
      VertexType* point = observation[j];

      //FLANN Neighbor search
      std::vector<int> k_indices;
      std::vector<float> k_squared_distances;
      kDTree.nearestKSearch(point, params.normalExtractionMaxNeighbors, k_indices, k_squared_distances);

      //prepare correspondences
      std::vector< VertexType* > neighbors;
      neighbors.reserve(k_indices.size());

      //store neighbors
      for(size_t i  = 0; i < k_indices.size(); ++i){
        if(k_squared_distances[i] < sqrDistance){
          VertexType* neighbor = observation[k_indices[i]];
          neighbors.push_back(neighbor);
        }
      }

      if((int) neighbors.size() >= params.normalExtractionMinNeighbors){
        size_t size = neighbors.size();
        //calc mean and covariance
        PointVector mean;
        Observation<VertexType>::getMean(neighbors, size, mean);
        point->covariance() = Observation<VertexType>:: getCovariance(neighbors, size, mean);
      } else {
        point->covariance()= PointMatrix::Identity();
      }
      point->_hasNormal = false;
    }
  }

  template <typename VertexType>
  typename Observation<VertexType>::PointTree Observation<VertexType>::getKDTree(){
    std::vector<VertexType* >& observation = *this;

    PointTree kDTree;
    kDTree.copyData(observation, true);
    kDTree.createKDTree();

    return kDTree;
  }
  
  template <typename VertexType>
  void Observation<VertexType>::calcMeanCov(double& distance){
    int maxNeighbors = 32;
    Observation<VertexType>& observation = *this;
    PointTree kdTree = observation.getKDTree();
    double sqrDistance = distance * distance;
    #pragma omp parallel for schedule(dynamic, 20) shared(observation, sqrDistance, kdTree) 
    for(unsigned int j = 0; j < observation.size(); ++j){
      VertexType* point = observation[j];

      //FLANN Neighbor search
      std::vector<int> k_indices;
      std::vector<float> k_squared_distances;
      kdTree.nearestKSearch(point, maxNeighbors, k_indices, k_squared_distances);

      //prepare correspondences
      std::vector< VertexType* > neighbors;
      neighbors.reserve(k_indices.size());

      for(size_t i  = 0; i < k_indices.size(); ++i){
        if(k_squared_distances[i] < sqrDistance){
          VertexType* correspondence = observation[k_indices[i]];
          neighbors.push_back(correspondence);
        }
      }

      if(neighbors.size() > maxNeighbors * 0.5){
        //calc mean and covariance
        PointVector mean = Observation<VertexType>::getMean(neighbors);
        point->covariance() = Observation<VertexType>::getCovariance(neighbors, mean);
      } else {
        point->covariance()= PointMatrix::Identity();
      }
      point->_hasNormal = false;
    }
  }

  template <typename VertexType>
  void Observation<VertexType>::calcMeanCovThreaded(SparseSurfaceAdjustmentParams& params){
    Observation<VertexType>& observation = *this;

    PointTree kdTree = observation.getKDTree();
    double sqrDistance = params.normalExtractionMaxNeighborDistance * params.normalExtractionMaxNeighborDistance;
    #pragma omp parallel for schedule(dynamic, 10) shared(observation, sqrDistance, kdTree) 
    for(unsigned int j = 0; j < observation.size(); ++j){
      VertexType* point = observation[j];

      //FLANN Neighbor search
      std::vector<int> k_indices;
      std::vector<float> k_squared_distances;
      kdTree.nearestKSearch(point, params.normalExtractionMaxNeighbors, k_indices, k_squared_distances);

      //prepare correspondences
      std::vector< VertexType* > neighbors;
      neighbors.reserve(k_indices.size());

      //store neighbors
      for(size_t i  = 0; i < k_indices.size(); ++i){
        if(k_squared_distances[i] < sqrDistance){
          VertexType* neighbor = observation[k_indices[i]];
          neighbors.push_back(neighbor);
        }
      }

      if((int) neighbors.size() >= params.normalExtractionMinNeighbors){
        //calc mean and covariance
        PointVector mean = Observation<VertexType>::getMean(neighbors);
        point->covariance() = Observation<VertexType>:: getCovariance(neighbors, mean);
      } else {
        point->covariance()= PointMatrix::Identity();
      }
//       point->_hasNormal = false;
    }
  }

  template <typename VertexType>
  void Observation<VertexType>::calcMeanCovThreadedAndCached(SparseSurfaceAdjustmentParams& params){
    Observation<VertexType>& observation = *this;

    if(!computedNeighbors){
      neighbors_.resize(observation.size());
      PointTree kdTree = observation.getKDTree();
      double sqrDistance = params.normalExtractionMaxNeighborDistance * params.normalExtractionMaxNeighborDistance;
      #pragma omp parallel for schedule(dynamic, 10) shared(observation, sqrDistance, kdTree) 
      for(unsigned int j = 0; j < observation.size(); ++j){
        VertexType* point = observation[j];

        //FLANN Neighbor search
        std::vector<int> k_indices;
        std::vector<float> k_squared_distances;
        kdTree.nearestKSearch(point, params.normalExtractionMaxNeighbors, k_indices, k_squared_distances);
  
        //store neighbors
        for(size_t i  = 0; i < k_indices.size(); ++i){
          if(k_squared_distances[i] < sqrDistance){
            VertexType* neighbor = observation[k_indices[i]];
            #pragma omp critical
            neighbors_[j].push_back(neighbor);
          }
        }
      }
      computedNeighbors = true;
    }

    #pragma omp parallel for schedule(dynamic, 10) shared(observation) 
    for(unsigned int j = 0; j < observation.size(); ++j){
      VertexType* point = observation[j];
      std::vector< VertexType* >& neighbors = neighbors_[j];
      if((int) neighbors.size() >= params.normalExtractionMinNeighbors){
        //calc mean and covariance
        PointVector mean = Observation<VertexType>::getMean(neighbors);
        point->covariance() = Observation<VertexType>:: getCovariance(neighbors, mean);
      } else {
        point->covariance()= PointMatrix::Identity();
      }
//       point->_hasNormal = false;
    }
  }
  
  template <typename VertexType>
  void Observation<VertexType>::calcMeanCovThreadedAndCached(std::vector<VertexType* >& observation, SparseSurfaceAdjustmentParams& params){
        
    if(!computedNeighbors){
      std::vector<int> k_indices(params.normalExtractionMaxNeighbors);
      std::vector<float> k_squared_distances(params.normalExtractionMaxNeighbors);    
      ///allocating memory for neighbors caching
      size_t pointCount = observation.size();
      allocNeighborCache(pointCount, params.normalExtractionMaxNeighbors);

      ///create kD-Tree
      PointTree kdTree;      
      kdTree.copyData(observation, true);
      kdTree.createKDTree();

      #pragma omp parallel for schedule(dynamic, 20) shared(observation, kdTree, params) firstprivate(k_indices, k_squared_distances) 
      for(unsigned int j = 0; j < observation.size(); ++j){
        ///FLANN neighbor search
        neighbors_count_[j] = kdTree.radiusSearch(observation[j], params.normalExtractionMaxNeighborDistance, k_indices, k_squared_distances, params.normalExtractionMaxNeighbors);
        //kdTree.nearestKSearch(observation[j], params.normalExtractionMaxNeighbors, k_indices, k_squared_distances);

        ///store neighbor pointer
        for(size_t i  = 0; i < neighbors_count_[j]; ++i){
          neighbors_[j][i] = observation[k_indices[i]];
        }
      }
      computedNeighbors = true;
    }

    double timing = g2o::get_time();
    ///calculate mean and covariance
    PointVector mean;
    VertexType* point = 0;
    #pragma omp parallel for schedule(dynamic, 20) shared(observation) firstprivate(mean, point)
    for(unsigned int j = 0; j < observation.size(); ++j){
      point = observation[j];
      if(neighbors_count_[j] >= (size_t) params.normalExtractionMinNeighbors){
        Observation<VertexType>::getMean(neighbors_[j], neighbors_count_[j], mean);
	PointMatrix covariance;
        Observation<VertexType>::getCovariance(neighbors_[j], neighbors_count_[j], mean, covariance);
	point->updateNormal(covariance);
      } else {
        point->covariance()= PointMatrix::Identity();
      }
//       point->_hasNormal = false;
    }
    std::cerr << "normal calculation took: " << (g2o::get_time() - timing) * 1000 << "ms." << std::endl;
  }

  template <typename VertexType>
  void Observation<VertexType>::clearCache(){
    computedNeighbors = false;
    neighbors_.clear();
  }

    template <typename VertexType>
  void Observation<VertexType>::allocNeighborCache(size_t& pointCount, int& maxNeighbors)
  {
    ///allocating memory for neighbors caching
    neighbors_.resize(1);
    neighbors_[0].resize(maxNeighbors);
    neighbors_.resize(pointCount, neighbors_[0]);
    neighbors_count_.resize(pointCount);
  }
  

  
} //end namespace
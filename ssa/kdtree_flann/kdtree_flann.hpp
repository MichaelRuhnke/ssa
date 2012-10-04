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

  template <typename PointVertexType>
  KDTreeFlannT<PointVertexType>::KDTreeFlannT() : data_(0), dim_(PointVertexType::Dimension), flann_index_(0), mode_(0)
  {
  };

  template <typename PointVertexType>
  KDTreeFlannT<PointVertexType>::KDTreeFlannT(int dim) : data_(0), dim_(dim), flann_index_(0), mode_(0)
  {
  };

  template <typename PointVertexType>
  KDTreeFlannT<PointVertexType>::~KDTreeFlannT()
  {
    this->clear();
  };

  template <typename PointVertexType>
  void KDTreeFlannT<PointVertexType>::copyData(std::vector<PointVertexType* >& vertices, bool useVectorIndices = false){

//     double timing = get_time();
    int num_of_points = vertices.size();

    data_ = (float*)malloc (num_of_points * dim_ * sizeof(float));
    index_mapping_.resize(num_of_points);

//     #pragma omp parallel for schedule(dynamic, 5) shared(num_of_points, useVectorIndices) 
    for (int v_index = 0; v_index < num_of_points; ++v_index)
    {
      const PointVertexType* v = vertices[v_index];
      if(useVectorIndices){
        index_mapping_[v_index] = v_index;
      } else {
        index_mapping_[v_index] = v->id();
      }
      for(int i = 0; i < dim_; ++i)
        data_[v_index * dim_ + i] = (float) v->estimate()[i];
    }
//    std::cerr << "FLANN data transfer took " << get_time() - timing << "ms."<< std::endl;
  }

  template <typename PointVertexType> 
  void KDTreeFlannT<PointVertexType>::createKDTree(){
//     double timing = get_time();
    if(index_mapping_.size() == 0){
      std::cerr << "Warning: you tried to create an empty kDTree." << std::endl;
      return;
    }
    // construct an randomized kd-tree index using 4 kd-trees
    flann_index_ = new flann::Index<flann::L2_Simple<float> >(flann::Matrix<float>(data_, index_mapping_.size(), dim_), flann::KDTreeSingleIndexParams(10));
    //flann::KDTreeSingleIndexParams(15) flann::AutotunedIndexParams(0.99, 0.01, 0, 0.1)
    flann_index_->buildIndex();
//     std::cerr << "FLANN kDTree construction took " << get_time() - timing << "ms."<< std::endl;
  }

  template <typename PointVertexType>
  void KDTreeFlannT<PointVertexType>::nearestKSearch(const PointVertexType* vertex, int k, 
                                       std::vector<int> &k_indices, 
                                       std::vector<float> &k_squared_distances){
    if(index_mapping_.size() == 0){
      std::cerr << "Warning: nearestKSearch on an empty kDTree." << std::endl;
      return;
    }

    std::vector<float> tmp (dim_);
    KDTreeFlannT<PointVertexType>::vectorize(dim_, vertex, tmp);

    k_indices.resize(k);
    k_squared_distances.resize(k);

    flann::Matrix<int> k_indices_mat (&k_indices[0], 1, k);
    flann::Matrix<float> k_distances_mat (&k_squared_distances[0], 1, k);

    // do a knn search
    flann_index_->knnSearch (flann::Matrix<float>(&tmp[0], 1, dim_), k_indices_mat, k_distances_mat, k, flann::SearchParams (-1));

    //change matrix indices to vertex indices
    for (size_t i = 0; i < k_indices.size (); ++i)
    {
      k_indices[i] = index_mapping_[k_indices[i]];
    }
  }

  template <typename PointVertexType>
  void KDTreeFlannT<PointVertexType>::nearestKSearch(const Eigen::Vector3f& point, int k, 
                                       std::vector<int> &k_indices, 
                                       std::vector<float> &k_squared_distances){

    if(index_mapping_.size() == 0){
      std::cerr << "Warning: nearestKSearch on an empty kDTree." << std::endl;
      return;
    }

    if(dim_ != 3)
      dim_ = 3;

    std::vector<float> tmp (dim_);
    KDTreeFlannT<PointVertexType>::vectorize(dim_, point, tmp);

    k_indices.resize(k);
    k_squared_distances.resize(k);

    flann::Matrix<int> k_indices_mat (&k_indices[0], 1, k);
    flann::Matrix<float> k_distances_mat (&k_squared_distances[0], 1, k);

    // do a knn search
    flann_index_->knnSearch (flann::Matrix<float>(&tmp[0], 1, dim_), k_indices_mat, k_distances_mat, k, flann::SearchParams (-1));

    //change matrix indices to vertex indices
    for (size_t i = 0; i < k_indices.size (); ++i)
    {
      k_indices[i] = index_mapping_[k_indices[i]];
    }
  }

  template <typename PointVertexType>
  int KDTreeFlannT<PointVertexType>::radiusSearch (const PointVertexType* vertex, double radius, std::vector<int> &k_indices,
                                          std::vector<float> &k_squared_distances, int max_nn) const
  {

    if(index_mapping_.size() == 0){
      std::cerr << "Warning: radiusSearch on an empty kDTree." << std::endl;
      return 0;
    }

    //limit maximum number of neighbors
    if(max_nn <= 0){ 
      max_nn = 1000;
    }

    //copy and cast 3d pose 
    std::vector<float> tmp (dim_);
    KDTreeFlannT<PointVertexType>::vectorize(dim_, vertex, tmp);
    radius *= radius; // flann uses squared radius

    int neighbors_in_radius = 0;
    k_indices.resize(1000);
    k_squared_distances.resize(1000);

    // if using preallocated vectors we ignore max_nn as we are sure to have enought space
    // to store all neighbors found in radius
      flann::Matrix<int> k_indices_mat (&k_indices[0], 1, k_indices.size());
      flann::Matrix<float> k_distances_mat (&k_squared_distances[0], 1, k_squared_distances.size());
      neighbors_in_radius = flann_index_->radiusSearch (flann::Matrix<float>(&tmp[0], 1, dim_),
          k_indices_mat, k_distances_mat, radius, flann::SearchParams (-1, 0.0001, true));

    if(neighbors_in_radius > max_nn){
      neighbors_in_radius = max_nn;
    }

    if (neighbors_in_radius == 0)
    {
      return (0);
    }

    // Do mapping to vertex indices
    for (int i = 0; i < neighbors_in_radius; ++i)
    {
      int& neighbor_index = k_indices[i];
      neighbor_index = index_mapping_[neighbor_index];
    }

    return (neighbors_in_radius);
  }


  template <typename PointVertexType>
  int KDTreeFlannT<PointVertexType>::neighborsInRadius(const PointVertexType* vertex, double radius) const
  {

    if(index_mapping_.size() == 0){
      std::cerr << "Warning: neighborsInRadius on an empty kDTree." << std::endl;
      return;
    }

    static flann::Matrix<int> indices_empty;
    static flann::Matrix<float> dists_empty;

    //copy and cast 3d pose 
    std::vector<float> tmp (dim_);
    KDTreeFlannT<PointVertexType>::vectorize(dim_, vertex, tmp);
    radius *= radius; // flann uses squared radius


    int neighbors_in_radius = flann_index_->radiusSearch (flann::Matrix<float>(&tmp[0], 1, dim_),
          indices_empty, dists_empty, radius, flann::SearchParams (-1, 0.0001, false));
    return neighbors_in_radius;
  }

  template <typename PointVertexType>
  void KDTreeFlannT<PointVertexType>::clear()
  {

    if(flann_index_ != 0)
    {
      delete flann_index_;
      flann_index_ = 0;
    }

    // Data array cleanup
    if(index_mapping_.size() > 0){
      free (data_);
      data_ = NULL;
      index_mapping_.clear();
    }

  }

  template <typename PointVertexType>
  void KDTreeFlannT<PointVertexType>::vectorize(unsigned int dim, const PointVertexType* vertex, std::vector<float>& vec)
  {
    unsigned int vertexDim = PointVertexType::Dimension;
    for(size_t i = 0; i < dim; ++i){
      if(i < vertexDim)
        vec[i] = (float) vertex->estimate()[i];
//       if(i >= vertexDim && i < 2*vertexDim)
//         vec[i] = (float) vertex->normal()[i-3];
    }
  }

  template <typename PointVertexType>
  void KDTreeFlannT<PointVertexType>::vectorize(unsigned int dim, const Eigen::Vector3f& point, std::vector<float>& vec)
  {
    for(size_t i = 0; i < dim; ++i)
      vec[i] = point(i);
  }

}
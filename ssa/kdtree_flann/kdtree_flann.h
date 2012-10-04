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

#ifndef __SSA_KDTREE_FLANN__
#define __SSA_KDTREE_FLANN__

#include <vector>
#include "ssa/core/ssa_graph.h"

#pragma GCC diagnostic ignored "-Wunused-parameter"  // Do not show warnings from PCL or FLANN
#include <flann/flann.hpp>
#pragma GCC diagnostic warning "-Wunused-parameter"
#include "g2o/stuff/timeutil.h"

namespace ssa {

  template <typename PointVertexType>
  class KDTreeFlannT{

    public:
    KDTreeFlannT();
    KDTreeFlannT(int dim);
    ~KDTreeFlannT();

    /** \brief copies point information from vertices to internal FLANN structure
        \param useVectorIndices 
         true: the index of the input vector is returned during search 
         false: the vertex id is returned during search
    */
    void copyData(std::vector<PointVertexType* >& vertices, bool useVectorIndices);

    void createKDTree();

    void nearestKSearch(const PointVertexType* vertex, int k, 
                                       std::vector<int> &k_indices, 
                                       std::vector<float> &k_squared_distances);

    void nearestKSearch(const Eigen::Vector3f& point, int k, 
                                       std::vector<int> &k_indices, 
                                       std::vector<float> &k_squared_distances);

    int radiusSearch(const PointVertexType* vertex, double radius, std::vector<int> &k_indices,
                                          std::vector<float> &k_squared_distances, int max_nn) const;

    int neighborsInRadius(const PointVertexType* vertex, double radius) const;

    void clear();

    static void vectorize(unsigned int dim, const PointVertexType* vertex, std::vector<float>& vec);

    static void vectorize(unsigned int dim, const Eigen::Vector3f& point, std::vector<float>& vec);

    /** \brief mapping between internal and external indices. */
    std::vector<int> index_mapping_;

    /** \brief Internal pointer to data. */
    float* data_;

    /** \brief Tree dimensionality (i.e. the number of dimensions per point). */
    int dim_;

    /** \brief A FLANN index object. */
    flann::Index<flann::L2_Simple<float> >* flann_index_;

    /** \brief descides which data to use for kdtree 
                0 metric infomation
                1 normals
                2 color
    */
    int mode_;

  };
}

#include "kdtree_flann.hpp"

#endif

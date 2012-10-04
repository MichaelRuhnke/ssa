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

#ifndef __SPARSE_SURFACE_ADJUSTMENT__
#define __SPARSE_SURFACE_ADJUSTMENT__

//#define NUMERIC_JACOBIAN_TWO_D_TYPES

#include <vector>
#include "g2o/stuff/timeutil.h"

#include "g2o/core/solver.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/optimization_algorithm_levenberg.h"

#include "ssa/data_association/data_association.h"
#include "ssa/data_association/normal_shooting_flann.h"
#include "ssa/core/parameter.h"

namespace ssa {

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  class SparseSurfaceAdjustmentT{

    typedef SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3> SSAGraph;
    typedef typename EdgeType2::VertexXiType   PoseVertex;
    typedef typename EdgeType2::VertexXjType   PointVertex;
    static const int Dj = EdgeType2::Dimension;
    typedef Eigen::Matrix<double, Dj, 1> PointVector;
    public:

    SparseSurfaceAdjustmentT();
    ~SparseSurfaceAdjustmentT();

    /** \brief inserts an existing graph */ 
    void setGraph(SSAGraph& graph);

    /** \brief direct access to graph */ 
    SSAGraph* graph();

    /** \brief set the sparse linear solver for optimization */ 
    void setSolver(g2o::BlockSolverX::LinearSolverType*& linearSolver);

    /** \brief show verbose outputs */ 
    void setVerbose(bool verbose);

    /** \brief set ssa params*/ 
    void setParams(SparseSurfaceAdjustmentParams& params);

    /** \brief direct access to params */ 
    SparseSurfaceAdjustmentParams& params();

    /** \brief perform optimization */
    void optimize(int level = 0);

    /** \brief perform optimization */
    void optimizeHierarchical(int startLevel = 0);

    /** \brief Optimize color information with mean filter */
    void optimizeColors();

    /** \brief Optimize color information with weighted mean filter. 
        \brief Weights are calculated from the incidence angle.
    */
    void optimizeColorsIncidenteWeight();

    void dumpConnectivityMatrix();

    SSAGraph graph_;

    protected:
      bool verbose_;
      SparseSurfaceAdjustmentParams params_;
  };

} //end namespace

#include "sparse_surface_adjustment.hpp"

#endif
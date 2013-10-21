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

#ifndef __SPARSE_SURFACE_ADJUSTMENT_GRAPH_3D__
#define __SPARSE_SURFACE_ADJUSTMENT_GRAPH_3D__

//Include vertex and edge types
#include "g2o/types/slam3d/types_slam3d.h"

#include "ssa/types_3d/vertex_point_xyzcov.h"
#include "ssa/types_3d/edge_se3_xyzcov.h"
#include "ssa/types_3d/edge_xyzcov_xyzcov.h"

#include "ssa/types_3d/vertex_point_xyz_normal.h"
#include "ssa/types_3d/edge_se3_point_xyz_normal.h"
#include "ssa/types_3d/edge_xyznormal_xyznormal.h"

//Include ssa_graph template
#include "ssa/core/ssa_graph.h"

namespace ssa{
  typedef SparseSurfaceAdjustmentGraphT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov> SparseSurfaceAdjustmentGraph3D;
  typedef SparseSurfaceAdjustmentGraphT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZNormal, ssa::EdgePointXYZNormalPointXYZNormal> SparseSurfaceAdjustmentGraphRGBA;
}
#endif
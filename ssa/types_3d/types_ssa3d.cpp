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

#include "vertex_point_xyzcov.h"
#include "edge_se3_xyzcov.h"
#include "edge_xyzcov_xyzcov.h"
#include "g2o/core/factory.h"
#include "g2o/stuff/macros.h"

namespace ssa {

  G2O_REGISTER_TYPE_GROUP(ssa3d);

  G2O_REGISTER_TYPE(VERTEX_XYZCOV, VertexPointXYZCov);
  G2O_REGISTER_TYPE(EDGE_SE3_XYZCOV, EdgeSE3PointXYZCov);
  G2O_REGISTER_TYPE(EDGE_XYZCOV_XYZCOV, EdgePointXYZCovPointXYZCov);

} // end namespace

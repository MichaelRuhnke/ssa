// Sparse Surface Optimization 2D
// Copyright (C) 2011 M. Ruhnke, R. Kuemmerle, G. Grisetti, W. Burgard
// 
// SSA2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// SSA2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "vertex_point_xycov.h"
#include "edge_se2_xycov.h"
#include "edge_xycov_xycov.h"
#include "g2o/core/factory.h"
#include "g2o/stuff/macros.h"

namespace ssa {

  G2O_REGISTER_TYPE_GROUP(ssa2d);

  G2O_REGISTER_TYPE(VERTEX_XYCOV, VertexPointXYCov);
  G2O_REGISTER_TYPE(EDGE_SE2_XYCOV, EdgeSE2PointXYCov);
  G2O_REGISTER_TYPE(EDGE_XYCOV_XYCOV, EdgePointXYCovPointXYCov);

//   void G2O_ATTRIBUTE_CONSTRUCTOR init_types_ssa(void)
//   {
//     Factory* factory = Factory::instance();
// 
//     factory->registerType("VERTEX_XYCOV", new HyperGraphElementCreator<ssa::VertexPointXYCov>);
//     factory->registerType("EDGE_SE2_XYCOV", new HyperGraphElementCreator<ssa::EdgeSE2PointXYCov>);
//     factory->registerType("EDGE_XYCOV_XYCOV", new HyperGraphElementCreator<ssa::EdgePointXYCovPointXYCov>);
// 
//   }

} // end namespace

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

#ifndef __SSA_IMPORT_POSEGRAPH__
#define __SSA_IMPORT_POSEGRAPH__

#include <istream>
#include <fstream>
#include <sstream>
#include <string>

#include "ssa/core/ssa_graph_2d.h"
//forward declaration
namespace g2o{
  class VertexSE2;
}

namespace ssa {
  //forward declarations

  class SSAPoseGraph2D {
  public:
  SSAPoseGraph2D();
  ~SSAPoseGraph2D();

  static void importPoseGraph2D(std::string logfile, SparseSurfaceAdjustmentGraph2D& ssaGraph);
  static void importRobotLaser(std::istream& is, g2o::VertexSE2* parent, SparseSurfaceAdjustmentGraph2D& ssaGraph);

  };

} //end namespace

#endif
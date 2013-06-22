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

//#include <Eigen/Eigenvalues>
#include "observation_point_xyz_rgb_normal.h"

namespace ssa {
  using namespace Eigen;

  ObservationXYZRGBNormal::ObservationXYZRGBNormal() : measurement(Eigen::Vector3d::Zero()), estimate(Eigen::Vector3d::Zero()), normal(Eigen::Vector3d::Zero()), rgb(Eigen::Vector3i::Zero()), level(0)
  {

  }

  ObservationXYZRGBNormal::ObservationXYZRGBNormal(Eigen::Vector3d& m, Eigen::Vector3d& e, Eigen::Vector3d& n, Eigen::Vector3i& color, int l) : measurement(m), estimate(e), normal(n), rgb(color), level(l)
  {

  };


  bool ObservationXYZRGBNormal::read(std::istream& is)
  {
    /// measurement in sensor frame
    for(int i = 0; i < 3; ++i)
      is >> measurement(i);

   /// estimate in global frame
    for(int i = 0; i < 3; ++i)
      is >> estimate(i);

    /// normal in local frame
    for(int i = 0; i < 3; ++i)
      is >> normal(i);

    /// color as r,g,b
    for(int i = 0; i < 3; ++i)
      is >> rgb(i);

    is >> level;
    return true;
  }

  bool ObservationXYZRGBNormal::write(std::ostream& os) const
  {
    /// measurement in sensor frame
    for(int i = 0; i < 3; ++i)
      os << measurement(i) << " ";

    /// estimate in global frame
    for(int i = 0; i < 3; ++i)
      os << estimate(i) << " ";

    /// normal in local frame
    for(int i = 0; i < 3; ++i)
      os << normal(i) << " ";

    /// color as r,g,b
    for(int i = 0; i < 3; ++i)
      os << rgb(i) << " ";

    os << level << " ";
    return os.good();
  }

} //end namespace


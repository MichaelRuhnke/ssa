#ifndef __SSA_OBSERVATION_XYZ_RGB_NORMAL__
#define __SSA_OBSERVATION_XYZ_RGB_NORMAL__

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

#include <Eigen/Core>

namespace ssa {

  /** \brief: SSA class for single observation of a surface point
   *
   * This point vertices represent the observed surface,
   * which will be adapted during optimization.
   *
   */

  class ObservationXYZRGBNormal
  {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

      Eigen::Vector3d measurement;
      Eigen::Vector3d estimate;
      Eigen::Vector3d normal;
      Eigen::Vector3i rgb;
      unsigned short level;


      ObservationXYZRGBNormal();
      ObservationXYZRGBNormal(Eigen::Vector3d& measurement, Eigen::Vector3d& estimate, Eigen::Vector3d& normal, Eigen::Vector3i& rgb, int level = 0);

      bool read(std::istream& is);
      bool write(std::ostream& os) const;

  };

}
#endif

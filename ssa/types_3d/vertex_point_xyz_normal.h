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

#ifndef __SSA_VERTEX_XYZ_NORMAL_3D__
#define __SSA_VERTEX_XYZ_NORMAL_3D__

#include <Eigen/Geometry>
#include "g2o/core/base_vertex.h"
#include "g2o/core/hyper_graph_action.h"
#include "g2o/types/slam3d/vertex_pointxyz.h"

//forward declaration
namespace g2o {
  class VertexSE3;
}

namespace ssa {
  //forward declaration
  class EdgeSE3PointXYZNormal;

  /** \brief: Vertex class for SSA surface points
   *
   * This point vertices represent the observed surface,
   * which will be adapted during optimization.
   *
   */

  class VertexPointXYZNormal : public g2o::VertexPointXYZ
  {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      explicit VertexPointXYZNormal();

      VertexPointXYZNormal(int id, Eigen::Vector3d& estimate, Eigen::Vector3d& normal);

      virtual bool read(std::istream& is);
      virtual bool write(std::ostream& os) const;

      /**
       * Getter for parent vertex pointer
       *
       * @return pointer to sensor pose
       */
      g2o::VertexSE3* parentVertex();
      g2o::VertexSE3* parentVertex() const;

      inline int& parentVertexId(){ return _parentVertexId; }

      /** Get pointer to parent vertices (makes sense in pruned graphs) */
      std::vector< g2o::VertexSE3* > parentVertices() const;

      /** Set the pointer to parent vertex */
      void setParentVertex(g2o::VertexSE3* pose);

      /** Normal in scan/observation coordinate frame */
      Eigen::Vector3d& normal();
      Eigen::Vector3d  normal() const;
      /** Normal in global coordinate frame */
      Eigen::Vector3d globalNormal() const;

      /** rgb color information per point */
      unsigned char cr;
      unsigned char cg;
      unsigned char cb;



    protected:

      g2o::VertexSE3* _parentVertex; /** Sensor pose vertex from which this point was observed */
      int             _parentVertexId;
      Eigen::Matrix3d _cov;     /** Covariance in sensor frame */
      Eigen::Vector3d _normal;      /** normal in sensor frame (not saved to disk) */
  };

}
#endif

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

#include "edge_xyznormal_xyznormal.h"

namespace ssa {
 using namespace Eigen;
 using namespace g2o;
 EdgePointXYZNormalPointXYZNormal::EdgePointXYZNormalPointXYZNormal() :
    BaseBinaryEdge<3, Eigen::Vector3d, VertexPointXYZNormal, VertexPointXYZNormal>()
  {
    identity.setIdentity();
  }

  void EdgePointXYZNormalPointXYZNormal::computeError()
  {
    const VertexPointXYZNormal* l1 = static_cast<const VertexPointXYZNormal*>(_vertices[0]);
    const VertexPointXYZNormal* l2 = static_cast<const VertexPointXYZNormal*>(_vertices[1]);
    _error = l1->estimate() - l2->estimate();

    /** this is the way Kurt computes the information matrix
      // re-define the information matrix
      // topLeftCorner<3,3>() is the rotation()
      const Matrix3d transform = ( vp0->estimate().inverse() *  vp1->estimate() ).matrix().topLeftCorner<3,3>();
      information() = ( cov0 + transform * cov1 * transform.transpose() ).inverse();
    **/
    Eigen::Vector3d normal = l1->globalNormal() + l2->globalNormal();
    normal.normalize();

    Eigen::Matrix3d rot;
    rot.col(2) = normal;
    double angleX = std::min(normal.dot(identity.row(0)), normal.dot(-1.0*identity.row(0)));
    double angleY = std::min(normal.dot(identity.row(1)), normal.dot(-1.0*identity.row(1)));
    double angleZ = std::min(normal.dot(identity.row(2)), normal.dot(-1.0*identity.row(2)));

    if((angleY <= angleX && angleY <= angleZ)){
      rot.col(0) = (identity.row(0).cross(normal)).normalized();
    } else {
      rot.col(0) = (identity.row(1).cross(normal)).normalized();
    }
    rot.col(1) = (normal.cross(rot.col(0))).normalized();

    information().setIdentity();
    information()(2,2) = 100000;
    information() = rot.transpose()*information()*rot;
  }

  bool EdgePointXYZNormalPointXYZNormal::read(std::istream& is)
  {
    /** the measurement should always be 0,0,0 and the information matrix is computed in computeError() */
    Eigen::Vector3d p(0.0,0.0,0.0);
//     is >> p[0] >> p[1] >> p[2];
    setMeasurement(p);
//
//     for (int i=0; i<3; i++)
//       for (int j=i; j<3; j++) {
//         is >> information()(i,j);
//         if (i!=j)
//           information()(j,i)=information()(i,j);
//       }
    return true;
  }

  bool EdgePointXYZNormalPointXYZNormal::write(std::ostream& os) const
  {
    /** the measurement should always be 0,0,0 and the information matrix is computed in computeError() */
//     for (int i=0; i<3; i++)
//       os << measurement()[i] << " ";
//     for (int i=0; i<3; i++)
//       for (int j=i; j<3; j++){
//         os << " " <<  information()(i,j);
//       }
    return os.good();
  }

  void EdgePointXYZNormalPointXYZNormal::initialEstimate(const OptimizableGraph::VertexSet& from, OptimizableGraph::Vertex* to)
  {
    assert(from.size() == 1 && from.count(_vertices[0]) == 1 && "Can not initialize VertexSE3 position by VertexPointXYZNormal");
    VertexPointXYZNormal* vi = static_cast<VertexPointXYZNormal*>(_vertices[0]);
    VertexPointXYZNormal* vj = static_cast<VertexPointXYZNormal*>(_vertices[1]);
    if (from.count(vi) > 0 && to == vj) {
      vj->setEstimate(vi->estimate());
    }
  }

#ifndef NUMERIC_JACOBIAN_THREE_D_TYPES
  void EdgePointXYZNormalPointXYZNormal::linearizeOplus()
  {
    _jacobianOplusXi( 0 , 0 ) = 1;
    _jacobianOplusXi( 0 , 1 ) = 0;
    _jacobianOplusXi( 0 , 2 ) = 0;
    _jacobianOplusXi( 1 , 0 ) = 0;
    _jacobianOplusXi( 1 , 1 ) = 1;
    _jacobianOplusXi( 1 , 2 ) = 0;
    _jacobianOplusXi( 2,  0 ) = 0;
    _jacobianOplusXi( 2 , 1 ) = 0;
    _jacobianOplusXi( 2 , 2 ) = 1;

    _jacobianOplusXj( 0 , 0 ) = -1;
    _jacobianOplusXj( 0 , 1 ) = 0;
    _jacobianOplusXj( 0 , 2 ) = 0;
    _jacobianOplusXj( 1 , 0 ) = 0;
    _jacobianOplusXj( 1 , 1 ) = -1;
    _jacobianOplusXj( 1 , 2 ) = 0;
    _jacobianOplusXj( 2 , 0 ) = 0;
    _jacobianOplusXj( 2 , 1 ) = 0;
    _jacobianOplusXj( 2 , 2 ) = -1;
  }
#endif

} //end namespace
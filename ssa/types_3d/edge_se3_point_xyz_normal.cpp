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

#include "edge_se3_point_xyz_normal.h"

namespace ssa {
  using namespace Eigen;
  using namespace g2o;

  EdgeSE3PointXYZNormal::EdgeSE3PointXYZNormal() :
    g2o::BaseBinaryEdge<3, Eigen::Vector3d, g2o::VertexSE3, VertexPointXYZNormal>() , level_(0)
  {
  }

  EdgeSE3PointXYZNormal::~EdgeSE3PointXYZNormal()
  {
  }

  bool EdgeSE3PointXYZNormal::read(std::istream& is)
  {
    /** read measurement*/
    Eigen::Vector3d p;
    is >> p[0] >> p[1] >> p[2];
    setMeasurement(p);

    /** read information matrix */
    for (int i=0; i<3; i++)
      for (int j=i; j<3; j++) {
        is >> information()(i,j);
        if (i!=j)
          information()(j,i)=information()(i,j);
      }

    /** read level of edge */
    int l = 0;
    is >> l;
    setLevel(l);

    return true;
  }

  bool EdgeSE3PointXYZNormal::write(std::ostream& os) const
  {
    /** write original measurement of sensor */
    for (int i=0; i<3; i++)
      os << measurement()[i] << " ";

    /** write sensor dependent information matrix */
    for (int i=0; i<3; i++)
      for (int j=i; j<3; j++){
        os <<  information()(i,j) << " ";
      }
    /** write level of edge (used for hierarchical ssa stuff)*/
    os << level_;

    return os.good();
  }

  void EdgeSE3PointXYZNormal::initialEstimate(const g2o::OptimizableGraph::VertexSet& from_, g2o::OptimizableGraph::Vertex* /*to_*/)
  {
    VertexSE3*       from = static_cast<VertexSE3*>(_vertices[0]);
    VertexPointXYZNormal* to = static_cast<VertexPointXYZNormal*>(_vertices[1]);
    if (from_.count(from) > 0)
      to->setEstimate(from->estimate() * measurement());
    else
      std::cerr << __PRETTY_FUNCTION__ << std::endl;
  }

}


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
#include "g2o/types/slam3d/types_slam3d.h"
#include "pcl/common/eigen.h"
#include "vertex_point_xyz_normal.h"
#include "edge_se3_point_xyz_normal.h"

namespace ssa {
  using namespace Eigen;
  using namespace g2o;

  VertexPointXYZNormal::VertexPointXYZNormal() : VertexPointXYZ(),
    _parentVertex(0),
    _normal(Eigen::Vector3d::Zero())
  {
    _id = 0;
  };

  VertexPointXYZNormal::VertexPointXYZNormal(int id, Eigen::Vector3d& estimate, Eigen::Vector3d& normal) : VertexPointXYZ(),
    _parentVertex(0),
    _normal(normal)
  {
    _id = id;
    _estimate = estimate;
  }

  bool VertexPointXYZNormal::read(std::istream& is)
  {
    ///read id of sensor pose vertex
    is >> _parentVertexId;

    ///read point pose
    Eigen::Vector3d p;
    is >> p[0] >> p[1] >> p[2];
    setEstimate(p);

    ///read normal
    is >> _normal(0) >> _normal(1) >> _normal(2);

    ///read color information
    uint r,g,b;
    is >> r >> g >> b;
    cr = (unsigned char) r;
    cg = (unsigned char) g;
    cb = (unsigned char) b;

    return true;
  }

  bool VertexPointXYZNormal::write(std::ostream& os) const
  {
    /// we store the id of the observation / parent vertex for faster lookups...
    /// this information is also available through the connected edges
    os << _parentVertex->id() << " ";

    /// estimated pose of point (global coordinate frame):
    /// this will be different from the measured depth after optimization!
    for (int i=0; i<3; i++)
      os << estimate()[i] << " ";

    /// normal of point (local coordinate frame):
    for (int i=0; i<3; i++)
      os << _normal(i) << " ";

    /// color information for the point (if available)
    uint r,g,b;
    r = (uint) cr;
    g = (uint) cg;
    b = (uint) cb;
    os << r << " " << g << " " << b << " ";

    return os.good();
  }

    g2o::VertexSE3* VertexPointXYZNormal::parentVertex() {
    if(_parentVertex != 0)
      return _parentVertex;
    for (g2o::OptimizableGraph::EdgeSet::iterator it=edges().begin(); it!=edges().end(); it++){
      EdgeSE3PointXYZNormal* e2=dynamic_cast<EdgeSE3PointXYZNormal*>(*it);
      if(e2){
        g2o::VertexSE3* p = dynamic_cast<g2o::VertexSE3* >(e2->vertices()[0]);
        if(p){
          _parentVertex = p;
        }
        return p;
      }
    }
    g2o::VertexSE3* p = 0;
    return p;
  }

  g2o::VertexSE3* VertexPointXYZNormal::parentVertex() const{
    if(_parentVertex != 0)
      return _parentVertex;
    for (g2o::OptimizableGraph::EdgeSet::iterator it=edges().begin(); it!=edges().end(); it++){
      EdgeSE3PointXYZNormal* e2=dynamic_cast<EdgeSE3PointXYZNormal*>(*it);
      if(e2){
        g2o::VertexSE3* p = dynamic_cast<g2o::VertexSE3* >(e2->vertices()[0]);
        return p;
      }
    }
    g2o::VertexSE3* p = 0;
    return p;
  }

  std::vector< g2o::VertexSE3* > VertexPointXYZNormal::parentVertices() const{
    std::vector< g2o::VertexSE3* > parents;
    for (g2o::OptimizableGraph::EdgeSet::iterator it=edges().begin(); it!=edges().end(); it++){
      ssa::EdgeSE3PointXYZNormal* e2=dynamic_cast<ssa::EdgeSE3PointXYZNormal*>(*it);
      if(e2){
        g2o::VertexSE3* p = dynamic_cast<g2o::VertexSE3* >(e2->vertices()[0]);
        if(p)
          parents.push_back(p);
      }
    }
    return parents;
  }

  void VertexPointXYZNormal::setParentVertex(VertexSE3* parent){
    _parentVertex = parent;
  }

  Eigen::Vector3d& VertexPointXYZNormal::normal(){
    return _normal;
  }

  Eigen::Vector3d VertexPointXYZNormal::globalNormal() const{
    Vector3d normal(0.0, 0.0, 0.0);
    if(parentVertex() > 0)
	    normal = parentVertex()->estimate().linear() * _normal;
    return normal;
  }

  Vector3d VertexPointXYZNormal::normal() const{
    return _normal;
  }

} //end namespace


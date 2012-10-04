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


#ifndef __SSA_EDGE_POINTXYCOV_2D_POINTXYCOV_2D__
#define __SSA_EDGE_POINTXYCOV_2D_POINTXYCOV_2D__

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include "vertex_point_xycov.h"
#include "g2o/types/slam2d/vertex_se2.h"
#include "g2o/core/base_binary_edge.h"

namespace ssa {
  using namespace Eigen;
  using namespace g2o;

class EdgePointXYCovPointXYCov : public BaseBinaryEdge<2, Eigen::Vector2d, VertexPointXYCov, VertexPointXYCov>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgePointXYCovPointXYCov();

    void computeError();
    
    static bool DataAssociationEdgeComp (EdgePointXYCovPointXYCov* i,EdgePointXYCovPointXYCov* j) { if(i && j){return (i->chi2()<j->chi2());} else {std::cerr << "DataAssociationEdgeComp: WARNING at least one edge not valid! This might happen due to wrong casting or edge deletion!"; return false;} }

    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;

    virtual void initialEstimate(const OptimizableGraph::VertexSet& from, OptimizableGraph::Vertex* to);
    virtual double initialEstimatePossible(const OptimizableGraph::VertexSet& from, OptimizableGraph::Vertex* to) { (void) to; return (from.count(_vertices[0]) == 1 ? 1.0 : -1.0);}
#ifndef NUMERIC_JACOBIAN_TWO_D_TYPES
    virtual void linearizeOplus();
#endif
    bool isMeshEdge;
};

}

// AIS_REGISTER_VERTEX("LANDMARK2BA", AISNavigation::VertexPointXYCov);
// AIS_REGISTER_EDGE("LANDMARKEDGE2BA", AISNavigation::EdgeSE2PointXYCov);
// AIS_REGISTER_EDGE("LAND2LANDEDGE2BA", AISNavigation::EdgePointXYCovPointXYCov);

#endif

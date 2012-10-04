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

#ifndef __SSA_VERTEX_XY_COV_2D__
#define __SSA_VERTEX_XY_COV_2D__
#include <Eigen/Geometry>

#include "g2o/core/base_vertex.h"
#include "g2o/core/hyper_graph_action.h"

// #define NUMERIC_JACOBIAN_TWO_D_TYPES

//forward declaration
namespace g2o {
  class VertexSE2;
}

namespace ssa {

  class VertexPointXYCov : public g2o::BaseVertex<2, Eigen::Vector2d>
  {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      VertexPointXYCov();

      virtual void setToOriginImpl() {
        _estimate.setZero();
      }

      virtual bool setEstimateDataImpl(const double* est){
        _estimate[0] = est[0];
        _estimate[1] = est[1];
        return true;
      }

      virtual bool getEstimateData(double* est) const{
        est[0] = _estimate[0];
        est[1] = _estimate[1];
        return true;
      }
      
      virtual int estimateDimension() const { 
        return 2;
      }

      virtual bool setMinimalEstimateDataImpl(const double* est){
        return setEstimateData(est);
      }

      virtual bool getMinimalEstimateData(double* est) const{
        return getEstimateData(est);
      }
      
      virtual int minimalEstimateDimension() const { 
        return 2;
      }

      virtual void oplusImpl(const double* update)
      {
        _estimate[0] += update[0];
        _estimate[1] += update[1];
      }

      virtual bool read(std::istream& is);
      virtual bool write(std::ostream& os) const;

      /** parent pointer TODO: could be removed... and defined edge dependent*/
    private:
      g2o::VertexSE2* _parentVertex;
      unsigned int _parentVertexId;
    public:
      g2o::VertexSE2* parentVertex() const;
      unsigned int parentVertexId() const;
      void setParentVertex(g2o::VertexSE2* pose);

      /** normal in sensor frame (not saved to disk) */
    private:
      Eigen::Vector2d    _normal; 
  
    public:
      Eigen::Vector2d& normal();
      Eigen::Vector2d normal() const;
      Eigen::Vector2d globalNormal();
      
      /** update normal and covariance based on point neighborhoods covariance (experimental) */
      void updateNormal(Eigen::Matrix2d& cov);
      
      bool            _hasNormal;
  
    /** Covariance in sensor frame */
    private:
      Eigen::Matrix2d _cov; 
    public: 
      Eigen::Matrix2d covariance() const;
      Eigen::Matrix2d& covariance();

    std::deque< VertexPointXYCov* >  neighbors;

    /** color information per point */
    unsigned char cr;
    unsigned char cg;
    unsigned char cb;

    /** ratio between eigenvalues */
    double ratio;
  };

}
#endif

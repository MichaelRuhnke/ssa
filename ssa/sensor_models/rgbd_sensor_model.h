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

#ifndef __SSA_RGBD_SENSOR_MODEL__
#define __SSA_RGBD_SENSOR_MODEL__

#include <Eigen/Core>

namespace ssa{

class RGBDSensorModel{

  public:

  inline RGBDSensorModel() : q_pix(8), baseline(0.075), focal_length(580)
  {
    ///default values are for Microsoft Kinects taken from http://www.ros.org/wiki/kinect_calibration/technical
    factor = q_pix * baseline * focal_length;
  };

  inline RGBDSensorModel(int q, double b, int f) : q_pix(q), baseline(b), focal_length(f)
  {
    factor = q_pix * baseline * focal_length;
  };

  inline ~RGBDSensorModel()
   {
//       std::cerr << __PRETTY_FUNCTION__ << "plop" << std::endl;
   };

  inline double depthUncertainty(double range){
    double result = (0.5 * factor * fabs((1.0 / ceil((factor/range) + 0.5)) - (1.0 / floor((factor/range) - 0.5))));
//     double depthUncertainty = 0.0016 * (range * range) + 0.0001; ///quadratic function used as approximation
//     if(fabs(result - depthUncertainty) > 0.00023)
//       std::cerr << fixed << result  << " " << depthUncertainty << " " << result - depthUncertainty << std::endl;
    return result;
  };

  inline void getInformationMatrix(Eigen::Vector3d& beam, Eigen::Matrix3d& information)
  {
      ///assume constant x,y error
      double cError = 0.001;
      double range = beam(2);

      double imagePlaneUncertainty = (range * cError);
      double depthUn = depthUncertainty(range);

      ///construct information matrix
      information = Eigen::Matrix3d::Identity();
      information(0,0) = 1.0 / (imagePlaneUncertainty * imagePlaneUncertainty);
      information(1,1) = information(0,0);
      information(2,2) = 1.0  / (depthUn*depthUn);

      ///rotate covariance
      Eigen::Isometry3d transformation;
      Eigen::Vector3d y_direction(0.0, 1.0, 0.0);

      Eigen::Vector3d tmp0 = (y_direction.cross(beam)).normalized();
      Eigen::Vector3d tmp1 = (beam.cross(tmp0)).normalized();
      Eigen::Vector3d tmp2 = beam.normalized();

      transformation(0,0)=tmp0[0]; transformation(0,1)=tmp0[1]; transformation(0,2)=tmp0[2]; transformation(0,3)=0.0f;
      transformation(1,0)=tmp1[0]; transformation(1,1)=tmp1[1]; transformation(1,2)=tmp1[2]; transformation(1,3)=0.0f;
      transformation(2,0)=tmp2[0]; transformation(2,1)=tmp2[1]; transformation(2,2)=tmp2[2]; transformation(2,3)=0.0f;
      transformation(3,0)=0.0f;    transformation(3,1)=0.0f;    transformation(3,2)=0.0f;    transformation(3,3)=1.0f;

      information = transformation.linear().matrix() * information * transformation.linear().matrix().transpose();
  }

  protected:
    double  factor;
    int     q_pix;
    double  baseline;
    int     focal_length;
};

}

#endif
#include <Eigen/Geometry>
#include "g2o/types/slam3d/se3quat.h"

namespace ssa
{

  inline g2o::SE3Quat Affine3fToSE3Quat(Eigen::Affine3f transformation){
    Eigen::Affine3d transformationd = transformation.cast<double>();
    Eigen::Quaterniond q;
    q = transformationd.rotation();
    SE3Quat optimizerTrans = SE3Quat(q, transformationd.translation());
    return optimizerTrans;
  }

  inline Eigen::Affine3f 
  SE3QuatToAffine3f(SE3Quat optimizerTrans){
    Eigen::Vector3f t = optimizerTrans.translation().cast<float>();
    Eigen::Quaternionf q = optimizerTrans.rotation().cast<float>();
    Eigen::Affine3f transformation = (Eigen::Affine3f)Eigen::Translation3f(t) * q;
    return transformation;
  }

}

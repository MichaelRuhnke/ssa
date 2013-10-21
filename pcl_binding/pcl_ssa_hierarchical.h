#ifndef __PCL_SSA_HIERARCHICAL__
#define __PCL_SSA_HIERARCHICAL__

#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "pcl/point_types.h"
#include "pcl/point_cloud.h"
#include "pcl/io/pcd_io.h"
#include "pcl/common/transforms.h"
#include "pcl/filters/filter.h"
#include "pcl/filters/voxel_grid.h"

#include "ssa/core/ssa_graph_3d.h"
#include "ssa/sensor_models/rgbd_sensor_model.h"
#include "ssa/core/parameter.h"
#include "g2o/stuff/string_tools.h"

#include "pcl_g2o_math_conversion.h"

// #warning "pcl_ssa_hierarchical included"
// Forward declarations
namespace pcl {
  class RangeImage;
}

namespace g2o{
  class EdgeSE3PointXYZCov;
  class EdgePointXYZCovPointXYZCov;
}

namespace ssa{

  template <typename PointCloudType>
  class PCLSSAHierarchicalT {
    typedef Eigen::Matrix<double, 6, 1> Vector6d;
    typedef typename PointCloudType::Ptr PointCloudPtr;

    public:

    PCLSSAHierarchicalT ();
    ~PCLSSAHierarchicalT ();

    enum Sensor{KINECT, LMS};

    /** \brief sets the input point cloud vector from wich the ssa graph will be constructed
        Make shure that each cloud has a properly set sensor_origin_ and sensor_orientation_
    */
    void
    setInput (std::vector< PointCloudType >* clouds);

    /** \brief sets the type of sensor from which the data are measured */
    void
    setSensor (Sensor sensor);

    /** \brief sets the number of hierarchy levels */
    void
    setLevels (int levels);

    /** \brief sets the resolution of one hierarchy level.
        Note level 0 is the lowest level with the highest resolution.
    */
    void
    setResolution (int level, double resolution);

    inline
    SparseSurfaceAdjustmentParams&
    params () { return params_; };

    /** \brief constructs ssa graph and optimizes the system. */
    void
    optimize (PointCloudType& result);

    /** \brief constructs ssa graph. This will be called from optimize. */
    void
    constructGraph (SparseSurfaceAdjustmentGraph3D& graph);

    /** \brief inserts cloud into ssa graph. */
    int
    insertCloud(SparseSurfaceAdjustmentGraph3D& graph, PointCloudType& cloud);

    int addCloudToGraph(SparseSurfaceAdjustmentGraph3D& graph, typename PointCloudType::Ptr cloud);

    private:
      std::vector< PointCloudType >* input_clouds_;
      Sensor                                            sensor_;
      int                                               levels_;
      std::map<int, double>                             resolution_;
      SparseSurfaceAdjustmentParams                     params_;

    public:

    /** Returns pcl pointcloud for landmark edges vector.*/
    static PointCloudType landmarksToPointCloud(std::vector<ssa::EdgeSE3PointXYZCov* >& landmarksFromVertex, bool useRawMeasurements = false);

    /** Returns pcl pointcloud for landmark edges vector.*/
    static PointCloudType landmarksToPointCloud(Observation<VertexPointXYZCov>& pointVertices);

    /** get pose information from cloud */
    static Eigen::Affine3f getPose(PointCloudType& cloud);
    static Eigen::Affine3f getPose(PointCloudPtr& cloudPtr);

    /** create pose vertex for Affine3f pose */
    static g2o::VertexSE3* createVertexSE3(Eigen::Affine3f pose);

    /** create pose vertex for cloud with valid sensor origin and orientation */
    static g2o::VertexSE3* createVertexSE3(PointCloudType& cloud);
    static g2o::VertexSE3* createVertexSE3(PointCloudPtr& cloud);



    /** RGB pcl to uint conversion*/
    static void pclRGBToRGB(float& rgb, uint& r, uint& g, uint& b);
    static void pclRGBToRGB(float& rgb, unsigned char& r, unsigned char& g, unsigned char& b);
    static void RGBToPclRGB(float& rgb, uint& r, uint& g, uint& b);
    static void RGBToPclRGB(float& rgb, unsigned char& r, unsigned char& g, unsigned char& b);
    static void RGBToPclRGB(uint32_t& rgb, unsigned char& r, unsigned char& g, unsigned char& b);

    inline double getMaxDistance(){return _maxDistanceForNormalCalulation;}
    inline void setMaxDistance(double value){_maxDistanceForNormalCalulation = value;}

    void rgbToHsv(Eigen::Vector3i rgb, Eigen::Vector3i& hsv);

    protected:
    double _maxDistanceForNormalCalulation;
  };

  #include "pcl_ssa_hierarchical.hpp"
}
#endif

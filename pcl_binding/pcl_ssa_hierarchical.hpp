
template <typename PointCloudType>
PCLSSAHierarchicalT<PointCloudType>::PCLSSAHierarchicalT() : input_clouds_(0), sensor_(KINECT), levels_(1)
{
}

template <typename PointCloudType>
PCLSSAHierarchicalT<PointCloudType>::~PCLSSAHierarchicalT()
{

}

template <typename PointCloudType>
void
PCLSSAHierarchicalT<PointCloudType>::setInput (std::vector< PointCloudType >* clouds)
{
  input_clouds_  = clouds;
}

template <typename PointCloudType>
void 
PCLSSAHierarchicalT<PointCloudType>::setSensor (Sensor sensor)
{
  sensor_ = sensor;
}

template <typename PointCloudType>
void 
PCLSSAHierarchicalT<PointCloudType>::setLevels (int levels)
{
  levels_ = levels;
}

template <typename PointCloudType>
void 
PCLSSAHierarchicalT<PointCloudType>::setResolution (int level, double resolution)
{
  resolution_[level] = resolution;
}

template <typename PointCloudType>
void 
PCLSSAHierarchicalT<PointCloudType>::optimize (PointCloudType & result)
{
  std::cerr << __PRETTY_FUNCTION__ << ": not yet implemented" << std::endl;
  std::cerr << PVAR(result.size()) << std::endl;
}

template <typename PointCloudType>
  int 
  PCLSSAHierarchicalT<PointCloudType>::addCloudToGraph(SparseSurfaceAdjustmentGraph3D& graph, typename PointCloudType::Ptr cloud)
  {
    if(cloud->size() == 0)
      return 0;
    PointCloudType cleanCloud;
    double bottomLevelTime = get_time();
    ssa::RGBDSensorModel  sensorModel;

    //Create pose vertex
    Eigen::Affine3f pose = getPose(cloud);
    VertexSE3* vertex = createVertexSE3(pose); 
    int id = graph.addVertex(vertex);

    Observation<VertexPointXYZCov> observation;
    Vector3d measurement(0.0, 0.0, 0.0);
    for(unsigned int j=0; j < cloud->points.size(); ++j){
      typename PointCloudType::PointType& p = cloud->points[j];

      /** skip invalid readings */
      if(pcl_isnan(p.z))
        continue;

      measurement(0) = p.x;
      measurement(1) = p.y;
      measurement(2) = p.z;
      ///create point vertex
      VertexPointXYZCov* vp = new VertexPointXYZCov(0, measurement); /** id 0 means auto assign */
      vp->setParentVertex(vertex);
      vp->cr = p.r;
      vp->cg = p.g;
      vp->cb = p.b;
      #pragma omp critical
      {
        observation.push_back(vp);
        cleanCloud.push_back(p);
      }
    }
    std::cerr << "construction of level 0 took " << (get_time() - bottomLevelTime) * 1000 << "ms so far...  " << observation.size() << " points created. " << std::endl;
    observation.calcMeanCovThreaded(params_);
    std::cerr << "construction of level 0 took " << (get_time() - bottomLevelTime) * 1000 << "ms so far...  calculated mean and cov. " << std::endl;

    std::vector<VertexPointXYZCov* > vertices;
    vertices.reserve(observation.size());

    for(unsigned int j=0; j < observation.size(); ++j){
      VertexPointXYZCov*& vp = observation[j];

      if(vp->covariance() != Eigen::Matrix3d::Identity()){
        EdgeSE3PointXYZCov* e = new EdgeSE3PointXYZCov;
        e->vertices()[0] = vertex;
        e->vertices()[1] = vp; 
        Vector3d measurement = vp->estimate();
        e->setMeasurement(measurement);
        e->setLevel(0);
        sensorModel.getInformationMatrix(measurement, e->information());

        //move endpoint in global coordinate frame
        vp->setEstimate(vertex->estimate() * vp->estimate());
        #pragma omp critical
        {
          graph.addVertex(vp);
          graph.addEdge(e);
          vertices.push_back(vp);
        }
      } else {
        delete vp;
      }
    }
    observation.clear();

    std::cerr << "construction of level 0 at resolution " << resolution_[0] << " took " << (get_time() - bottomLevelTime) * 1000 << "ms." << std::endl;

    if(levels_ > 1 && vertices.size() > 0){
      KDTreeFlannT<VertexPointXYZCov>  kDTree;
      kDTree.copyData(vertices, true);
      kDTree.createKDTree();

      std::vector<int> k_indices;
      std::vector<float> k_squared_distances;

      for(int i=1; i < levels_; ++i){
        double timing = get_time();
        double resolution = resolution_[i];
        int levelCount = 0;
        pcl::VoxelGrid< typename PointCloudType::PointType > sor;
        sor.setInputCloud (boost::make_shared<PointCloudType >(cleanCloud));
        sor.setLeafSize (resolution, resolution, resolution); 
        sor.setDownsampleAllData(false);
        PointCloudType tmpCloud;
        sor.filter(tmpCloud);

	for(size_t j=0; j < tmpCloud.size(); ++j){

          Eigen::Vector3f point = pose * tmpCloud[j].getVector3fMap();
          kDTree.nearestKSearch(point, 1, k_indices, k_squared_distances);

          if(k_indices.size() > 0 && sqrt(k_squared_distances[0]) < resolution){
            VertexPointXYZCov*& vp = vertices[k_indices[0]];
            for(g2o::OptimizableGraph::EdgeSet::iterator it=vp->edges().begin(); it!=vp->edges().end(); it++){
              EdgeSE3PointXYZCov* e = dynamic_cast<EdgeSE3PointXYZCov* >(*it);
              if(e && e->vertices()[0]->id() == vertex->id()){
                e->setLevel(i);
                levelCount++;
                if(e->level() != i)
                  std::cerr << "setLevel failed " << e->level() << " != " <<  i << std::endl;
              } else {
                std::cerr << "something went wrong with casting ssa::EdgeSE3PointXYZCov!" << std::endl;
              }
            }
          }
        }
        std::cerr << "construction of level " << i << " at resolution " << resolution << " took " << (get_time() - timing) * 1000 << "ms. " << levelCount << std::endl;
      }
    }
    return id;
  }


  template <>
  pcl::PointCloud<pcl::PointXYZRGB>
  PCLSSAHierarchicalT< pcl::PointCloud<pcl::PointXYZRGB> >::landmarksToPointCloud(std::vector<EdgeSE3PointXYZCov* >& landmarksFromVertex)
  {
    pcl::PointCloud<pcl::PointXYZRGB> cloud;
    for(unsigned int k = 0; k < landmarksFromVertex.size(); k++)
    {
      EdgeSE3PointXYZCov*& edge = landmarksFromVertex[k];
      VertexSE3* pose  = dynamic_cast<VertexSE3* >(edge->vertices()[0]);
      VertexPointXYZCov* point = dynamic_cast<VertexPointXYZCov* >(edge->vertices()[1]);
      pcl::PointXYZRGB p; 

      if(point->covariance() == Eigen::Matrix3d::Identity())
        continue;

      p.x = (pose->estimate() * edge->measurement())[0];
      p.y = (pose->estimate() * edge->measurement())[1];
      p.z = (pose->estimate() * edge->measurement())[2];

//       p.x = point->estimate()[0];
//       p.y = point->estimate()[1];
//       p.z = point->estimate()[2];

      RGBToPclRGB(p.rgb, point->cr, point->cg, point->cb);
      cloud.push_back(p);
    }
    return cloud;
  }

  template <>
  pcl::PointCloud<pcl::PointXYZRGBA>
  PCLSSAHierarchicalT< pcl::PointCloud<pcl::PointXYZRGBA> >::landmarksToPointCloud(std::vector<EdgeSE3PointXYZCov* >& landmarksFromVertex)
  {
    pcl::PointCloud<pcl::PointXYZRGBA> cloud;
    for(unsigned int k = 0; k < landmarksFromVertex.size(); k++)
    {
      EdgeSE3PointXYZCov*& edge = landmarksFromVertex[k];
      VertexSE3* pose  = dynamic_cast<VertexSE3* >(edge->vertices()[0]);
      VertexPointXYZCov* point = dynamic_cast<VertexPointXYZCov* >(edge->vertices()[1]);
      pcl::PointXYZRGBA p; 

      if(point->covariance() == Eigen::Matrix3d::Identity())
        continue;

      p.x = (pose->estimate() * edge->measurement())[0];
      p.y = (pose->estimate() * edge->measurement())[1];
      p.z = (pose->estimate() * edge->measurement())[2];

      RGBToPclRGB(p.rgba, point->cr, point->cg, point->cb);
      cloud.push_back(p);
    }
    return cloud;
  }

  template <>
  pcl::PointCloud<pcl::PointXYZRGB> 
  PCLSSAHierarchicalT< pcl::PointCloud<pcl::PointXYZRGB> >::landmarksToPointCloud(Observation<VertexPointXYZCov>& pointVertices){
    pcl::PointCloud<pcl::PointXYZRGB> cloud;
    for(unsigned int k = 0; k < pointVertices.size(); k++)
    {
      VertexPointXYZCov* point = pointVertices[k];
      if(point->covariance() == Eigen::Matrix3d::Identity())
        continue;
      pcl::PointXYZRGB p; 

      p.x = point->estimate()[0];
      p.y = point->estimate()[1];
      p.z = point->estimate()[2];

      RGBToPclRGB(p.rgb, point->cr, point->cg, point->cb);
      cloud.push_back(p);
    }
    return cloud;
  }

  template <>
  pcl::PointCloud<pcl::PointXYZRGBA> 
  PCLSSAHierarchicalT< pcl::PointCloud<pcl::PointXYZRGBA> >::landmarksToPointCloud(Observation<VertexPointXYZCov>& pointVertices){
    pcl::PointCloud<pcl::PointXYZRGBA> cloud;
    for(unsigned int k = 0; k < pointVertices.size(); k++)
    {
      VertexPointXYZCov* point = pointVertices[k];
      if(point->covariance() == Eigen::Matrix3d::Identity())
        continue;
      pcl::PointXYZRGBA p; 

      p.x = point->estimate()[0];
      p.y = point->estimate()[1];
      p.z = point->estimate()[2];

      RGBToPclRGB(p.rgba, point->cr, point->cg, point->cb);
      if(p.r != point->cr)
        std::cerr << "ERROR red channel "  << p.r << " != " << point->cr << std::endl;
      if(p.g != point->cg)
        std::cerr << "ERROR green channel "  << p.g << " != " << point->cg << std::endl;
      if(p.b != point->cb)
        std::cerr << "ERROR blue channel "  << p.b << " != " << point->cb << std::endl;
      cloud.push_back(p);
    }
    return cloud;
  }

  template <>
  pcl::PointCloud<pcl::PointXYZRGBNormal> 
  PCLSSAHierarchicalT< pcl::PointCloud<pcl::PointXYZRGBNormal> >::landmarksToPointCloud(std::vector<EdgeSE3PointXYZCov* >& landmarksFromVertex)
  {

    pcl::PointCloud<pcl::PointXYZRGBNormal> cloud;
    for(unsigned int k = 0; k < landmarksFromVertex.size(); k++)
    {
      EdgeSE3PointXYZCov*& edge = landmarksFromVertex[k];
      VertexSE3* pose  = dynamic_cast<VertexSE3* >(edge->vertices()[0]);
      VertexPointXYZCov* point = dynamic_cast<VertexPointXYZCov* >(edge->vertices()[1]);
      pcl::PointXYZRGBNormal p; 

      if(point->covariance() == Eigen::Matrix3d::Identity())
        continue;

      p.x = (pose->estimate() * edge->measurement())[0];
      p.y = (pose->estimate() * edge->measurement())[1];
      p.z = (pose->estimate() * edge->measurement())[2];

//       p.x = point->estimate()[0];
//       p.y = point->estimate()[1];
//       p.z = point->estimate()[2];

      RGBToPclRGB(p.rgb, point->cr, point->cg, point->cb);
      cloud.push_back(p);
    }
    return cloud;
  }

  template <>
  pcl::PointCloud<pcl::PointXYZRGBNormal> 
  PCLSSAHierarchicalT< pcl::PointCloud<pcl::PointXYZRGBNormal> >::landmarksToPointCloud(Observation<VertexPointXYZCov>& pointVertices){
    pcl::PointCloud<pcl::PointXYZRGBNormal> cloud;
    for(unsigned int k = 0; k < pointVertices.size(); k++)
    {
      VertexPointXYZCov* point = pointVertices[k];
      if(point->covariance() == Eigen::Matrix3d::Identity())
        continue;
      pcl::PointXYZRGBNormal p; 

      p.x = point->estimate()[0];
      p.y = point->estimate()[1];
      p.z = point->estimate()[2];

      p.normal_x = point->normal().x();
      p.normal_y = point->normal().y();
      p.normal_z = point->normal().z();


      RGBToPclRGB(p.rgb, point->cr, point->cg, point->cb);
      cloud.push_back(p);
    }
    return cloud;
  }

  template <typename PointCloudType>
  void 
  PCLSSAHierarchicalT<PointCloudType>::pclRGBToRGB(float& rgb, uint& r, uint& g, uint& b)
  {
    //convert color information
    int irgb = *reinterpret_cast<int*>(&rgb);
    r = ((irgb >> 16) & 0xff);
    g = ((irgb >> 8) & 0xff);
    b = (irgb & 0xff);
  }

  template <typename PointCloudType>
  void 
  PCLSSAHierarchicalT<PointCloudType>::pclRGBToRGB(float& rgb, unsigned char& r, unsigned char& g, unsigned char& b)
  {
    //convert color information
    int irgb = *reinterpret_cast<int*>(&rgb);
    r = ((irgb >> 16) & 0xff);
    g = ((irgb >> 8) & 0xff);
    b = (irgb & 0xff);
  }

  template <typename PointCloudType>
  void 
  PCLSSAHierarchicalT<PointCloudType>::RGBToPclRGB(float& rgb, uint& r, uint& g, uint& b)
  {
    //convert color information
    unsigned char* p = (unsigned char*) &rgb;
    p[0] = (unsigned char) b;
    p[1] = (unsigned char) g;
    p[2] = (unsigned char) r;
    p[3] = 0;
  }

  template <typename PointCloudType>
  void 
  PCLSSAHierarchicalT<PointCloudType>::RGBToPclRGB(float& rgb, unsigned char& r, unsigned char& g, unsigned char& b)
  {
    //convert color information
    unsigned char* p = (unsigned char*) &rgb;
    p[0] = b;
    p[1] = g;
    p[2] = r;
    p[3] = 0;
  }

  template <typename PointCloudType>
  void 
  PCLSSAHierarchicalT<PointCloudType>::RGBToPclRGB(uint32_t& rgb, unsigned char& r, unsigned char& g, unsigned char& b)
  {
    //convert color information
    unsigned char* p = (unsigned char*) &rgb;
    p[0] = b;
    p[1] = g;
    p[2] = r;
    p[3] = 0;
  }

  template <typename PointCloudType>
  Eigen::Affine3f 
  PCLSSAHierarchicalT<PointCloudType>::getPose(PointCloudType& cloud){
    Eigen::Affine3f pose = (Eigen::Affine3f)Eigen::Translation3f(cloud.sensor_origin_(0), cloud.sensor_origin_(1), cloud.sensor_origin_(2)) * cloud.sensor_orientation_;
    return pose;
  }

  template <typename PointCloudType>
  Eigen::Affine3f 
  PCLSSAHierarchicalT<PointCloudType>::getPose(PointCloudPtr& cloudPtr){
    Eigen::Affine3f pose = (Eigen::Affine3f)Eigen::Translation3f(cloudPtr->sensor_origin_(0), cloudPtr->sensor_origin_(1), cloudPtr->sensor_origin_(2)) * cloudPtr->sensor_orientation_;
    return pose;
  }


  template <typename PointCloudType>
  g2o::VertexSE3* 
  PCLSSAHierarchicalT<PointCloudType>::createVertexSE3(Eigen::Affine3f pose){
    ///Create pose vertex
    VertexSE3* vertex = new VertexSE3(); 
    vertex->setId(0);
    vertex->setEstimate(Affine3fToSE3Quat(pose));
    return vertex;
  }

  template <typename PointCloudType>
  g2o::VertexSE3* 
  PCLSSAHierarchicalT<PointCloudType>::createVertexSE3(PointCloudType& cloud){
    return createVertexSE3(getPose(cloud));
  }


  template <typename PointCloudType>
  g2o::VertexSE3* 
  PCLSSAHierarchicalT<PointCloudType>::createVertexSE3(PointCloudPtr& cloudPtr){
    return createVertexSE3(getPose(cloudPtr));
  }


  template <typename PointCloudType>
  void PCLSSAHierarchicalT<PointCloudType>::rgbToHsv(Eigen::Vector3i rgb, Eigen::Vector3i& hsv)
  {
    float min;
     hsv(2) = std::max (rgb(0), std::max (rgb(1), rgb(2)));
     min = std::min (rgb(0), std::min (rgb(1), rgb(2)));

     if (hsv(2) != 0)
       hsv(1) = (hsv(2) - min) / hsv(2);
     else
     {
       hsv(1) = 0;
       hsv(0) = -1;
       return;
     }
 
     if (rgb(0) == hsv(2))
       hsv(0) = (rgb(1) - rgb(2)) / (hsv(2) - min);
     else if (rgb(1) == hsv(2))
       hsv(0) = 2 + (rgb(2) - rgb(0)) / (hsv(2) - min);
     else 
       hsv(0) = 4 + (rgb(0) - rgb(1)) / (hsv(2) - min);
     hsv(0) *= 60;
     if (hsv(0) < 0)
       hsv(0) += 360;
  }
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

#ifndef __SPARSE_SURFACE_ADJUSTMENT_GRAPH__
#define __SPARSE_SURFACE_ADJUSTMENT_GRAPH__

//#define NUMERIC_JACOBIAN_TWO_D_TYPES

#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "observation.h"
#include "ssa/data_association/correspondence.h"

#include "ssa_sparse_optimizer.h"

#include "local_neighborhood.h"

#include "g2o/core/optimizable_graph.h"
#include "g2o/core/factory.h"
#include "g2o/stuff/macros.h"
#include "g2o/stuff/timeutil.h"
#include "g2o/types/icp/types_icp.h"
//#include "g2o/apps/g2o_hierarchical/edge_labeler.h"

namespace ssa {

  //forward declaration
  class SSASparseOptimizer;

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  class SparseSurfaceAdjustmentGraphT{

    public:
    typedef EdgeType1 SLAMEdgeType;
    typedef EdgeType2 SensorEdgeType;
    typedef EdgeType3 DataAssociationEdgeType;

    typedef typename SensorEdgeType::VertexXiType   PoseVertex;
    static const int Di = SensorEdgeType::Dimension;
    typedef Eigen::Matrix<double, Di, Di>  PoseMatrix;

    typedef typename SensorEdgeType::VertexXjType   PointVertex;
    static const int Dj = SensorEdgeType::Dimension;
    typedef Eigen::Matrix<double, Dj, 1> PointVector;
    typedef Eigen::Matrix<double, Dj, Dj> PointMatrix;

    //typedef KDTree::KDTree<Dj, PointVertex*, KdTreeAccess>    PointTree;


    SparseSurfaceAdjustmentGraphT();
    ~SparseSurfaceAdjustmentGraphT();
  
    void load(std::string filename);
    void save(std::string filename);

    void clear();

    /** This method sets the parent vertex pointer of all 
        PointVertex vertices and fills the pointer into 
        type specific vectors for faster and easier access. */
    void linkNodesToVertices();
    void fillObservations();
  
    //!returns a map with the ids of all scans
    std::map<int, int> getScanIds();
    std::vector<int>   getPoseIds();
  
    //--------------------------
    // Mean and covariance stuff
    //--------------------------
    void calcMeanCov(SparseSurfaceAdjustmentParams& params);
    //! calculate mean and covariance of for scan with the given scanId
    void calcMeanCov(SparseSurfaceAdjustmentParams& params, int scandId);

    /** Fill local observation neighbor cache */
    void fillNeighborCache(SparseSurfaceAdjustmentParams& params);

    /** Fill local observation neighbor cache of given scan */
    void fillNeighborCache(SparseSurfaceAdjustmentParams& params, int scanId);

    /** Access to the local observation neighbor cache */
    inline LocalNeighborhoodT<PointVertex>& neighborCache(){ return cached_neighbors_;};

    /** get the id of the highest level in the hierarchcal graph */
    inline int getMaxLevel(){return max_level_;}

    /** Returns a subset of edges from the current graph for the 
      optimization. The subset contains all edge types for point 
      vertices  that have a data association to another point 
      vertex from a different observation.*/
    g2o::OptimizableGraph::EdgeSet getEdgeset();
    g2o::OptimizableGraph::EdgeSet getEdgesetFast();
    g2o::OptimizableGraph::EdgeSet getEdgesetMesh();
  
    /** move points that have no data association */
    void moveUnoptimizedPoints();

    /** hide not well connected points*/
    void filterOutlier(int minNeighbors);

    //!returns a vector with all observation vertices without the given scan observations
    std::vector<PointVertex* > getPointVerticesOfScans(int& scanId);

    std::vector<PointVertex* > getPointVerticesOfScan(int& scanId);
    Observation<PointVertex>&  getObservationOfScan(int& scanId);
  
    int addVertex(PoseVertex*& v);
    int addVertex(PointVertex*& v);
    int addVertex(PointVertex*& v, std::vector< PointVertex* >& neighbors);

    void removeVertex(PoseVertex*& v);
    void removeVertex(PointVertex*& v);

    void addEdge(SLAMEdgeType*& edge);
    void addEdge(SensorEdgeType*& edge);
    void addEdge(PointVertex*& from, PointVertex*& to);
    void addEdge(DataAssociationEdgeType*& edge);
    void addEdge(g2o::Edge_V_V_GICP*& edge);
    DataAssociationEdgeType* createEdge(PointVertex*& from, PointVertex*& to);

    void removeEdge(SLAMEdgeType* edge);
    void removeEdge(SensorEdgeType* edge);
    void removeEdge(DataAssociationEdgeType* edge);
  
    void dropDataAssociation();
    void fixDataAssociation();
    void dropSensorEdges();
    void dropSLAMEdges();

    void unfixVertices();
    void setFixedVertices(bool fixed);

    void pruneGraph(double minDistance);
    void pruneDataAssociation(int maxEdges);
    void marginalizeEdges(PoseVertex* v);

    bool connected(PoseVertex*& v1, PoseVertex*& v2);
    bool connected(PointVertex*& v1, PointVertex*& v2);

    bool haveCommonSensorVertex(PointVertex*& v1, PointVertex*& v2); 

    void checkGraph();

    Observation<PointVertex>&  observation(int i);

    std::vector<PointVertex* >&  getPointVertices(int poseId, int level);

    /** Applies some heuristics for a sparser (faster) connectivity */
    void checkConnectivity();

    bool DataAssociationEdgeComp (DataAssociationEdgeType* i,DataAssociationEdgeType* j) { return (i->chi2()<j->chi2()); }

    /** Optimizer */
    SSASparseOptimizer _optimizer;
    SSASparseOptimizer*   optimizer() { return &_optimizer;};

    /** Direct access to differnt edgetypes */
  
    std::vector<PoseVertex* >                   _verticies_poses;
    std::vector<PointVertex* >                  _verticies_points;

    /** Caching of point vertices hierarchcal structure. 
        Will be filled during first call of getPointVertices(poseIdm level).
        First index is pose vertex id 
        second index is the level.*/
    std::map< int, std::map< int, std::vector<PointVertex* > > >  _hierarchical_point_vertex_index;
    
    /** Caching of correspondences */
    std::vector<CorrespondenceT<PointVertex, EdgeType3> >correspondences;
    //deprecated 
    std::map< int, Observation<PointVertex> >   _verticies_observations;

    std::vector< SLAMEdgeType* >                _edges_odometry;
    std::vector< SensorEdgeType* >              _edges_observations;
    std::vector<DataAssociationEdgeType* >      _edges_data_association;
    std::vector<g2o::Edge_V_V_GICP* >           _edges_data_association_gicp;
    std::vector<DataAssociationEdgeType* >      _edges_surface_mesh;
  
    std::map<int, std::map<int, bool> >         _fixedDataAssociation;
    int maxVertexId;
    
    /** vertices/edges of points that are not in the current optimization set 
        are moved accordingly to the resulting observation position */
    std::vector<SensorEdgeType* >               _edges_points_to_move;
    std::vector<PointVector, Eigen::aligned_allocator<PointVector> >   _edges_points_to_move_transforms;


    private:
    /** caching local neighborhood */
    LocalNeighborhoodT<PointVertex>             cached_neighbors_;
    int                                         max_level_;
  };

} //end namespace

#include "ssa_graph.hpp"

#endif
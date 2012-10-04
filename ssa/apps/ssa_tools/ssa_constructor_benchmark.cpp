// Sparse Surface Optimization
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

#include <iostream>
#include <signal.h>

#include "ssa/core/ssa_graph_3d.h"

#include "ssa/core/allocate_solver.h"
#include "ssa/core/sparse_surface_adjustment.h"
#include "g2o/types/sba/types_sba.h"

using namespace std;
using namespace ssa;
using namespace g2o;

const char *message[]={
  "ssa_constructor_benchmark: constructs n vertex objects of different types",
  "usage ssa_constructor_benchmark [options]",
  "options:",
  "-n     number of objects to create (default: 1000000)",
  0
};

int main(int argc, char **argv)
{

  const char**v=message;
  while (*v){
    cout << *v << endl;
    v++;
  }

  unsigned int n = 1000000;
  int c=1;
  while (c<argc){
    if (!strcmp(argv[c],"-n")){
      c++;
      n = (unsigned int) atoi(argv[c]);
      c++;
    }
  }
  // HACK reset numeric locale, we need C
  setlocale (LC_NUMERIC,"C");

//   cerr << "Benchmarking HyperGraphElement" << endl;
//   std::vector<HyperGraph::HyperGraphElement* > hyperGraphElements;
//   hyperGraphElements.reserve(n);
  double timing = get_time();
//   for(unsigned int i = 0; i < n; ++i)
//   {
//     HyperGraph::HyperGraphElement* v = new HyperGraph::HyperGraphElement();
//     hyperGraphElements.push_back(v);
//   }
//   cerr << "creation of " << n << " HyperGraph::HyperGraphElement objects took " << (get_time() - timing) * 1000 << "  ms." << endl;
// 
//   timing = get_time();
//   for(unsigned int i = 0; i < n; ++i)
//   {
//     HyperGraph::HyperGraphElement*& v = hyperGraphElements[i];
//     delete v;
//   }
//   cerr << "deletion of " << n << " HyperGraph::HyperGraphElement objects took " << (get_time() - timing) * 1000 << "  ms." << endl;


  cerr << "Benchmarking HyperGraph::Vertex" << endl;
  std::vector<HyperGraph::Vertex* > hyperGraphVertices;
  hyperGraphVertices.reserve(n);
  timing = get_time();
  for(unsigned int i = 0; i < n; ++i)
  {
    HyperGraph::Vertex* v = new HyperGraph::Vertex;
    hyperGraphVertices.push_back(v);
  }
  cerr << "creation of " << n << " HyperGraph::Vertex objects took " << (get_time() - timing) * 1000 << "  ms." << endl;

  timing = get_time();
  for(unsigned int i = 0; i < n; ++i)
  {
    HyperGraph::Vertex*& v = hyperGraphVertices[i];
    delete v;
  }
  cerr << "deletion of " << n << " HyperGraph::Vertex objects took " << (get_time() - timing) * 1000 << "  ms." << endl;

  cerr << "Benchmarking VertexSE3" << endl;
  std::vector<VertexSE3* > poseVertices;
  poseVertices.reserve(n);
  timing = get_time();
  for(unsigned int i = 0; i < n; ++i)
  {
    VertexSE3* v = new VertexSE3;
    poseVertices.push_back(v);
  }
  cerr << "creation of " << n << " VertexSE3 objects took " << (get_time() - timing) * 1000 << "  ms." << endl;

  timing = get_time();
  for(unsigned int i = 0; i < n; ++i)
  {
    VertexSE3*& v = poseVertices[i];
    delete v;
  }
  cerr << "deletion of " << n << " VertexSE3 objects took " << (get_time() - timing) * 1000 << "  ms." << endl;

  cerr << "Benchmarking VertexPointXYZCov" << endl;
  std::vector<VertexPointXYZCov* > baseVertices;
  baseVertices.reserve(n);
  timing = get_time();
  for(unsigned int i = 0; i < n; ++i)
  {
    VertexPointXYZCov* v = new VertexPointXYZCov;
    baseVertices.push_back(v);
  }
  cerr << "creation of " << n << " ssa::VertexPointXYZCov objects took " << (get_time() - timing) * 1000 << "  ms." << endl;

  timing = get_time();
  for(unsigned int i = 0; i < n; ++i)
  {
    VertexPointXYZCov*& v = baseVertices[i];
    delete v;
  }
  cerr << "deletion of " << n << " ssa::VertexPointXYZCov objects took " << (get_time() - timing) * 1000 << "  ms." << endl;
  
//   cerr << "Benchmarking VertexPointXYZCov Threaded " << endl;
//   omp_set_num_threads(2);
//   timing = get_time();
//   std::vector<VertexPointXYZCov* > baseVertices2;
//   baseVertices2.resize(n);
//   
//   #pragma omp parallel for shared(baseVertices2) 
//   for(unsigned int i = 0; i < n; ++i)
//   {
//     baseVertices2[i] = new VertexPointXYZCov;
//   }
//   cerr << "creation of " << n << " ssa::VertexPointXYZCov objects took " << (get_time() - timing) * 1000 << "  ms." << endl;
// 
//   timing = get_time();
//   for(unsigned int i = 0; i < n; ++i)
//   {
//     VertexPointXYZCov*& v = baseVertices2[i];
//     delete v;
//   }
//   cerr << "deletion of " << n << " ssa::VertexPointXYZCov objects took " << (get_time() - timing) * 1000 << "  ms." << endl;
  

  
  cerr << "Benchmarking VertexPointXYZCov vector" << endl;
  timing = get_time();
  std::vector<VertexPointXYZCov> pointVerticesVector;
  pointVerticesVector.resize(n);
  cerr << "vector creation of " << n << " g2o::VertexPointXYZCov objects took " << (get_time() - timing) * 1000 << "  ms." << endl;

  timing = get_time();
  {
    std::vector<VertexPointXYZCov> tmp;
    std::swap(pointVerticesVector, tmp);
  }
  cerr << "vector deletion of " << n << " g2o::VertexPointXYZCov objects took " << (get_time() - timing) * 1000 << "  ms." << endl;
  
  cerr << "Benchmarking VertexPointXYZ" << endl;
  std::vector<VertexPointXYZ* > pointVertices;
  pointVertices.reserve(n);
  timing = get_time();
  for(unsigned int i = 0; i < n; ++i)
  {
    VertexPointXYZ* v = new VertexPointXYZ;
    pointVertices.push_back(v);
  }
  cerr << "creation of " << n << " g2o::VertexPointXYZ objects took " << (get_time() - timing) * 1000 << "  ms." << endl;

  timing = get_time();
  for(unsigned int i = 0; i < n; ++i)
  {
    VertexPointXYZ*& v = pointVertices[i];
    delete v;
  }
  cerr << "deletion of " << n << " g2o::VertexPointXYZ objects took " << (get_time() - timing) * 1000 << "  ms." << endl;
  
  cerr << "Benchmarking Eigen structs" << endl;
  std::vector<Eigen::Matrix3d* > matrices;
  matrices.reserve(n);
  timing = get_time();
  for(unsigned int i = 0; i < n; ++i)
  {
    Eigen::Matrix3d* v = new Eigen::Matrix3d;
    matrices.push_back(v);
  }
  cerr << "creation of " << n << " Eigen::Matrix3d objects took " << (get_time() - timing) * 1000 << "  ms." << endl;

  timing = get_time();
  for(unsigned int i = 0; i < n; ++i)
  {
    Eigen::Matrix3d*& v = matrices[i];
    delete v;
  }
  cerr << "deletion of " << n << " Eigen::Matrix3d objects took " << (get_time() - timing) * 1000 << "  ms." << endl;

  cerr << "Benchmarking stl containers" << endl;
// <EstimateType,  Eigen::aigned_allocator<EstimateType> >
  std::vector<std::deque<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> >* > deques;
  deques.reserve(n);
  timing = get_time();
  for(unsigned int i = 0; i < n; ++i)
  {
    std::deque<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> >* v = new std::deque<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> >;
    deques.push_back(v);
  }
  cerr << "creation of " << n << " std::deque<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > objects took " << (get_time() - timing) * 1000 << "  ms." << endl;

  timing = get_time();
  for(unsigned int i = 0; i < n; ++i)
  {
    std::deque<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> >*& v = deques[i];
    delete v;
  }
  cerr << "deletion of " << n << " std::deque<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > objects took " << (get_time() - timing) * 1000 << "  ms." << endl;
}
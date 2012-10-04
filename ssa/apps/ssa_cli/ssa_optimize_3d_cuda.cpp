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

#include <iostream>
#include <signal.h>

#include "ssa/core/ssa_graph_3d.h"

#include "ssa/core/allocate_solver.h"
#include "ssa/core/sparse_surface_adjustment.h"
#include "ssa/solvers/pcg_cuda/linear_solver_pcg_cuda.h"
//#include "ssa/solvers/pcg_cusp/linear_solver_pcg_cusp.h"

using namespace std;
using namespace ssa;

const char *message[]={
  "ssa_optimize: optimizes a given ssa file",
  "usage ssa_optimize [options] <ssa_file>",
  "options:",
  "-ini filename  provide parameter set in an ini file.",
  "-o [filename]  save result in given file",
  0
};

static bool running = true;
unsigned int breakCounter = 0;

void sighandler(int sig)
{
    cout << endl << "Signal " << sig << " caught..." << endl;

   running = false;
   breakCounter++;
   if(breakCounter >= 3){
      exit(1);
   }
}


int main(int argc, char **argv)
{

  const char**v=message;
  while (*v){
    cout << *v << endl;
    v++;
  }

  // HACK reset numeric locale, we need C
  setlocale (LC_NUMERIC,"C");
  int c=1;
  const char* logfile=0;
  const char* outfile=0;
  bool saveOutput = false;

  bool useConfigFile = false;
  const char* configFile=0;
 
  while (c<argc){
    if (!strcmp(argv[c],"-o")){
      saveOutput = true;
      c++;
      outfile=argv[c];
      c++;
    } else
     if (!strcmp(argv[c],"-ini")){
      useConfigFile=true;
      c++;
      configFile = argv[c];
      c++;
    } else 
    if (! logfile){
      logfile=argv[c];
      c++;
      break;
    }
  }

  if(!logfile){
    cerr << "Error: no input ssa file given." << endl;
    exit(0);
  }
  if(!outfile){
    cerr << "Warning: missing output file." << endl;
  }

  //initialize optimizer
  SparseSurfaceAdjustmentT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov> ssa;
  ssa.setVerbose(true);
  //g2o::BlockSolverX::LinearSolverType* linearSolver = new g2o::LinearSolverPCG<g2o::BlockSolverX::PoseMatrixType>;
  g2o::BlockSolverX::LinearSolverType* linearSolver = new g2o::LinearSolverCUDACG<g2o::BlockSolverX::PoseMatrixType>;
  //g2o::BlockSolverX::LinearSolverType* linearSolver = new g2o::LinearSolverCuspCG<g2o::BlockSolverX::PoseMatrixType>;

  ssa.setSolver(linearSolver);

  if (! logfile){

  } else {
    ssa.graph()->load(logfile);
  }

  if(useConfigFile){
    ifstream configStream(configFile);
    ssa.params().readParams(configStream);
    configStream.close();
    ssa.params().printParams();
  }
  //omp_set_num_threads(1);
  double startTime = get_time();
  ssa.optimize(); 
  cerr << "calculation took " << get_time() - startTime << " seconds." << endl;

  if(saveOutput){
    cerr << "wrote resulting graph to " << outfile << endl;
    ssa.graph()->save(outfile);
  }

  return 0;
}

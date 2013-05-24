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

#include <iostream>
#include <signal.h>
#include <QApplication>
#include <QtGui/QMainWindow>
#include <boost/thread.hpp>
#include "ui_base_main_window.h"

#include "ssa/core/allocate_solver.h"
#include "ssa/core/ssa_graph_2d.h"
#include "ssa/core/sparse_surface_adjustment.h"
#include "ssa/data_association/local_refinement.h"
//#include "ssa/solvers/pcg_cuda/linear_solver_pcg_cuda.h"
// #include "ssa/solvers/pcg_cusp/linear_solver_pcg_cusp.h"

using namespace std;
using namespace ssa;

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
  signal(SIGINT, &sighandler);
  QApplication qapp(argc, argv);
  // HACK reset numeric locale, we need C
  setlocale (LC_NUMERIC,"C");
  //bool debug = false;
  bool restoreViewerState = false;

  const char* logfile=0;
  const char* viewerStateFile=0;
  bool useConfigFile = false;
  const char* configFile=0;
  int c=1;
  while (c<argc){
    if (!strcmp(argv[c],"-cam")){
      restoreViewerState=true;
      c++;
      viewerStateFile = argv[c];
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

  //initialize optimizer
  SparseSurfaceAdjustmentT<g2o::EdgeSE2, ssa::EdgeSE2PointXYCov, ssa::EdgePointXYCovPointXYCov> ssa;
  ssa.setVerbose(true);
  g2o::BlockSolverX::LinearSolverType* linearSolver = AllocateLinearSolver<g2o::BlockSolverX>(1); //CHOLMOD
  ssa.setSolver(linearSolver);

  if(useConfigFile){
    ifstream configStream(configFile);
    ssa.params().readParams(configStream);
    configStream.close();
    ssa.params().printParams();
  }

  if (! logfile){

  } else {
    ssa.graph()->load(logfile);
  }



  QMainWindow* mw = new QMainWindow;
  Ui_MainWindow umw;
  umw.setupUi(mw);

  umw.ssaGLWidget->setMap(ssa.graph());
  mw->show();
  umw.ssaGLWidget->setBackgroundColor(qRgb(211, 211, 211));
  if(restoreViewerState){
    umw.ssaGLWidget->setStateFileName(viewerStateFile);
    umw.ssaGLWidget->restoreStateFromFile();
  }

  umw.ssaGLWidget->_iterations = ssa.params().ssaIterations;
  umw.doubleSpinBox_3->setValue(ssa.params().normalShooting.stepSize);
  umw.doubleSpinBox_4->setValue(ssa.params().normalShooting.steps);
  umw.doubleSpinBox_5->setValue(ssa.params().normalShooting.maxAngleDifference);

  QObject::connect( &(ssa.graph()->_optimizer), SIGNAL(iterationDone()), umw.ssaGLWidget, SLOT(Gen3DObjectList_update()));
  QObject::connect( &(ssa.graph()->_optimizer), SIGNAL(iterationDone()), umw.ssaGLWidget, SLOT(updateGL()));

  if(ssa.graph()->_verticies_poses.size() > 0){
    g2o::VertexSE2* v = ssa.graph()->_verticies_poses[0];
    qglviewer::Vec initialNodePose = qglviewer::Vec(v->estimate().translation()(0), v->estimate().translation()(1), 20.0);
    umw.ssaGLWidget->camera()->setPosition(initialNodePose);
    initialNodePose = qglviewer::Vec(v->estimate().translation()(0), v->estimate().translation()(1), 0.0f);
    umw.ssaGLWidget->camera()->setSceneCenter(initialNodePose);
    umw.ssaGLWidget->camera()->setUpVector(qglviewer::Vec(0.0, 0.0, 1.0));
    umw.ssaGLWidget->camera()->lookAt(initialNodePose);
  }
  cerr << ssa.graph()->_verticies_poses.size() << endl;
  ssa.graph()->setFixedVertices(false);

  boost::mutex mutex;
  bool updateGL = true;
  umw.ssaGLWidget->show();
  qapp.processEvents();
  ssa.graph()->fillNeighborCache(ssa.params());
  omp_set_nested(1);
  QPalette palette;
  #pragma omp parallel
  {
    ///Gui thread
    if(omp_get_thread_num() == 0){
      while(mw->isVisible() && running) {
	qapp.processEvents();
	usleep(1000);

	if(updateGL){
	  mutex.lock();
	  updateGL = false;
	  umw.ssaGLWidget->Gen3DObjectList_update();
	  umw.ssaGLWidget->updateGL();
	  qapp.processEvents();
	  mutex.unlock();
	}
	if(umw.ssaGLWidget->_optimize){
	  umw.stateLabel->setText(QString("running"));
	  palette.setColor(umw.stateLabel->backgroundRole(), Qt::yellow);
	  umw.stateLabel->setPalette(palette);
	  umw.stateLabel->setAutoFillBackground(true);
	  mw->repaint();
	} else {
	  umw.stateLabel->setText(QString("ready"));
	  palette.setColor(umw.stateLabel->backgroundRole(), Qt::green);
	  umw.stateLabel->setPalette(palette);
	  umw.stateLabel->setAutoFillBackground(true);
	}
      }
    }
    if(omp_get_thread_num() == 1){
      while(mw->isVisible() && running) {
	usleep(1000);
	if(umw.ssaGLWidget->_optimize){
	  ///filling selected gui params into ssa params
	  ssa.params().ssaIterations = umw.ssaGLWidget->_iterations;
	  //ssa.params().g2oIterations = 6;
	  ssa.params().normalShooting.stepSize = umw.doubleSpinBox_3->value();
	  ssa.params().normalShooting.steps = (int) umw.doubleSpinBox_4->value();
	  ssa.params().normalShooting.maxAngleDifference = umw.doubleSpinBox_5->value();
	  ssa.params().printParams();
	  cerr << "optimizing level " << umw.spinBox_2->value() << endl;
	  ssa.optimize(umw.spinBox_2->value());
	  //ssa.optimize(ssa.graph()->getMaxLevel());

	  mutex.lock();
	    updateGL = true;
	  mutex.unlock();
	  umw.ssaGLWidget->_optimize = false;
	}
      }
    }
  }


  return 0;
}

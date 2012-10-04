INCLUDE(${QT_USE_FILE})

# ${QT_INCLUDE_DIR} does not appear to work; therefore have to
# manually rig it; also have to add QtCore directly because qglviewer
# doesn't seem to handle the Qt directories properly

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CHOLMOD_INCLUDE_DIR} ${CSPARSE_INCLUDE_DIR} ${QT_QTCORE_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${QT_QTCORE_INCLUDE_DIR}/..)
 
ADD_EXECUTABLE(ssa_optimize_3d
  ssa_optimize_3d.cpp
)

TARGET_LINK_LIBRARIES(ssa_optimize_3d ssa_core types_3d g2o_stuff g2o_types_slam3d g2o_types_icp g2o_core g2o_solver_cholmod ${CHOLMOD_LIBRARIES} ${BLAS_LIBRARIES} g2o_csparse_extension g2o_solver_csparse flann ${CSPARSE_LIBRARY} ${QT_QTXML_LIBRARY} ${QT_QTCORE_LIBRARY})

ADD_EXECUTABLE(ssa_optimize_hierarchical_3d
  ssa_optimize_hierarchical_3d.cpp
)

TARGET_LINK_LIBRARIES(ssa_optimize_hierarchical_3d ssa_core types_3d g2o_stuff g2o_types_slam3d g2o_types_icp g2o_core g2o_solver_cholmod ${CHOLMOD_LIBRARIES} ${BLAS_LIBRARIES} g2o_solver_csparse g2o_csparse_extension flann ${CSPARSE_LIBRARY} ${QT_QTXML_LIBRARY} ${QT_QTCORE_LIBRARY})

IF(CUDA_FOUND AND SSA_USE_CUDA)

INCLUDE_DIRECTORIES(${CUDA_SDK_ROOT_DIR}/shared/inc/ ${CUDA_SDK_ROOT_DIR}/C/common/inc/ ${CUDA_INCLUDE_DIRS})
ADD_EXECUTABLE(ssa_optimize_3d_cuda
  ssa_optimize_3d_cuda.cpp
)

TARGET_LINK_LIBRARIES(ssa_optimize_3d_cuda ssa_core types_3d g2o_stuff g2o_types_slam3d g2o_types_icp g2o_core g2o_solver_cholmod ${CHOLMOD_LIBRARIES} ${BLAS_LIBRARIES} g2o_solver_csparse g2o_csparse_extension flann ${CSPARSE_LIBRARY} ${QT_QTXML_LIBRARY} ${QT_QTCORE_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_cusparse_LIBRARY} solver_pcg_cuda pcg_cuda ${CUDA_SDK_LIBRARIES})
ENDIF(CUDA_FOUND AND SSA_USE_CUDA)

INSTALL(TARGETS ssa_optimize_3d
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

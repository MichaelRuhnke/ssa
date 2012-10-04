INCLUDE(${QT_USE_FILE})

# ${QT_INCLUDE_DIR} does not appear to work; therefore have to
# manually rig it; also have to add QtCore directly because qglviewer
# doesn't seem to handle the Qt directories properly

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CHOLMOD_INCLUDE_DIR} ${CSPARSE_INCLUDE_DIR} ${QT_QTCORE_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${QT_QTCORE_INCLUDE_DIR}/..)

# ADD_EXECUTABLE(ssa_dump_matrix
#   ssa_dump_matrix.cpp
# )
# 
# TARGET_LINK_LIBRARIES(ssa_dump_matrix ssa_core types_3d g2o_stuff g2o_types_slam3d g2o_core g2o_solver_cholmod ${CHOLMOD_LIBRARIES} ${BLAS_LIBRARIES} g2o_solver_csparse flann ${CSPARSE_LIBRARY} ${QT_QTCORE_LIBRARY})
# 
# INSTALL(TARGETS ssa_dump_matrix
#   RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
#   LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
#   ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
# )

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CHOLMOD_INCLUDE_DIR} ${CSPARSE_INCLUDE_DIR} ${QT_QTCORE_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${QT_QTCORE_INCLUDE_DIR}/..)

ADD_EXECUTABLE(ssa_constructor_benchmark
  ssa_constructor_benchmark.cpp
)

TARGET_LINK_LIBRARIES(ssa_constructor_benchmark ssa_core types_3d g2o_stuff g2o_types_slam3d g2o_types_sba g2o_core g2o_solver_cholmod ${CHOLMOD_LIBRARIES} ${BLAS_LIBRARIES} g2o_solver_csparse flann ${CSPARSE_LIBRARY} ${QT_QTCORE_LIBRARY})

INSTALL(TARGETS ssa_constructor_benchmark
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)


INCLUDE(${QT_USE_FILE})

# ${QT_INCLUDE_DIR} does not appear to work; therefore have to
# manually rig it; also have to add QtCore directly because qglviewer
# doesn't seem to handle the Qt directories properly

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CHOLMOD_INCLUDE_DIR} ${QT_QTCORE_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${QT_QTCORE_INCLUDE_DIR}/..)


ADD_EXECUTABLE(pcd_to_ssa
  pcd_to_ssa.cpp
)

TARGET_LINK_LIBRARIES(pcd_to_ssa  ssa_core types_3d g2o_types_slam3d g2o_core flann sensor_models ${Boost_SYSTEM_LIBRARY} ${PCL_COMMON_LIBRARY} ${PCL_IO_LIBRARY} ${PCL_FILTERS_LIBRARY} )

ADD_EXECUTABLE(pcd_to_ssa_rgb
  pcd_to_ssa_rgb.cpp
)

TARGET_LINK_LIBRARIES(pcd_to_ssa_rgb  ssa_core types_3d g2o_types_slam3d g2o_core flann sensor_models ${Boost_SYSTEM_LIBRARY} ${PCL_COMMON_LIBRARY} ${PCL_IO_LIBRARY} ${PCL_FILTERS_LIBRARY} )

#ADD_EXECUTABLE(ply_to_ssa_rgb
  #ply_to_ssa_rgb.cpp
#)

#TARGET_LINK_LIBRARIES(ply_to_ssa_rgb  ssa_core types_3d g2o_types_slam3d g2o_core flann sensor_models ${PCL_COMMON_LIBRARY} ${PCL_IO_LIBRARY} ${PCL_FILTERS_LIBRARY} )


ADD_EXECUTABLE(ssa_to_pcd
  ssa_to_pcd.cpp
)

TARGET_LINK_LIBRARIES(ssa_to_pcd  ssa_core types_3d g2o_types_slam3d g2o_core flann sensor_models ${Boost_SYSTEM_LIBRARY} ${PCL_COMMON_LIBRARY} ${PCL_IO_LIBRARY} ${PCL_FILTERS_LIBRARY} )


INSTALL(TARGETS pcd_to_ssa
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

FILE(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

INSTALL(FILES ${headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/pcl_binding)

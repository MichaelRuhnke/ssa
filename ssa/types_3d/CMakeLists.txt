ADD_LIBRARY(types_3d SHARED
types_ssa3d.cpp
vertex_point_xyzcov.h
vertex_point_xyzcov.cpp
edge_se3_xyzcov.h
edge_se3_xyzcov.cpp
edge_xyzcov_xyzcov.h
edge_xyzcov_xyzcov.cpp
observation_point_xyz_rgb_normal.h
observation_point_xyz_rgb_normal.cpp
)

SET_TARGET_PROPERTIES(types_3d PROPERTIES OUTPUT_NAME ${LIB_PREFIX}types_3d)

TARGET_LINK_LIBRARIES(types_3d)

INSTALL(TARGETS types_3d
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

FILE(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

INSTALL(FILES ${headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ssa/types_3d)

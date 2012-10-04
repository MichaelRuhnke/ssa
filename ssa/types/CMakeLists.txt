ADD_LIBRARY(types SHARED
types_ssa2d.cpp
vertex_point_xycov.h
vertex_point_xycov.cpp
edge_se2_xycov.h
edge_se2_xycov.cpp
edge_xycov_xycov.h
edge_xycov_xycov.cpp
)

SET_TARGET_PROPERTIES(types PROPERTIES OUTPUT_NAME ${LIB_PREFIX}types)

TARGET_LINK_LIBRARIES(types)

INSTALL(TARGETS types
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

FILE(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

INSTALL(FILES ${headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ssa/types)

INCLUDE(${QT_USE_FILE})

# ${QT_INCLUDE_DIR} does not appear to work; therefore have to
# manually rig it; also have to add QtCore directly because qglviewer
# doesn't seem to handle the Qt directories properly

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CHOLMOD_INCLUDE_DIR} ${QT_QTCORE_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${QT_QTCORE_INCLUDE_DIR}/..)


ADD_EXECUTABLE(ssa_from_posegraph2d
  posegraph2ssa.cpp
  import_posegraph2d.cpp
)

TARGET_LINK_LIBRARIES(ssa_from_posegraph2d ssa_core types sensor_models g2o_core g2o_types_slam2d ${QT_QTXML_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTCORE_LIBRARY})


INSTALL(TARGETS ssa_from_posegraph2d
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

ADD_EXECUTABLE(ssa_to_posegraph2d
  ssa2posegraph.cpp
)

TARGET_LINK_LIBRARIES(ssa_to_posegraph2d ssa_core types sensor_models g2o_core g2o_types_slam2d ${QT_QTXML_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTCORE_LIBRARY})


INSTALL(TARGETS ssa_to_posegraph2d
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

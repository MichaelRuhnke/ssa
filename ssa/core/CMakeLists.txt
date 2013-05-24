INCLUDE(${QT_USE_FILE})

QT4_WRAP_CPP(MOC_SOURCES ssa_sparse_optimizer.h)

# set(sparse_optimizer_SRCS
#     ssa_sparse_optimizer.cpp
# )
#
# qt4_automoc(${sparse_optimizer_SRCS})


# ${QT_INCLUDE_DIR} does not appear to work; therefore have to
# manually rig it; also have to add QtCore directly because qglviewer
# doesn't seem to handle the Qt directories properly

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR} ${QT_QTCORE_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${QT_QTCORE_INCLUDE_DIR}/..)

ADD_LIBRARY(ssa_core SHARED
ssa_sparse_optimizer.h
ssa_sparse_optimizer.cpp
${MOC_SOURCES}
)

SET_TARGET_PROPERTIES(ssa_core PROPERTIES OUTPUT_NAME ${LIB_PREFIX}ssa_core)

TARGET_LINK_LIBRARIES(ssa_core ${QT_QTXML_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTCORE_LIBRARY})

INSTALL(TARGETS ssa_core
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

FILE(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

INSTALL(FILES ${headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ssa/core)

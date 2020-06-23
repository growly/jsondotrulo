#include(ExternalProject)

set(GFLAGS_IS_SUBPROJECT TRUE)
add_subdirectory(gflags)

#add_subdirectory(glog)

add_subdirectory(json)

#ExternalProject_Add(
#  glog-project
#  DEPENDS gflags
#  PREFIX ${CMAKE_BINARY_DIR}/glog
#  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/glog
#  CMAKE_ARGS -Dgflags_DIR=${CMAKE_CURRENT_BINARY_DIR}/gflags
#)

set(WITH_GFLAGS FALSE)
set(gflags_FOUND TRUE CACHE BOOL "")
set(gflags_DIR ${CMAKE_CURRENT_BINARY_DIR}/gflags CACHE FILEPATH "")
set(gflags_NAMESPACE "gflags" CACHE INTERNAL "Namespace for gflags")
set(gflags_VERSION_STRING "2.2.2" CACHE STRING "")
set(gflags_VERSION_MAJOR  2)
set(gflags_VERSION_MINOR  2)
set(gflags_VERSION_PATCH  0)
set(gflags_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/gflags/include")
set(gflags_LIBRARIES "gflags_nothreads_static")
export(TARGETS ${gflags_LIBRARIES}
       NAMESPACE gflags::
       FILE "${CMAKE_CURRENT_BINARY_DIR}/gflags/gflags-targets.cmake")
add_subdirectory(glog ${CMAKE_BINARY_DIR}/glog EXCLUDE_FROM_ALL)

include(ExternalProject)

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

set(gflags_DIR ${CMAKE_CURRENT_BINARY_DIR}/gflags CACHE FILEPATH "")
set(gflags_NAMESPACE "gflags" CACHE INTERNAL "Namespace for gflags")
add_subdirectory(glog ${CMAKE_BINARY_DIR}/glog EXCLUDE_FROM_ALL)

set(NAME core)

set(HEADERS
  Edges.hpp
  Faces.hpp
  Graph.hpp
  HalfEdges.hpp
  Heap.hpp
  Geometry.hpp
  Optimization.hpp
  Partition.hpp
  PolygonMesh.hpp
  PolygonMeshTest.hpp
  Variable.hpp
) # HEADERS    

set(SOURCES
  Edges.cpp
  Faces.cpp
  Graph.cpp
  Geometry.cpp
  HalfEdges.cpp
  Heap.cpp
  Optimization.cpp
  Partition.cpp
  PolygonMesh.cpp
  PolygonMeshTest.cpp
  Variable.cpp
) # SOURCES

add_library(${NAME}
    ${HEADERS}
    ${SOURCES}
    ${FORMS}
    ${ASSETS}
) # add_library

find_package(Eigen3 REQUIRED)


target_compile_features(${NAME} PRIVATE cxx_lambdas)

target_link_libraries(${NAME} ${LIB_LIST} Eigen3::Eigen)


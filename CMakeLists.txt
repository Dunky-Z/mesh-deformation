cmake_minimum_required (VERSION 3.8)

project ("laplacian-deformation")

#Eigen
find_package(Eigen3 REQUIRED)

#SURFACE MESH
include_directories("$ENV{Surface_Mesh_DIR}/include")
include_directories("D:/Lib/pmp-library-debug/include")

add_executable(laplacian-deformation src/main.cpp )

target_link_libraries(laplacian-deformation Eigen3::Eigen)
target_link_libraries(laplacian-deformation optimized "D:/Lib/pmp-library-release/lib/pmp.lib")


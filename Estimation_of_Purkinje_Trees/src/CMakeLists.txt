
cmake_minimum_required(VERSION 2.8)

PROJECT(ImageUtilites)

# Find ITK
find_package( ITK REQUIRED )
include( ${USE_ITK_FILE} )
#SET(ITK_LIBRARIES)

# Find VTK
find_package( VTK REQUIRED )
include( ${USE_VTK_FILE} )
#SET(VTK_LIBRARIES vtkCommon vtkIO vtkGraphics)

SET(Required_Packages
  ITK VTK
)


SET(Libraries
${ITK_LIBRARIES} ${VTK_LIBRARIES})




ADD_EXECUTABLE(Voxelize Voxelize.cpp) 
TARGET_LINK_LIBRARIES(Voxelize ${Libraries})

ADD_EXECUTABLE(ProjectPointsToSurface ProjectPointsToSurface.cxx)
TARGET_LINK_LIBRARIES(ProjectPointsToSurface ${Libraries})

add_executable( SmoothPolyLine SmoothPolyLine.cpp )
target_link_libraries( SmoothPolyLine ${Libraries})

add_executable(ApplyScalarsToVolume ApplyScalarsToVolume.cxx )
target_link_libraries( ApplyScalarsToVolume ${Libraries})

add_executable(DetectGradientSingularities DetectGradientSingularities.cpp DToptimo3d.cpp)
target_link_libraries( DetectGradientSingularities ${Libraries})

add_executable(GroupPoints GroupPoints.cxx)
target_link_libraries( GroupPoints ${Libraries})

add_executable(ComputePKTree ComputePKTree.cpp DToptimo3d.cpp)
target_link_libraries( ComputePKTree ${Libraries})

ADD_EXECUTABLE(edgedetect3D edgedetect3D.cpp)
TARGET_LINK_LIBRARIES(edgedetect3D ${Libraries})








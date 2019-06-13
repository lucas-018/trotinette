cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(Shapes)

find_package(VTK COMPONENTS
  vtkCommonColor
  vtkCommonCore
  vtkCommonDataModel
  vtkFiltersCore
  vtkInteractionStyle
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingOpenGL2 QUIET)
if (NOT VTK_FOUND)
  message("Skipping Shapes: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(Shapes MACOSX_BUNDLE app.cpp )
  target_link_libraries(Shapes PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(Shapes MACOSX_BUNDLE app.cpp )
  target_link_libraries(Shapes PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS Shapes
    MODULES ${VTK_LIBRARIES}
    )
endif ()

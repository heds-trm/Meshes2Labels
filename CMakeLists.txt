PROJECT(Meshes2Labels)

cmake_minimum_required(VERSION 3.10)

cmake_policy(SET CMP0054 NEW)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
include(heds_apps_macros)

find_package(VTK REQUIRED)
if(NOT VTK_FOUND)
  message(FATAL_ERROR "VTK not found!")
endif()

find_package(ITK REQUIRED)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
else(ITK_FOUND)
  message(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
endif(ITK_FOUND)

# for command line management
FIND_PATH( TCLAP_INCLUDE_DIR
 tclap
)

OPTION(USE_CUSTOM_OBJ_WRITER "If true, custom obj writer is used" OFF) 
MARK_AS_ADVANCED(USE_CUSTOM_OBJ_WRITER)

# installation of binaries
SET(INSTALL_PATH "C:/Programming/Exec" CACHE PATH "Install path m2l of executable" )

IF(USE_CUSTOM_OBJ_WRITER)
add_executable(m2l 3rdparty/gdiam/gdiam.h 3rdparty/gdiam/gdiam.cxx m2l.cxx io.cxx io.h vtkOBJWriter.h vtkOBJWriter.cxx)
ELSE(USE_CUSTOM_OBJ_WRITER)
add_executable(m2l 3rdparty/gdiam/gdiam.h 3rdparty/gdiam/gdiam.cxx m2l.cxx io.cxx io.h)
ENDIF(USE_CUSTOM_OBJ_WRITER)

target_include_directories(m2l SYSTEM PRIVATE  ${TCLAP_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/3rdparty/gdiam)
disable_usual_warnings(m2l)

set(vtk_req_targets VTK::CommonTransforms VTK::ImagingColor VTK::FiltersModeling VTK::FiltersGeometry VTK::FiltersSources VTK::FiltersParallel VTK::IOLegacy VTK::CommonColor VTK::IOGeometry VTK::IOPLY VTK::RenderingCore VTK::ImagingMorphological VTK::ImagingStencil)
target_link_libraries(m2l ${vtk_req_targets} ${ITK_LIBRARIES})

# with vtk > 8.2 we must use the following
# otherwise it may fail at *runtime*
vtk_module_autoinit(
  TARGETS m2l
  MODULES ${VTK_LIBRARIES})

if( EXISTS ${INSTALL_PATH} )
    install (TARGETS m2l DESTINATION ${INSTALL_PATH})
endif()


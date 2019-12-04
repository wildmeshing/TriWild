################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.3)
project(TriWild)

################################################################################
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

################################################################################
# Build static binaries
#set(BUILD_SHARED_LIBS OFF)
set(BUILD_STATIC_RUNTIME OFF)
set(BUILD_DYNAMIC_RUNTIME ON)

# Use folder in Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Export compile flags (used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Generate position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#ADD_DEFINITIONS("-D_DLL")

include(TriWildDownloadExternal)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	# If we use static library, we link with the static C++ runtime.
	foreach(config ${CMAKE_CONFIGURATION_TYPES})
		string(TOUPPER ${config} config)
		string(REPLACE /MD /MT CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
		string(REPLACE /MD /MT CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
		string(REPLACE /MDd /MTd CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
		string(REPLACE /MDd /MTd CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
	endforeach()
endif()

triwild_download_aabbcc()
triwild_download_geogram()
triwild_download_igl()
triwild_download_nlopt()

################################################################################
# 3rd party libraries
################################################################################

# libigl library
# libigl
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_XML               "Use XML"            OFF)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       OFF)

# gui
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  OFF)

find_package(LIBIGL REQUIRED QUIET)
find_package(GMPTriWild)

IF(NOT ${GMP_FOUND})
        MESSAGE(FATAL_ERROR "Cannot find GMP")
ENDIF()


################################################################################

set(SOURCE_FILES
		extern/aabbcc/src/AABB.h extern/aabbcc/src/AABB.cc
		extern/pymesh/MshLoader.h extern/pymesh/MshLoader.cpp
		extern/pymesh/MshSaver.h extern/pymesh/MshSaver.cpp
		extern/pymesh/Exception.h
		extern/CLI11.hpp

		src/triwild/Logger.h
		src/triwild/Args.cpp src/triwild/Args.h
		src/triwild/meshio.hpp src/triwild/meshio.cpp

		src/triwild/auto_p_bases.hpp src/triwild/auto_p_bases.cpp
		src/triwild/auto_det_checker.hpp src/triwild/auto_det_checker.cpp
		src/triwild/Curves.h src/triwild/Curves.cpp
		src/triwild/CurvedTriUntangler.hpp src/triwild/CurvedTriUntangler.cpp

		src/triwild/Rational.h
		src/triwild/Point_2.h src/triwild/Point_2f.h
		src/triwild/TrimeshElements.h
		src/triwild/triangulation.cpp src/triwild/triangulation.h
		src/triwild/AMIPS.h src/triwild/AMIPS.cpp
		src/triwild/optimization.cpp src/triwild/optimization.h
		src/triwild/edge_collapsing.cpp src/triwild/edge_collapsing.h
		src/triwild/edge_splitting.cpp src/triwild/edge_splitting.h
		src/triwild/edge_swapping.cpp src/triwild/edge_swapping.h
		src/triwild/vertex_smoothing.cpp src/triwild/vertex_smoothing.h

		src/triwild/FeatureElements.cpp src/triwild/FeatureElements.h
		src/triwild/feature_preprocessing.cpp src/triwild/feature_preprocessing.h
		src/triwild/feature.cpp src/triwild/feature.h

		src/triwild/reference_triangle.cpp src/triwild/reference_triangle.h
		src/triwild/do_triwild.cpp src/triwild/do_triwild.h)

add_definitions(-D_REF_VS=\"${CMAKE_CURRENT_SOURCE_DIR}/src/triwild/reference_triangle_vertices.txt\")
add_definitions(-D_REF_FS=\"${CMAKE_CURRENT_SOURCE_DIR}/src/triwild/reference_triangle_faces.txt\")

option(NLOPT_PYTHON OFF)
option(NLOPT_GUILE OFF)
option(NLOPT_SWIG OFF)
option(NLOPT_OCTAVE OFF)
option(NLOPT_MATLAB OFF)
option(NLOPT_TESTS OFF)
option(NLOPT_FORTRAN OFF)
option(BUILD_SHARED_LIBS OFF)
add_subdirectory(${THIRD_PARTY_DIR}/nlopt)


add_library(triwild_lib  ${SOURCE_FILES})

# Geogram library
include(geogram)


if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	target_compile_definitions(triwild_lib PUBLIC -D_ENABLE_EXTENDED_ALIGNED_STORAGE)

	# If we use static library, we link with the static C++ runtime.
	foreach(config ${CMAKE_CONFIGURATION_TYPES})
		string(TOUPPER ${config} config)
		string(REPLACE /MD /MT CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
		string(REPLACE /MD /MT CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
		string(REPLACE /MDd /MTd CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
		string(REPLACE /MDd /MTd CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
	endforeach()
endif()


target_include_directories(triwild_lib PUBLIC ${GMP_INCLUDE_DIRS})
target_include_directories(triwild_lib PUBLIC "src/triwild")
target_link_libraries(triwild_lib PUBLIC geogram igl::core ${GMP_LIBRARIES} nlopt)

if(NOT TARGET json)
	triwild_download_json()
	add_library(json INTERFACE)
	target_include_directories(json SYSTEM INTERFACE ${THIRD_PARTY_DIR}/json/include)
endif()
target_link_libraries(triwild_lib PUBLIC geogram igl::core ${GMP_LIBRARIES} nlopt json)




#add_executable(TriWild main.cpp UIPanel.hpp UIPanel.cpp)
add_executable(TriWild src/main.cpp)
target_link_libraries(TriWild triwild_lib)

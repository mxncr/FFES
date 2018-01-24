cmake_minimum_required(VERSION 3.0)
project(Fast-Finite-Element-Sampling LANGUAGES C CXX)

# Generate a compile_commands.json file for autocompletion
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")

# C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")


### Third parties ###
include_directories(${CMAKE_SOURCE_DIR}/third_party)
add_subdirectory(${CMAKE_SOURCE_DIR}/third_party/geogram_basic)

set(GLFW_FOLDER  ${CMAKE_SOURCE_DIR}/third_party/glfw)
set(GLM_FOLDER ${CMAKE_SOURCE_DIR}/third_party/glm)   # projection matrix
set(GLSW_FOLDER ${CMAKE_SOURCE_DIR}/third_party/glsw) # shader loading code
set(IMGUI_FOLDER ${CMAKE_SOURCE_DIR}/third_party/imgui)
set(GL3W_FOLDER  ${CMAKE_SOURCE_DIR}/third_party/gl3w)
set(IMGUI_HEADER ${IMGUI_FOLDER}/imgui_impl_glfw_gl3.h)

include_directories(${GLFW_FOLDER})
include_directories(${IMGUI_FOLDER})
include_directories(${GL3W_FOLDER})
include_directories(${GLSW_FOLDER})
include_directories(${GLM_FOLDER})
include_directories(${CMAKE_SOURCE_DIR})

find_library(GLFW_LIB glfw)
if(NOT GLFW_LIB)
    message(STATUS "glfw not found, will be built")
    add_subdirectory(${CMAKE_SOURCE_DIR}/third_party/glfw)
else()
    message(STATUS "glfw found")
endif()

add_library(glsw STATIC third_party/glsw/glsw.c third_party/glsw/bstrlib.c)
add_library(imgui STATIC third_party/imgui/imgui_impl_glfw_gl3.cpp third_party/imgui/imgui.cpp
    third_party/imgui/imgui_draw.cpp third_party/gl3w/GL/gl3w.c)
target_link_libraries(imgui PUBLIC glfw)

# OpenGL libraries
set(GL_LIBRARIES GL)
if(WIN32)
	set(GL_LIBRARIES opengl32)
endif()

set(DL_LIBRARIES dl)
if(WIN32)
    set(DL_LIBRARIES ${CMAKE_DL_LIBS})
endif()

find_package(Threads)

### Fast Finite Element Sampling ###
include_directories(${CMAKE_SOURCE_DIR}/src)
# compute-distance executable
add_executable(ffes-distance src/bin/compute_distance.cpp src/ffes/distance.cpp src/ffes/field_sampling.cpp src/ffes/field_io.cpp src/ffes/utils.cpp)
target_link_libraries(ffes-distance ${GL_LIBRARIES} ${DL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} glfw glsw imgui geogram_basic)
# copy the shaders from src/shaders to build_directory/shaders
add_custom_command(TARGET ffes-distance POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E remove_directory $<TARGET_FILE_DIR:ffes-distance>/shaders)
add_custom_command(TARGET ffes-distance POST_BUILD
               COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/src/shaders $<TARGET_FILE_DIR:ffes-distance>/shaders)
if(WIN32)
	target_compile_definitions(ffes-distance PUBLIC -D_USE_MATH_DEFINES)
endif()


# convert-from-mfem: executable to convert a mfem::Mesh and mfem::GridFunction to
#                    FieldObject binary format
#  only if mfem is found in ./third_party/mfem (can be a symlink)
if(EXISTS ${CMAKE_SOURCE_DIR}/third_party/mfem)
    set(MFEM_INCLUDES ${CMAKE_SOURCE_DIR}/third_party/mfem)
    set(MFEM_LIB ${CMAKE_SOURCE_DIR}/third_party/mfem/libmfem.a)
    message(STATUS "mfem folder found: ${MFEM_INCLUDES}")
    if(EXISTS ${MFEM_LIB})
        include_directories(${MFEM_INCLUDES})
        add_executable(convert-from-mfem src/bin/convert_from_mfem.cpp)
        target_link_libraries(convert-from-mfem ${MFEM_LIB})
    else()
        message(STATUS "libmfem.a not found, please build mfem")
    endif()
endif()

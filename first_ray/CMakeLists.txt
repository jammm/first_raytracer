cmake_minimum_required(VERSION 3.10)

# set the project name
project(first_ray VERSION 0.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -DAVX2 -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/")

find_package(OpenGL REQUIRED)

#add GLEW as STATIC library
add_library(GLEW STATIC dependencies/glew-2.1.0/src/glew.c dependencies/glew-2.1.0/include/GL/glew.h)

# add dependencies
# add_subdirectory(dependencies/cpp-taskflow-2.1.0)

FILE(GLOB MY_CPP_FILES *.cpp)

# some more variables
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

#add include directories
include_directories(
                    ${CMAKE_CURRENT_SOURCE_DIR}
                    dependencies/cpp-taskflow-2.1.0
                    dependencies/assimp-5.0.1/include/
                    dependencies/glfw-3.3.2/include/
                    dependencies/glew-2.1.0/include/ 
                    )

# set some variables
set(BUILD_SHARED_LIBS ON)
set(ASSIMP_BUILD_TESTS OFF)
set(ASSIMP_BUILD_TESTS OFF)
set(ASSIMP_BUILD_ASSIMP_VIEW OFF)
set(ASSIMP_DOUBLE_PRECISION ON)
set(ASSIMP_BUILD_ASSIMP_TOOLS OFF)

add_subdirectory(dependencies/assimp-5.0.1)
add_subdirectory(dependencies/glfw-3.3.2)

set(extLibs
    ${extLibs}
    ${GLFW_LIBRARIES}
    )

# add the executable
add_executable(main ${MY_CPP_FILES})

target_link_libraries(main assimp GLEW glfw ${OPENGL_LIBRARIES} ${GLEW_LIBRARIES} ${GLFW_LIBRARIES} ${ASSIMP_LIBRARIES} -lpthread)

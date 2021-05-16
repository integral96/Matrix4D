cmake_minimum_required(VERSION 3.5)

project(Matrix4D ${PROJECT_NAME} LANGUAGES CXX )

#set(CMAKE_CXX_COMPILER "/opt/intel/oneapi/compiler/2021.2.0/linux/bin/dpcpp")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(Boost_USE_MULTITHREADED ON)



find_package(Boost ${Boost_VERSION_MAJOR} REQUIRED timer regex )


find_package(TBB REQUIRED tbb)

find_package(OpenGL REQUIRED)

find_package(glfw3 REQUIRED)

find_package(GLUT REQUIRED)

if(NOT DEFINED OpenGL_GL_PREFERENCE)
  set(OpenGL_GL_PREFERENCE "LEGACY")
endif()

file(GLOB SOURCES
    include/*.hpp
    sources/*.cpp)

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ltbb -pedantic -pthread -Wall -Wextra -lGLEW -lglut -lGL -lGLU")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth=2900")

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries (${PROJECT_NAME} ${TBB_IMPORTED_TARGETS} ${Boost_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})

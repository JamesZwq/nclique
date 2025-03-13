cmake_minimum_required(VERSION 3.10)
project(pivoter CXX)

if(APPLE)
    set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
    set(TBB_ROOT "/Users/zhangwenqian/myBin/tbb")
    set(TBB_LIBRARY "${TBB_ROOT}/lib/libtbb.dylib")
else()
    set(CMAKE_CXX_COMPILER "/usr/bin/g++")
    set(TBB_ROOT "~/bin/tbb")
    set(TBB_LIBRARY "${TBB_ROOT}/lib/libtbb.so")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# add Boost
find_package(Boost REQUIRED COMPONENTS serialization)
include_directories(${Boost_INCLUDE_DIRS})
#print boost include directories
message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")

# 添加 TBB 的头文件搜索路径
include_directories("${TBB_ROOT}/include")
include_directories(src)

file(GLOB_RECURSE SOURCES src/*.cpp src/*/*.cpp src/*/*.h)

add_executable(degeneracy_cliques ${SOURCES})

# 指定 TBB 库的路径（根据平台，macOS 一般是 .dylib）
target_link_libraries(degeneracy_cliques PRIVATE ${TBB_LIBRARY} ${Boost_LIBRARIES})
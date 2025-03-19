cmake_minimum_required(VERSION 3.10)
project(pivoter CXX)

if(APPLE)
    set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
    set(TBB_ROOT "/Users/zhangwenqian/myBin/tbb")
    set(TBB_LIBRARY "${TBB_ROOT}/lib/libtbb.dylib")
else()
    set(CMAKE_CXX_COMPILER "/usr/bin/g++")
    get_filename_component(HOME_DIR "~" REALPATH)
    set(TBB_ROOT "${HOME_DIR}/bin/tbb")
    set(TBB_LIBRARY "${TBB_ROOT}/lib/libtbb.so")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# add Boost
#find_package(Boost REQUIRED COMPONENTS serialization)
find_package(Boost REQUIRED CONFIG COMPONENTS serialization)
include_directories(${Boost_INCLUDE_DIRS})
#print boost include directories
message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")

# 添加 TBB 的头文件搜索路径
include_directories("${TBB_ROOT}/include")
include_directories(src)

file(GLOB_RECURSE ALL_SOURCES src/*.cpp)
set(EXCLUDE_SRCS
        ${CMAKE_SOURCE_DIR}/src/main.cpp
        ${CMAKE_SOURCE_DIR}/src/degeneracy_cliques.cpp
)
list(REMOVE_ITEM ALL_SOURCES ${EXCLUDE_SRCS})
set(COMMON_SOURCES ${ALL_SOURCES})

# 指定每个可执行文件
add_executable(degeneracy_cliques ${COMMON_SOURCES} src/degeneracy_cliques.cpp)
add_executable(main ${COMMON_SOURCES} src/main.cpp)

# 链接库
target_link_libraries(degeneracy_cliques PRIVATE ${TBB_LIBRARY} ${Boost_LIBRARIES})
target_link_libraries(main PRIVATE ${TBB_LIBRARY} ${Boost_LIBRARIES})
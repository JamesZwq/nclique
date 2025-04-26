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

    cmake_policy(SET CMP0074 NEW)

    # 设置 Boost 的目录
    set(BOOST_ROOT "/home/wenqianz/bin/boost")
    set(BOOST_INCLUDEDIR "/home/wenqianz/bin/boost/include")
    set(BOOST_LIBRARYDIR "/home/wenqianz/bin/boost/lib")
    set(BOOST_NO_SYSTEM_PATHS TRUE)
    message(STATUS "Linux Boost root: ${BOOST_ROOT}")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

find_package(Boost REQUIRED COMPONENTS serialization)
#find_package(Boost REQUIRED CONFIG COMPONENTS serialization)
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

#if debug
if (CMAKE_BUILD_TYPE MATCHES "Debug")
    message(STATUS "Debug build type detected.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    add_compile_options(
            -Wall
            -Wextra
    #        -Werror
    #        -fsanitize=address
            -pedantic-errors
            -Wconversion
            -Wsign-conversion
            -Wdouble-promotion
            -Wcast-align
            -Wformat=2
            -Wuninitialized
            -Wnull-dereference
            -Wnon-virtual-dtor
            -Woverloaded-virtual
    #        -Wold-style-cast
            -Wzero-as-null-pointer-constant
            -Wno-self-assign-overloaded
            -Wno-self-move
            -Wsuggest-override
            -fstack-protector-strong
            -ferror-limit=1000
    )
    set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -fsanitize=address -fno-omit-frame-pointer -g")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

if (CMAKE_BUILD_TYPE MATCHES "Release")
    message(STATUS "Release build type detected.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    add_compile_options(
            -march=native
    )
    set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -O3")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O3")
endif()

message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost libraries: ${Boost_LIBRARIES}")
# 指定每个可执行文件
add_executable(degeneracy_cliques ${COMMON_SOURCES} src/degeneracy_cliques.cpp
        src/tree/NucleusCoreDecomposition.cpp
        src/tree/NucleusCoreDecomposition.h
        src/dataStruct/CliqueHashMap.h
        src/graph/EdgeList.h
        src/graph/EdgeList.cpp
        src/graph/DynamicGraph.cpp
        src/graph/DynamicGraph.h)
add_executable(main ${COMMON_SOURCES} src/main.cpp)

# 链接库
target_link_libraries(degeneracy_cliques PRIVATE ${TBB_LIBRARY} ${Boost_LIBRARIES})
target_link_libraries(main PRIVATE ${TBB_LIBRARY} ${Boost_LIBRARIES})
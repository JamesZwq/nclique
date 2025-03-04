cmake_minimum_required(VERSION 3.10)
project(pivoter CXX)

# 使用 C++23 标准
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 让 CMake 自动创建 `bin` 目录（防止链接错误）
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# 头文件搜索路径
include_directories(src)

# 获取所有的 .cpp 文件（不需要手动列举）
file(GLOB_RECURSE SOURCES src/*.cpp)

# 生成可执行文件（包含所有 .cpp 文件）
add_executable(degeneracy_cliques ${SOURCES})
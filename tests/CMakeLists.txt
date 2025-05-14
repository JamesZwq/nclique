# ───────────────────────────────────────────────────────────────
# tests/CMakeLists.txt   —— 仅负责编译 / 运行测试
# ───────────────────────────────────────────────────────────────

# 0️⃣ 与顶层保持一致的编译器 & 第三方库路径
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

    # Boost 手工路径（与顶层保持一致）
    set(BOOST_ROOT       "/home/wenqianz/bin/boost")
    set(BOOST_INCLUDEDIR "/home/wenqianz/bin/boost/include")
    set(BOOST_LIBRARYDIR "/home/wenqianz/bin/boost/lib")
    set(BOOST_NO_SYSTEM_PATHS TRUE)
    message(STATUS "Linux Boost root: ${BOOST_ROOT}")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1️⃣ 查找 Boost.Serialization
find_package(Boost REQUIRED COMPONENTS serialization)
include_directories(${Boost_INCLUDE_DIRS})
message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")

# 2️⃣ 把 TBB 和工程头文件加入搜索路径（与主工程一致）
include_directories("${TBB_ROOT}/include")
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)

# 3️⃣ 组装测试可执行文件
#    除了 test_bronkerbosch.cpp，还把运行时需要的实现文件列进来
add_executable(test_bronkerbosch
        test_bronkerbosch.cpp
        ${PROJECT_SOURCE_DIR}/src/Global/Global.cpp
        ${PROJECT_SOURCE_DIR}/src/tree/BronKerbosch.cpp
        ${PROJECT_SOURCE_DIR}/src/graph/DynamicGraph.cpp
)

# 4️⃣ 链接 Boost & TBB 动态库（与主程序相同写法）
target_link_libraries(test_bronkerbosch
        PRIVATE
        ${TBB_LIBRARY}
        ${Boost_LIBRARIES}
)

# 5️⃣ 可选：打开 ASan / 调试选项（保持你原来的习惯）
target_compile_options(test_bronkerbosch PRIVATE -fsanitize=address -g)
target_link_options   (test_bronkerbosch PRIVATE -fsanitize=address)

# 6️⃣ 注册到 CTest
add_test(
        NAME BronKerboschTest
        COMMAND $<TARGET_FILE:test_bronkerbosch>
)
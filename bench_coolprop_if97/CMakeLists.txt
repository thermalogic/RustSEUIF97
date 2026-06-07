cmake_minimum_required(VERSION 3.15)
project(bench_if97 LANGUAGES CXX)

# 设置默认构建类型为 Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 添加可执行文件
add_executable(bench_if97 bench_if97.cpp)

# 添加编译优化选项
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # GCC/Clang 优化选项
    target_compile_options(bench_if97 PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>
        $<$<CONFIG:Release>:-flto>
        $<$<CONFIG:Release>:-funroll-loops>
        $<$<CONFIG:Release>:-ffast-math>
    )
    target_link_options(bench_if97 PRIVATE
        $<$<CONFIG:Release>:-flto>
    )
elif(MSVC)
    # MSVC 优化选项
    target_compile_options(bench_if97 PRIVATE
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/GL>
        $<$<CONFIG:Release>:/arch:AVX2>
        $<$<CONFIG:Release>:/fp:fast>
    )
    target_link_options(bench_if97 PRIVATE
        $<$<CONFIG:Release>:/LTCG>
    )
endif()

# 包含头文件目录
include_directories(${PROJECT_SOURCE_DIR})

# 链接 Rust SEUIF97 共享库
# Windows 下使用 MSVC 编译时，库文件为 seuif97.lib
if(WIN32)
    # 设置库搜索路径
    set(SEUIF97_LIB_PATH "${PROJECT_SOURCE_DIR}/../target/release")
    
    # 查找并链接 seuif97 库
    find_library(SEUIF97_LIBRARY 
        NAMES seuif97 libseuif97
        PATHS ${SEUIF97_LIB_PATH}
        REQUIRED
    )
    
    target_link_libraries(bench_if97 PRIVATE ${SEUIF97_LIBRARY})
    
    # 将 DLL 复制到输出目录（运行时需要）
    add_custom_command(TARGET bench_if97 POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${SEUIF97_LIB_PATH}/seuif97.dll"
            $<TARGET_FILE_DIR:bench_if97>
        COMMENT "Copying seuif97.dll to output directory"
    )
else()
    # Linux 系统设置
    set(SEUIF97_LIB_PATH "${PROJECT_SOURCE_DIR}/../target/release")
    find_library(SEUIF97_LIBRARY 
        NAMES seuif97
        PATHS ${SEUIF97_LIB_PATH}
        REQUIRED
    )
    target_link_libraries(bench_if97 PRIVATE ${SEUIF97_LIBRARY} m)
    set_target_properties(bench_if97 PROPERTIES
        INSTALL_RPATH "${SEUIF97_LIB_PATH}"
        BUILD_WITH_INSTALL_RPATH ON
    )
endif()

# 设置输出目录
set_target_properties(bench_if97 PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)
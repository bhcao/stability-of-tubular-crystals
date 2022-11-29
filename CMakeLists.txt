cmake_minimum_required(VERSION 3.16)
project(nano VERSION 0.3)

# 文件头
include_directories(include)

# 是否使用 KOKKOS
option(USE_KOKKOS "是否使用 KOKKOS（默认否）" OFF)

if(${USE_KOKKOS})
    message("KOKKOS 选项已开启")
    # Kokkos 文件夹，将第一个替换成自己的文件
    set(Kokkos_DIR /home/bovera/libkokkos/lib/cmake/Kokkos)
    find_package(Kokkos REQUIRED)
    add_compile_options(-DUSE_KOKKOS)
endif()

set(CMAKE_BUILD_TYPE "Debug")

# 源文件
file(GLOB SOURCES "src/*.cpp")
add_executable(nano ${SOURCES})

if(${USE_KOKKOS})
    target_link_libraries(nano Kokkos::kokkos)
endif()
cmake_minimum_required(VERSION 3.16)
project(nano VERSION 0.2.0)

# 文件头
include_directories(include)

# 编译选项，是否使用 KOKKOS 和是否生成 Python 库
option(USE_KOKKOS "是否使用 Kokkos" OFF)
option(PYTHON_LIB "是否编译成 Python 库" OFF)
option(DEBUG "是否调试" OFF)

if(${USE_KOKKOS})
    message("Kokkos 选项已开启")
    # Kokkos 文件夹，将第一个替换成自己的文件
    set(Kokkos_DIR /home/bovera/libkokkos/lib/cmake/Kokkos)
    find_package(Kokkos REQUIRED)
    add_compile_options(-DUSE_KOKKOS)
endif()

if(${PYTHON_LIB})
    message("Python 库选项已开启")
    find_package(Python3 COMPONENTS Interpreter Development)
    add_compile_options(-fPIC) # -fPIC 选项还需要在 Kokkos 编译时开启
endif()

if(${DEBUG})
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# 源文件，选择生成可执行文件件还是 Python 链接库
file(GLOB SOURCES "src/*.cpp")
if(${PYTHON_LIB})
    add_library(nano MODULE ${SOURCES} pythonc.cpp)
else()
    add_executable(nano ${SOURCES} main.cpp)
endif()

# 根据是否调用 Kokkos 和是否生成 Python 库决定链接库
if(${USE_KOKKOS} AND ${PYTHON_LIB})
    target_link_libraries(nano Kokkos::kokkos Python3::Python)
elseif(${USE_KOKKOS})
    target_link_libraries(nano Kokkos::kokkos)
elseif(${PYTHON_LIB})
    target_link_libraries(nano Python3::Python)
endif()
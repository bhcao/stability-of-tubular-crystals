cmake_minimum_required(VERSION 3.16)
project(nano VERSION 0.3)

# 文件头
include_directories(include)

# 源文件
file(GLOB SOURCES "src/*.cpp")
add_executable(nano ${SOURCES})
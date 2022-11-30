cmake_minimum_required(VERSION 3.16)
project(nanodoc VERSION 0.2.0)

find_package(Doxygen)

# set input and output files
set(DOXYGEN_IN ${CMAKE_CURRENT_LIST_DIR}/Doxyfile.in)
set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

# request to configure the file
configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

# note the option ALL which allows to build the docs together with the application
add_custom_target( docs ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CUREENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM
)
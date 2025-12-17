cmake_minimum_required(VERSION 3.25.1)

# -----------------------------------------
# Name of the project
set(PROJECT F2-Cmax-revisited)
# Version of the project
set(VERSION_PROJECT 1.0.0 )
# Language of the project
set(PROJECT_LANGUAGES CXX)
project(${PROJECT} VERSION ${VERSION_PROJECT} LANGUAGES ${PROJECT_LANGUAGES})

# -----------------------------------------
# files .cpp
# -----------------------------------------
# source for PROJECT (create new source for new lib)
set(SRCS_F2_Cmax
        ${PROJECT_SOURCE_DIR}/src/Parser.cpp
        ${PROJECT_SOURCE_DIR}/src/Instance.cpp
)
set(SRCS_PROJECT_MAIN
        ${PROJECT_SOURCE_DIR}/src/main.cpp
)
# ---------------------------------------

# ---------------------------------------
# files .h
# -----------------------------------------
set(HEADERS_PROJECT
        ${PROJECT_SOURCE_DIR}/include/
)
# ---------------------------------------

# Flags compilation
# compilation parameters
if (UNIX)
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-deprecated-copy -Wuninitialized -Wmaybe-uninitialized -pthread ")
    set(SANITIZER_FLAGS
            -fsanitize=address
            -fsanitize=leak
            -fsanitize=undefined
            -fsanitize=signed-integer-overflow
            -fsanitize=integer-divide-by-zero
            -fsanitize=float-divide-by-zero
            -fsanitize=alignment
            -fsanitize=bounds
            -fsanitize=return
            -fsanitize=shift
            -fsanitize=unreachable
            -fsanitize=vla-bound
            -fsanitize=null
            -fsanitize=pointer-compare
            -fsanitize=pointer-subtract
            -fsanitize=bool
            -fsanitize=enum
            -fsanitize=nonnull-attribute
            -fsanitize=float-cast-overflow
            -fno-sanitize-recover=all
            -fno-omit-frame-pointer
            -fsanitize-address-use-after-scope)

    # Convert the list of sanitizer in string to use as flag
    string(JOIN " " SANITIZER_FLAGS_STR ${SANITIZER_FLAGS})

    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -MMD -DNDEBUG -march=native")
    set(CMAKE_CXX_FLAGS_DEBUG "-g ${SANITIZER_FLAGS_STR}")
    set(LINKER_OPTIONS -flto -Wl,--no-a s-needed)
    set(LIB_OUTPUT_DIR ${PROJECT_SOURCE_DIR}/lib/${PROJECT}) # define our path for the lib
endif ()

# ---------------------------------------
# executable path
# ---------------------------------------
message(STATUS ${CMAKE_BUILD_TYPE} " mode" )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}")
# We indicate that we want an executable named PROJECT compiled from the files described by the SRCS and HEADERS variables.
add_executable(${PROJECT} ${SRCS_PROJECT_MAIN})

# include files
include_directories(
        ${HEADERS_PROJECT}
)
set_property(TARGET ${PROJECT} PROPERTY CXX_STANDARD 20)

# If we want to defined a library
# ---------------------------------------------------------
# Our F2Cmax Lib LIBRAIRY
# ---------------------------------------------------------
add_library(libF2Cmax STATIC ${SRCS_F2_Cmax}) # add the library to the sources
add_dependencies(${PROJECT} libF2Cmax) # PROJECT depends on LIB_NAME
set_property(TARGET libF2Cmax PROPERTY CXX_STANDARD 20)

#  Specify the output directory for the library file (.a)
set_target_properties(libF2Cmax PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${LIB_OUTPUT_DIR}
)

# --------------------------------------------------------
#                   ADD LIBRARIES
# --------------------------------------------------------


set(LIB_PATH "${PROJECT_SOURCE_DIR}/lib/")

# -----------------------------------------------------------
# ------- Inclusion of JSON library --------

message(STATUS "**************************************")
message(STATUS "**** ADDING external JSON library ****")
message(STATUS "**************************************")
add_subdirectory(${LIB_PATH}/json-3.11.3)


# specify target link of F2_Cmax

target_link_libraries(libF2Cmax
        nlohmann_json::nlohmann_json
)
# add link for library
target_link_libraries(${PROJECT}
        debug libF2Cmax
        optimized libF2Cmax
)
#need to specify the library's header file location.
target_include_directories(${PROJECT} PUBLIC
        ${PROJECT_SOURCE_DIR}/include) # main header

message(STATUS "*******************")
message(STATUS "**** TERMINATE ****")
message(STATUS "*******************")
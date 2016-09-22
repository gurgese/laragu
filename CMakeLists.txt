cmake_minimum_required(VERSION 3.5)
project(laragu)

set(CMAKE_CXX_STANDARD  14)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(./include)

set(SOURCE_FILES src/laragu.cpp)

add_executable(laragu ${SOURCE_FILES})

# Vienna RNA static library should be installed
find_library(RNA libRNA.a PATHS)

target_link_libraries(laragu ${RNA})
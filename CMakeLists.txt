cmake_minimum_required(VERSION 3.5.2)

# Project info
PROJECT(9815_Final_Project) 
SET(MAINFILE main)

file(GLOB_RECURSE SOURCES src/*.cpp)

# g++ config
SET(CMAKE_CXX_COMPILER "g++") 
SET(CMAKE_CXX_FLAGS "-std=c++14 -O3") 
SET(CMAKE_BUILD_TYPE Debug) 

# clang++ config
#SET(CMAKE_CXX_COMPILER "clang++") 
#SET(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ -Werror -Weverything -Wno-deprecated-declarations -Wno-disabled-macro-expansion -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-global-constructors -Wno-exit-time-destructors -Wno-missing-prototypes -Wno-padded -Wno-old-style-cast")
#SET(CMAKE_EXE_LINKER_FLAGS "-lc++ -lc++abi") 
#SET(CMAKE_BUILD_TYPE Debug)

ADD_EXECUTABLE("${MAINFILE}.out" ${SOURCES} "${MAINFILE}.cpp")


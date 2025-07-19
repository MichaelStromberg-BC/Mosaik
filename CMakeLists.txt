cmake_minimum_required(VERSION 3.10)

project(MOSAIK 
    VERSION 1.0.1388
    DESCRIPTION "Reference-guided assembler for next-generation sequencing data"
    LANGUAGES C CXX
)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wno-char-subscripts")
set(CMAKE_C_FLAGS_RELEASE "-O3 -w")

# Platform-specific definitions
if(UNIX)
    add_definitions(-D_FILE_OFFSET_BITS=64)
endif()

# Find required packages
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# Include directories for common source
include_directories(
    src/CommonSource/Config
    src/CommonSource/DataStructures 
    src/CommonSource/ExternalReadFormats
    src/CommonSource/MosaikReadFormat
    src/CommonSource/PairwiseAlignment
    src/CommonSource/Utilities
    src/CommonSource/AssemblyFormats
)

# Add subdirectories
add_subdirectory(src)
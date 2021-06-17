cmake_minimum_required(VERSION 3.10.2)
project (mipfinder VERSION 0.1
                   DESCRIPTION "Tool to find microProteins"
                   LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#Specify g++-10 compiler.
if (CMAKE_HOST_UNIX)
  set(CMAKE_CXX_COMPILER "g++-10")
endif()

#Set C++ standard to 20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#Set flags for RELEASE and BUILD 
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

#Set output dir for the binary
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

include_directories(include)
include_directories(lib)

add_executable(mipfinder
               ${CMAKE_CURRENT_SOURCE_DIR}/src/configuration.cc
               ${CMAKE_CURRENT_SOURCE_DIR}/src/easylogging++.cc
               ${CMAKE_CURRENT_SOURCE_DIR}/src/fasta.cc
               ${CMAKE_CURRENT_SOURCE_DIR}/src/file.cc
               ${CMAKE_CURRENT_SOURCE_DIR}/src/go.cc
               ${CMAKE_CURRENT_SOURCE_DIR}/src/helpers.cc
               ${CMAKE_CURRENT_SOURCE_DIR}/src/hmmer.cc
               ${CMAKE_CURRENT_SOURCE_DIR}/src/interpro.cc
               ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cc
               ${CMAKE_CURRENT_SOURCE_DIR}/src/mipfinder.cc
               ${CMAKE_CURRENT_SOURCE_DIR}/src/proteome.cc
               ${CMAKE_CURRENT_SOURCE_DIR}/src/protein.cc
               ${CMAKE_CURRENT_SOURCE_DIR}/src/printer.cc
)

#Set warnings for g++
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(mipfinder PRIVATE -Wall -Wextra -pedantic)
endif()

#Disable EasyLogging default log file
add_definitions("-DELPP_NO_DEFAULT_LOG_FILE")

#Enable the use of filesystem header with g++-8
target_link_libraries(mipfinder -lstdc++fs)
#add_library(mipfinder STATIC ${SOURCEFILES})




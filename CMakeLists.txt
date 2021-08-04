cmake_minimum_required(VERSION 3.14)

project (mipfinder VERSION 0.2
                   DESCRIPTION "Tool to find microProteins"
                   LANGUAGES CXX)

if (MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wpedantic -Wextra)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#Enable building tests
add_subdirectory(tests)

#Set flags for RELEASE and DEBUG builds 
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

#Disable EasyLogging default log file
add_compile_definitions("-DELPP_NO_DEFAULT_LOG_FILE")

#Set output dir for the binary
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(mipfinder
               ${CMAKE_CURRENT_SOURCE_DIR}/src/argparse.cc
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

target_compile_features(mipfinder PUBLIC cxx_std_20)
target_include_directories(mipfinder PUBLIC include)
target_include_directories(mipfinder PUBLIC lib)
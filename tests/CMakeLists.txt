cmake_minimum_required(VERSION 3.14)

project (mipfinder_tests VERSION 0.2
                         DESCRIPTION "Tool for testing mipfinder"
                         LANGUAGES CXX)

#For gtest (Copied from https://github.com/google/googletest/blob/master/googletest/README.md)
include(FetchContent)
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(FASTA_TEST_FILES
     ${CMAKE_SOURCE_DIR}/src/fasta.cc
     ${CMAKE_SOURCE_DIR}/src/file.cc
     ${CMAKE_CURRENT_SOURCE_DIR}/src/test_fasta.cc
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/bin)

add_executable(fasta_test
               ${FASTA_TEST_FILES}
)

target_compile_features(fasta_test PUBLIC cxx_std_20)
target_include_directories(fasta_test PUBLIC ${CMAKE_SOURCE_DIR}/include
                                             ${CMAKE_SOURCE_DIR}/lib
                                             ${CMAKE_CURRENT_SOURCE_DIR}/lib)
target_link_libraries(fasta_test gtest)
cmake_minimum_required(VERSION 3.10)
project(toga_cutils C)
set(CMAKE_C_STANDARD 99)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_executable(chain_score_filter chain_score_filter.c)
add_executable(chain_filter_by_id chain_filter_by_id.c)

add_library(chain_coords_converter_slib SHARED chain_coords_converter_slib.c)
add_library(extract_subchain_slib SHARED extract_subchain_slib.c)
add_library(chain_bst_lib SHARED chain_bst_lib.c)

set_target_properties(chain_coords_converter_slib PROPERTIES OUTPUT_NAME "chain_coords_converter_slib")
set_target_properties(extract_subchain_slib PROPERTIES OUTPUT_NAME "extract_subchain_slib")
set_target_properties(chain_bst_lib PROPERTIES OUTPUT_NAME "chain_bst_lib")

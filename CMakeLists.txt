cmake_minimum_required(VERSION 3.16)
project(pebble C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

add_library(pebble_core STATIC
    src/pebble.c
    src/pebble_io.c
)

target_include_directories(pebble_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_options(pebble_core PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)

add_executable(pebble src/main.c)
target_link_libraries(pebble PRIVATE pebble_core)
target_compile_options(pebble PRIVATE -Wall -Wextra -Wpedantic)

add_executable(pebble_tests tests/test_pebble.c)
target_link_libraries(pebble_tests PRIVATE pebble_core m)
target_compile_options(pebble_tests PRIVATE -Wall -Wextra -Wpedantic)

add_executable(pebble_io_tests tests/test_io.c)
target_link_libraries(pebble_io_tests PRIVATE pebble_core m)
target_compile_options(pebble_io_tests PRIVATE -Wall -Wextra -Wpedantic)

enable_testing()
add_test(NAME pebble_tests COMMAND pebble_tests)
add_test(NAME pebble_io_tests COMMAND pebble_io_tests)
set_tests_properties(pebble_io_tests PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

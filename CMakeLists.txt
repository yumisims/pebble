cmake_minimum_required(VERSION 3.16)
project(pebble C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Conda and other environments often prepend incompatible lib dirs to
# LD_LIBRARY_PATH, which makes dynamically linked HPC binaries segfault.
# Static Linux executables avoid that libc mismatch.
option(PEBBLE_STATIC "Statically link executables on Linux" ON)
option(PEBBLE_BIGWIG "Enable BigWig output (requires libBigWig and zlib)" OFF)

function(pebble_link_zlib target)
    if(PEBBLE_STATIC AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
        find_library(PEBBLE_ZLIB_LIBRARY NAMES libz.a z.a
            HINTS /usr/lib /usr/lib64 /usr/lib/x86_64-linux-gnu /lib/x86_64-linux-gnu)
        if(NOT PEBBLE_ZLIB_LIBRARY)
            message(FATAL_ERROR
                "PEBBLE_STATIC with PEBBLE_BIGWIG requires libz.a "
                "(on Debian/Ubuntu: apt install zlib1g-dev)")
        endif()
        target_link_libraries(${target} PUBLIC ${PEBBLE_ZLIB_LIBRARY})
    else()
        find_package(ZLIB REQUIRED)
        target_link_libraries(${target} PUBLIC ZLIB::ZLIB)
    endif()
endfunction()

set(PEBBLE_CORE_SOURCES
    src/pebble.c
    src/pebble_io.c
)

if(PEBBLE_BIGWIG)
    set(PEBBLE_LIBBIGWIG_DIR "${CMAKE_SOURCE_DIR}/third_party/libBigWig")
    if(NOT EXISTS "${PEBBLE_LIBBIGWIG_DIR}/bigWig.h")
        include(FetchContent)
        FetchContent_Declare(
            libBigWig
            GIT_REPOSITORY https://github.com/dpryan79/libBigWig.git
            GIT_TAG        0.4.8
        )
        FetchContent_GetProperties(libBigWig)
        if(NOT libBigWig_POPULATED)
            FetchContent_Populate(libBigWig)
        endif()
        set(PEBBLE_LIBBIGWIG_DIR "${libbigwig_SOURCE_DIR}")
    endif()

    add_library(libBigWig STATIC
        ${PEBBLE_LIBBIGWIG_DIR}/io.c
        ${PEBBLE_LIBBIGWIG_DIR}/bwValues.c
        ${PEBBLE_LIBBIGWIG_DIR}/bwRead.c
        ${PEBBLE_LIBBIGWIG_DIR}/bwStats.c
        ${PEBBLE_LIBBIGWIG_DIR}/bwWrite.c
    )
    target_include_directories(libBigWig PUBLIC "${PEBBLE_LIBBIGWIG_DIR}")
    target_compile_definitions(libBigWig PUBLIC NOCURL)
    target_compile_options(libBigWig PRIVATE -Wno-sign-compare)
    pebble_link_zlib(libBigWig)
    target_link_libraries(libBigWig PUBLIC m)

    list(APPEND PEBBLE_CORE_SOURCES src/pebble_bigwig.c)
else()
    list(APPEND PEBBLE_CORE_SOURCES src/pebble_bigwig_stub.c)
endif()

add_library(pebble_core STATIC ${PEBBLE_CORE_SOURCES})

target_include_directories(pebble_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_options(pebble_core PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)

if(PEBBLE_BIGWIG)
    target_compile_definitions(pebble_core PUBLIC PEBBLE_BIGWIG=1)
    target_link_libraries(pebble_core PUBLIC libBigWig)
endif()

function(pebble_apply_linux_static target)
    if(PEBBLE_STATIC AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_options(${target} PRIVATE -static)
    endif()
endfunction()

add_executable(pebble src/main.c)
target_link_libraries(pebble PRIVATE pebble_core)
target_compile_options(pebble PRIVATE -Wall -Wextra -Wpedantic)
if(PEBBLE_BIGWIG)
    target_compile_definitions(pebble PRIVATE PEBBLE_BIGWIG=1)
endif()
pebble_apply_linux_static(pebble)

add_executable(pebble_tests tests/test_pebble.c)
target_link_libraries(pebble_tests PRIVATE pebble_core m)
target_compile_options(pebble_tests PRIVATE -Wall -Wextra -Wpedantic)
pebble_apply_linux_static(pebble_tests)

add_executable(pebble_io_tests tests/test_io.c)
target_link_libraries(pebble_io_tests PRIVATE pebble_core m)
target_compile_options(pebble_io_tests PRIVATE -Wall -Wextra -Wpedantic)
if(PEBBLE_BIGWIG)
    target_compile_definitions(pebble_io_tests PRIVATE PEBBLE_BIGWIG=1)
    target_include_directories(pebble_io_tests PRIVATE "${PEBBLE_LIBBIGWIG_DIR}")
    target_link_libraries(pebble_io_tests PRIVATE libBigWig)
endif()
pebble_apply_linux_static(pebble_io_tests)

enable_testing()
add_test(NAME pebble_tests COMMAND pebble_tests)
add_test(NAME pebble_io_tests COMMAND pebble_io_tests)
set_tests_properties(pebble_io_tests PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

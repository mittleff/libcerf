cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(PreventInSourceBuilds)

if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build as shared library" ON)
endif()
if(NOT DEFINED LIB_MAN)
    option(LIB_MAN "Build libcerf man pages" ON)
endif()
if(NOT DEFINED LIB_RUN)
    option(LIB_RUN "Build executables for command-line computation" ON)
endif()
if(NOT DEFINED CERF_INTROSPECT)
    option(CERF_INTROSPECT "Provide variables cerf_algorithm and cerf_nofterms to inform about used computation" OFF)
endif()
if(NOT DEFINED CERF_CPP)
    option(CERF_CPP "Build libcerf as native C++" OFF)
endif()
if(NOT DEFINED CERF_IEEE754)
    option(CERF_IEEE754 "Assume IEEE754 floating-point layout" ON)
    # This switch allows us to fall back to frexp from libc instead of using our own frexp2.
    # It is unclear though whether the remainder of the code will work on non-IEEE754 machines.
endif()

if(CERF_CPP)
    message("Build C++ library libcerfcpp")
    project(cerfcpp CXX)
    set(CMAKE_CXX_STANDARD 11)
    add_compile_definitions(CERF_AS_CPP)
else()
    message("Build C library libcerf")
    project(cerf C)
    set(CMAKE_C_STANDARD 11)
endif()

set(CERF_SOVERSION                 2) # API version
set(CERF_VERSION ${CERF_SOVERSION}.5) # minor version

include(GNUInstallDirs)

SET(CERF_COMPILE_OPTIONS "" CACHE
    STRING "User-supplied compiler options, will be appended to hard-coded ones. Not for Windows.")

if (CERF_INTROSPECT)
    add_compile_definitions(CERF_INTROSPECT)
endif()

if (NOT CERF_IEEE754)
    add_compile_definitions(CERF_NO_IEEE754)
endif()

if(MSVC)
    if (NOT CERF_CPP)
        message(FATAL_ERROR "Under MSVC, only CERF_CPP=ON is supported")
    endif()
    add_compile_options("/MP") # parallel compilation
    add_compile_options("$<$<CONFIG:Debug>:/MTd>") # multithreaded, debug, required for post-build under vcpkg
    add_compile_options("/wd4244")
    # /wd4018 /wd4068 /wd4101 /wd4244 /wd4267 /wd4305 /wd4715") # warnings off
    add_compile_options("-D_CRT_SECURE_NO_WARNINGS -D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING")
    set(CTEST_CONFIGURATION_TYPE "${JOB_BUILD_CONFIGURATION}")
    if(BUILD_SHARED_LIBS)
        set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    endif()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

elseif(CERF_COMPILE_OPTIONS)
    add_compile_options(${CERF_COMPILE_OPTIONS})

else()
    option(PEDANTIC "Compile with pedantic warnings" ON)
    option(WERROR "Treat warnings as errors" OFF)
    if(PEDANTIC)
        add_compile_options(-pedantic -Wall)
    endif()
    if(WERROR)
        add_compile_options(-Werror)
    endif()

    add_compile_options(-Wno-sign-compare -fno-omit-frame-pointer)
    if (CMAKE_BUILD_TYPE MATCHES Debug)
        add_compile_options(-g)
    else()
        add_compile_options(-O3)
    endif()
endif()

include(CTest)

configure_file(var/libcerf.pc.in libcerf.pc @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libcerf.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

add_subdirectory(lib)
add_subdirectory(test)
if(LIB_RUN)
    add_subdirectory(run)
endif()
if(LIB_MAN)
    add_subdirectory(man)
endif()

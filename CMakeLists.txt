project(UVAT)
CMAKE_MINIMUM_REQUIRED(VERSION 3.11)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Print the build type

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

### Configure the compiler
# This is a basic, decent setup that should do something sane on most compilers

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++17 -Wall -Wextra ")
  SET(DISABLED_WARNINGS " -Wdeprecated-copy -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-narrowing -Wno-sign-compare -Wno-sequence-point")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3 -fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Setting gcc-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wdeprecated-copy -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()

  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-O3  -DNDEBUG")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()


# == Deps
add_subdirectory(deps/geometry-central)
add_subdirectory(deps/polyscope)
add_subdirectory(deps/libigl)
SET(BUILD_EXAMPLES OFF CACHE BOOL "DGtal examples")
add_subdirectory(deps/DGtal)
file(GLOB SRCFILES *.cpp)

##DGtal Specific settings
include("${PROJECT_BINARY_DIR}/deps/DGtal/DGtalConfig.cmake")
include_directories("${DGTAL_INCLUDE_DIRS}")

FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  ADD_DEFINITIONS(${BOOST_DEFINITIONS} -DBOOST_ALL_NO_LIB)
  include_directories( ${Boost_INCLUDE_DIRS} )
endif( Boost_FOUND )

# == Build our project stuff

add_executable(at-uv "src/at-uv.cpp" ${SRCFILES})
target_include_directories(at-uv PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(at-uv PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/")
target_sources(at-uv PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/UV.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/AT.cpp")
target_link_libraries(at-uv geometry-central polyscope DGtal ${DGTAL_LIBRARIES}   ${Boost_PROGRAM_OPTIONS_LIBRARY} Boost::program_options igl::core)


# To mark ".ih" as _headers in Xcode
FILE(GLOB_RECURSE UV_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/src/*h)
add_custom_target(_headers SOURCES ${UV_HEADERS})

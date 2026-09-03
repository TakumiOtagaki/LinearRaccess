cmake_minimum_required(VERSION 3.20)

project(LinearRaccess VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LINEARRACCESS_WITH_RACCESS
       "Build the optional adapter for a separately obtained Raccess source tree"
       OFF)
set(RACCESS_INCLUDE_DIR "" CACHE PATH
    "Path containing raccess/ and util/ headers (normally <raccess>/src)")

set(LINEARRACCESS_SOURCES
    src/FileReader.cpp
    src/accessibility.cpp
    src/api.cpp
    src/beam_inside_outside.cpp
    src/debug.cpp
    src/debugging.cpp
    src/energy_linearcapr.cpp
    src/io.cpp
    src/legacy_energy_param.cpp
    src/lincapr.cpp
)

if(LINEARRACCESS_WITH_RACCESS)
  if(NOT EXISTS "${RACCESS_INCLUDE_DIR}/raccess/energy_model_api.hpp")
    message(FATAL_ERROR
      "LINEARRACCESS_WITH_RACCESS=ON requires RACCESS_INCLUDE_DIR to point "
      "to an authorized Raccess include tree (normally <raccess>/src)")
  endif()
  list(APPEND LINEARRACCESS_SOURCES src/energy_raccess.cpp)
endif()

add_library(linearraccess STATIC ${LINEARRACCESS_SOURCES})
target_include_directories(linearraccess
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  PRIVATE
    src
)
target_compile_features(linearraccess PUBLIC cxx_std_17)

if(LINEARRACCESS_WITH_RACCESS)
  target_include_directories(linearraccess PRIVATE "${RACCESS_INCLUDE_DIR}")
  target_compile_definitions(linearraccess PUBLIC LINEARRACCESS_WITH_RACCESS=1)
  # Older Raccess releases use standard-library adapters removed in C++17.
  target_compile_definitions(linearraccess PRIVATE
    _LIBCPP_ENABLE_CXX17_REMOVED_BINDERS
    _LIBCPP_ENABLE_CXX17_REMOVED_RANDOM_SHUFFLE
    _LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION
  )
endif()

add_executable(LinRacc src/linracc_main.cpp)
target_include_directories(LinRacc PRIVATE src)
target_link_libraries(LinRacc PRIVATE linearraccess)

include(CTest)
if(BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)

  add_executable(linearraccess_core_test tests/core_test.cpp)
  target_include_directories(linearraccess_core_test PRIVATE src)
  target_link_libraries(linearraccess_core_test PRIVATE linearraccess)
  add_test(NAME core COMMAND linearraccess_core_test)

  set(LINEARRACCESS_CLI_TEST_ARGS
      --executable "$<TARGET_FILE:LinRacc>")
  if(LINEARRACCESS_WITH_RACCESS)
    list(APPEND LINEARRACCESS_CLI_TEST_ARGS --with-raccess)
  endif()
  add_test(
    NAME cli
    COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_smoke_test.py"
            ${LINEARRACCESS_CLI_TEST_ARGS}
  )
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
install(TARGETS linearraccess LinRacc
        EXPORT LinearRaccessTargets
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(DIRECTORY include/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(EXPORT LinearRaccessTargets
        FILE LinearRaccessTargets.cmake
        NAMESPACE LinearRaccess::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/LinearRaccess")
configure_package_config_file(
  cmake/LinearRaccessConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/LinearRaccessConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/LinearRaccess"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/LinearRaccessConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/LinearRaccessConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/LinearRaccessConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/LinearRaccess")

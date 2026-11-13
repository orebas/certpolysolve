cmake_minimum_required(VERSION 3.20)
project(CertPolySolve VERSION 0.5.1 LANGUAGES C CXX)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)
include(CTest)
if(BUILD_TESTING AND CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_GENERATOR MATCHES "^(Ninja|Unix Makefiles)$")
    # Set before targets are created; the compiled module check requires their commands.
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()
# Build-tree executables find libcertpolysolve through $ORIGIN, not an absolute
# build path, so a frozen copy of workers plus library is self-contained and
# a later rebuild of this tree cannot change what a frozen worker loads.
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)
include(CMakePackageConfigHelpers)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/Numerics.cmake")
find_package(FLINT 3.6 REQUIRED)

option(CERTPOLYSOLVE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(CERTPOLYSOLVE_TSAN "Enable ThreadSanitizer (the tsan preset; meaningful inside FLINT only with a FLINT built with it)" OFF)
option(CERTPOLYSOLVE_CLANG_TIDY "Check move/borrow lifetime rules during C++ builds" OFF)
if(CERTPOLYSOLVE_CLANG_TIDY)
    find_program(CERTPOLYSOLVE_CLANG_TIDY_EXECUTABLE NAMES clang-tidy REQUIRED)
    set(CMAKE_CXX_CLANG_TIDY "${CERTPOLYSOLVE_CLANG_TIDY_EXECUTABLE}")
endif()
option(CERTPOLYSOLVE_BUILD_STRESS "Build replayable independent exact-checking campaigns" OFF)
# The 128-bit integers and checked arithmetic of src/core/wide_int.hpp as a compiler without __int128 and the GCC
# overflow built-ins (MSVC) gets them, forced on GCC and Clang so that the suite and the gate worker can run on
# them here. Every target of this directory gets the definition: a header with a 128-bit member
# (node_products.hpp) must mean the same type in the library and in every test that includes it.
option(CERTPOLYSOLVE_PORTABLE_INT128 "Use the portable 128-bit integers and checked arithmetic of src/core/wide_int.hpp instead of the compiler's" OFF)
if(CERTPOLYSOLVE_PORTABLE_INT128)
    add_compile_definitions(CERTPOLYSOLVE_PORTABLE_INT128)
endif()
find_package(Threads REQUIRED)
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules.cmake")
add_library(certpolysolve SHARED ${certpolysolve_api_sources})
add_library(CertPolySolve::certpolysolve ALIAS certpolysolve)
target_compile_features(certpolysolve PUBLIC cxx_std_20)
target_include_directories(certpolysolve PRIVATE src)
target_include_directories(certpolysolve PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(certpolysolve PRIVATE certpolysolve_core FLINT::FLINT Threads::Threads)
target_compile_definitions(certpolysolve PRIVATE CERTPOLYSOLVE_BUILDING_LIBRARY)
set_target_properties(certpolysolve PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 0
    CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN YES)
if(MSVC)
    target_compile_options(certpolysolve PRIVATE /W4 /permissive-)
    certpolysolve_strict_floating_point(certpolysolve)
else()
    # One source gives the same doubles on every CPU: a product and a sum are fused only where the source
    # writes std::fma (Clang's default fuses them within an expression wherever the target has the
    # instruction: AArch64, or x86-64 built for a recent -march).
    target_compile_options(certpolysolve PRIVATE -Wall -Wextra -Wpedantic -Wconversion)
    certpolysolve_strict_floating_point(certpolysolve)
endif()
if(CERTPOLYSOLVE_SANITIZERS)
    if(MSVC)
        message(FATAL_ERROR "This sanitizer preset currently requires GCC or Clang")
    endif()
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined)
    target_compile_options(certpolysolve PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(certpolysolve PRIVATE -fsanitize=address,undefined)
endif()
# ThreadSanitizer (the tsan preset, dev/docs/thread-safety-audit-20261018.md F1). It sees the memory accesses of
# instrumented code only: against an uninstrumented FLINT it checks CertPolySolve's own code and nothing inside FLINT,
# so the preset takes FLINT from FLINT_ROOT, a FLINT built with -fsanitize=thread.
if(CERTPOLYSOLVE_TSAN)
    if(MSVC OR CERTPOLYSOLVE_SANITIZERS)
        message(FATAL_ERROR "ThreadSanitizer needs GCC or Clang and excludes CERTPOLYSOLVE_SANITIZERS")
    endif()
    if(NOT FLINT_ROOT AND NOT DEFINED ENV{FLINT_ROOT})
        message(WARNING "CERTPOLYSOLVE_TSAN with the system FLINT: races inside FLINT are not seen; set FLINT_ROOT to a FLINT built with -fsanitize=thread")
    endif()
    add_compile_options(-fsanitize=thread -fno-omit-frame-pointer)
    add_link_options(-fsanitize=thread)
    target_compile_options(certpolysolve PRIVATE -fsanitize=thread -fno-omit-frame-pointer)
    target_link_options(certpolysolve PRIVATE -fsanitize=thread)
endif()

add_executable(certpolysolve_cli examples/solve.cpp)
target_link_libraries(certpolysolve_cli PRIVATE CertPolySolve::certpolysolve)

if(BUILD_TESTING OR CERTPOLYSOLVE_BUILD_STRESS)
    # Header-only Boost is test-only: an exact rational reference independent
    # of FLINT's root-finding and enclosure arithmetic.
    find_path(CERTPOLYSOLVE_BOOST_INCLUDE_DIR boost/multiprecision/cpp_int.hpp REQUIRED)
endif()
if(BUILD_TESTING)
    include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/Tests.cmake")
endif()


if(CERTPOLYSOLVE_BUILD_STRESS)
    add_executable(certpolysolve_stress tests/stress_core.cpp)
    target_include_directories(certpolysolve_stress PRIVATE src "${CERTPOLYSOLVE_BOOST_INCLUDE_DIR}")
    target_link_libraries(certpolysolve_stress PRIVATE certpolysolve_core CertPolySolve::certpolysolve FLINT::FLINT)
    if(BUILD_TESTING)
        add_test(NAME stress_smoke COMMAND certpolysolve_stress --cases 20)
        set_tests_properties(stress_smoke PROPERTIES TIMEOUT 120)
        # Each explicit method, with exact transforms of every planted-root case
        # (rotation by (3+4i)/5, shift by a 64-bit sqrt(2), conjugation, pi multiple).
        # --method takes the cps_method value: 1 classical, 2 secular, 3 piecewise.
        add_test(NAME stress_transforms_classical COMMAND certpolysolve_stress --cases 8 --method 1 --transforms 1)
        add_test(NAME stress_transforms_secular COMMAND certpolysolve_stress --cases 8 --method 2 --transforms 1)
        add_test(NAME stress_transforms_piecewise COMMAND certpolysolve_stress --cases 8 --method 3 --transforms 1)
        set_tests_properties(stress_transforms_classical stress_transforms_secular stress_transforms_piecewise
                             PROPERTIES TIMEOUT 300)
    endif()
endif()

install(TARGETS certpolysolve EXPORT CertPolySolveTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
# certpolysolve_testing.h is the tests' header (a stop at an exact cancellation poll), not part of the installed API.
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} PATTERN "certpolysolve_testing.h" EXCLUDE)
install(EXPORT CertPolySolveTargets NAMESPACE CertPolySolve:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CertPolySolve)
configure_package_config_file(cmake/CertPolySolveConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/CertPolySolveConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CertPolySolve)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/CertPolySolveConfigVersion.cmake"
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMinorVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/CertPolySolveConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/CertPolySolveConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CertPolySolve)
install(FILES COPYING COPYING.LESSER DESTINATION ${CMAKE_INSTALL_DATADIR}/certpolysolve)

install(DIRECTORY licenses/ DESTINATION ${CMAKE_INSTALL_DATADIR}/certpolysolve/licenses)

# Test purpose, independent of whether a research/benchmark build is enabled.
# Labels allow a short ordinary edit cycle without removing mathematical tests.
get_property(certpolysolve_tests DIRECTORY PROPERTY TESTS)
foreach(test IN LISTS certpolysolve_tests)
    if(test MATCHES "^tooling_")
        set_property(TEST ${test} PROPERTY LABELS "tooling;quick")
    elseif(test MATCHES "^(portfolio|local_contraction|campaign_.*|default_consolidation_tools|benchmark_.*|corpus_.*|mpsolve_battery|mps_certificate|acceptance_gate)$")
        set_property(TEST ${test} APPEND PROPERTY LABELS research)
    elseif(test MATCHES "^stress")
        set_property(TEST ${test} APPEND PROPERTY LABELS stress)
    else()
        set_property(TEST ${test} APPEND PROPERTY LABELS correctness)
    endif()
    if(test MATCHES "^(core|c_abi|methods|solve_state|allocation_faults|secular_cli|resources|no_backend_solvers|guide_first_solve_c|guide_first_solve_cpp|examples_quoted)$")
        set_property(TEST ${test} APPEND PROPERTY LABELS api)
    endif()
    if(test MATCHES "^(core|c_abi|methods|allocation_faults|approximation|ownership|root_index|route_invariants|no_backend_solvers|doc_links|seed_equivariance|equivalence_properties|guide_first_solve_c|guide_first_solve_cpp|examples_quoted)$")
        set_property(TEST ${test} APPEND PROPERTY LABELS quick)
    endif()
endforeach()

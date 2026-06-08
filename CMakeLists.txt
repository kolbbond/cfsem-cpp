# cfsem-cpp — C-ABI shim library over the cfsem Rust crate.
#
# Layout:
#   cfsem-py/    Rust crate `cfsem` (git submodule, upstream, untouched)
#   cfsem-capi/  Rust crate `cfsem-capi` -> libcfsem_capi.a (extern "C" wrappers)
#   include/     cfsem.h (C ABI) + cfsem.hpp (C++ namespace wrapper)
#
# Targets:
#   cfsem_cpp   INTERFACE library bundling the Rust static lib + headers.
#               This is what downstream consumers (e.g. goose) link against:
#                   target_link_libraries(goose PRIVATE cfsem_cpp)
#   cfsem-cpp   smoke-test executable (validates the library builds/links/runs).
cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Emit compile_commands.json so clangd / IDEs resolve include/ correctly.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(
  cfsem-cpp
  VERSION 0.1
  LANGUAGES C CXX)

# Leverage Corrosion to build the Rust crates through cmake/cargo.
include(FetchContent)
FetchContent_Declare(
  Corrosion
  GIT_REPOSITORY https://github.com/corrosion-rs/corrosion.git
  GIT_TAG v0.5)
FetchContent_MakeAvailable(Corrosion)

# Build the C-ABI shim crate. Corrosion creates an imported target named after
# the crate's lib name (hyphens normalized to underscores): `cfsem_capi`. It
# links libcfsem_capi.a plus the system libraries the Rust staticlib requires.
corrosion_import_crate(MANIFEST_PATH cfsem-capi/Cargo.toml)

# Public library target consumers link against: Rust static lib + headers.
add_library(cfsem_cpp INTERFACE)
target_include_directories(cfsem_cpp
                           INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cfsem_cpp INTERFACE cfsem_capi)

# Smoke-test executable. Not the deliverable — proves the library links and runs.
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE cfsem_cpp)

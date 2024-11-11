cmake_minimum_required(VERSION 3.12)
include(FetchContent)
project(run_systematics)

add_compile_options(-Wall -Werror -Wno-error=maybe-uninitialized -Ofast)

# Find packages
find_package(ROOT REQUIRED)
find_package(sbnanaobj)

FetchContent_Declare(
    tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG v3.4.0
)
FetchContent_MakeAvailable(tomlplusplus)

set(SBNANAOBJ_INCLUDE_DIRS "$ENV{MRB_SOURCE}/sbnanaobj/")

file(GLOB SYSINC "include/*.h")

# Library for configuration details
add_library(configuration SHARED src/configuration.cc)
target_link_libraries(configuration PRIVATE tomlplusplus::tomlplusplus)
target_include_directories(configuration PRIVATE include/ tomlplusplus/include)

# Add the executable target
add_executable(run_systematics src/main.cc ${SYSINC})

# Link the ROOT libraries to the target
target_link_libraries(run_systematics ${ROOT_LIBRARIES} ${sbnanaobj_LIBRARY_DIRS}/libsbnanaobj_StandardRecord.so tomlplusplus::tomlplusplus configuration)

# Include the ROOT headers
include_directories(${ROOT_INCLUDE_DIRS} ${SBNANAOBJ_INCLUDE_DIRS} include/ tomlplusplus/include)

# Add ROOT definitions
add_definitions(${ROOT_CXX_FLAGS})
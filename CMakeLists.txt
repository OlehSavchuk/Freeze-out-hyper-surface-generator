cmake_minimum_required(VERSION 3.10)
project(smash_sampling LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(hypersurface_builder
    src/main.cpp
    src/hypersurface_builder.cpp
    src/cornelius.cpp
)
target_include_directories(hypersurface_builder PRIVATE include)

target_link_libraries(hypersurface_builder)

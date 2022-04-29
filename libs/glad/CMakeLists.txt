project(glad
    VERSION 1.0.0
    LANGUAGES C
)

add_library(glad)
target_sources(glad PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/glad.c")

target_include_directories(glad PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include" )
add_library(glad::glad ALIAS glad)
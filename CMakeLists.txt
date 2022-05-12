cmake_minimum_required(VERSION 3.23.1)
project(WorldGenerator
    VERSION 1.0.0
    LANGUAGES CXX)

set(EXECUTABLE_NAME "WorldGenerator")

include(cmake/cpm.cmake)

add_subdirectory( "${CMAKE_CURRENT_SOURCE_DIR}/libs/glad" )

CPMAddPackage("gh:g-truc/glm#0.9.9.1")
CPMAddPackage("gh:glfw/glfw#3.3.4")
CPMAddPackage(
  NAME imgui
  GIT_REPOSITORY https://github.com/ocornut/imgui.git
  VERSION 1.87
  DOWNLOAD_ONLY YES
)

if (imgui_ADDED)
  FILE(GLOB imgui_sources ${imgui_SOURCE_DIR}/*.cpp)
  list(APPEND imgui_sources "${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp;${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp")
  message(${imgui_sources})
  add_library(imgui)
  target_sources(imgui PRIVATE ${imgui_sources})
  target_link_libraries(imgui
      PUBLIC glad::glad
      PUBLIC glfw)
  target_include_directories(imgui
    PUBLIC ${imgui_SOURCE_DIR}
    PUBLIC glad::glad
    PUBLIC glfw
  )
endif()

# GLFW options
set(GLFW_BUILD_EXAMPLES OFF CACHE INTERNAL "Build the GLFW example programs")
set(GLFW_BUILD_TESTS OFF CACHE INTERNAL "Build the GLFW test programs")
set(GLFW_BUILD_DOCS OFF CACHE INTERNAL "Build the GLFW documentation")
set(GLFW_INSTALL OFF CACHE INTERNAL "Generate installation target")

set(GLM_TEST_ENABLE OFF CACHE INTERNAL "Build the GLM test programs")

add_executable(${EXECUTABLE_NAME})
target_compile_features(${EXECUTABLE_NAME} PUBLIC cxx_std_20)
target_sources(${EXECUTABLE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
target_sources(${EXECUTABLE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/perlin.cpp)
target_sources(${EXECUTABLE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/stb.cpp)

target_compile_definitions(${EXECUTABLE_NAME} PUBLIC -DGLFW_INCLUDE_NONE -DGLM_ENABLE_EXPERIMENTAL)

target_link_libraries(${EXECUTABLE_NAME} 
    PUBLIC glad::glad 
    PUBLIC glm 
    PUBLIC glfw
    PUBLIC imgui)

target_include_directories(${EXECUTABLE_NAME} 
    PUBLIC glad::glad 
    PUBLIC glm 
    PUBLIC glfw
    PUBLIC imgui)

target_include_directories(${EXECUTABLE_NAME} 
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
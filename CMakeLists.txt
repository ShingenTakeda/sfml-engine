cmake_minimum_required(VERSION 3.16.0)

project(sfml_engine)

#SFML
find_package(SFML 2.5 REQUIRED network audio graphics window system)
if(SFML_FOUND)
  include_directories(${SFML_INCLUDE_DIR})
endif()

#IMGUI_SFML
find_package(ImGui-SFML REQUIRED)
include_directories(${IMGUI_SFML_INCLUDE_DIR})

#LUA
include_directories(${PROJECT_SOURCE_DIR}/extern/lua-5.4.4/src/)
file(GLOB_RECURSE LUA_SRC "${PROJECT_SOURCE_DIR}/extern/lua-5.4.4/src/*.c" "${PROJECT_SOURCE_DIR}/extern/lua-5.4.4/src/*.h" "${PROJECT_SOURCE_DIR}/extern/lua-5.4.4/src/*.hpp")
add_library(lua STATIC ${LUA_SRC})

#SQUIRREL
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/squirrel3/")
include_directories(${PROJECT_SOURCE_DIR}/extern/squirrel3/include)

#WREN
include_directories(${PROJECT_SOURCE_DIR}/extern/wren/src/include/)
file(GLOB_RECURSE WREN_SRC "${PROJECT_SOURCE_DIR}/extern/wren/src/vm/*.c" "${PROJECT_SOURCE_DIR}/extern/wren/src/vm/*.hpp")
add_library(wren ${WREN_SRC})

#ENGINE
include_directories(${PROJECT_SOURCE_DIR}/engine/src/)
file(GLOB_RECURSE ENGINE_SRC "${PROJECT_SOURCE_DIR}/engine//src/*.cpp" "${PROJECT_SOURCE_DIR}/engine//src/*.hpp")
add_library(sfml_engine STATIC ${ENGINE_SRC})
target_link_libraries(sfml_engine  sfml-graphics sfml-audio sfml-window sfml-network sfml-system lua)

#EXAMPLES
file(GLOB_RECURSE SIMPLE_SRC "${PROJECT_SOURCE_DIR}/examples/simple/*.cpp" "${PROJECT_SOURCE_DIR}/examples/simple/*.hpp")
add_executable(simple ${SIMPLE_SRC})
target_link_libraries(simple sfml_engine)

file(GLOB_RECURSE IMGUITEST_SRC "${PROJECT_SOURCE_DIR}/examples/imgui/*.cpp")
add_executable(imgui ${IMGUITEST_SRC})
target_link_libraries(imgui PRIVATE sfml_engine  ImGui-SFML::ImGui-SFML)

#add_executable(lua_example ${PROJECT_SOURCE_DIR}/examples/lua/main.cpp)
#target_link_libraries(lua_example lua)
#install(TARGETS sfml_engine RUNTIME DESTINATION bin)

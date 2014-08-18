project(hump_changing_height)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp lumped_projection.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

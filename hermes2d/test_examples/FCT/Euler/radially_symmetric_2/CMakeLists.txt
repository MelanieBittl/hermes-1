project(radially_symmetric_2)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp lumped_projection.cpp euler_util.cpp initial_condition.cpp boundary_condition.cpp)



set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

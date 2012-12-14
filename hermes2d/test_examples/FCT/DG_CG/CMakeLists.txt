project(DG_CG)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp lumped_projection.cpp hp_adapt.cpp prev_solution.cpp l2_semi_cg_space.cpp )

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

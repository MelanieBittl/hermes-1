project(DG_CG_stationary_convection_diffusion)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp )

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

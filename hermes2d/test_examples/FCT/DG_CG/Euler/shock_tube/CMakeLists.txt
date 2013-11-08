project(shock_tube_CGDG)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp euler_util.cpp initial_condition.cpp boundary_condition.cpp interface.cpp numerical_flux.cpp lumped_projection.cpp hp_adapt.cpp boundary_convection.cpp)
#include (${hermes2d_SOURCE_DIR}/CMake.common)
#set_common_target_properties(${PROJECT_NAME})


set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

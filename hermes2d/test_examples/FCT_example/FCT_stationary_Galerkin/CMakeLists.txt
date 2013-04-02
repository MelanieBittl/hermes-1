project(FCT_stationary_Galerkin)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp )
#include (${hermes2d_SOURCE_DIR}/CMake.common)
#set_common_target_properties(${PROJECT_NAME})


set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

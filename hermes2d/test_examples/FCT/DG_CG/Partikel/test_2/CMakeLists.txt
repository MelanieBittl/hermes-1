project(partikel_2)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp euler_util.cpp initial_condition.cpp interface.cpp numerical_flux.cpp euler_flux.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

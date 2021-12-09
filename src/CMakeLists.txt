find_package(GSL REQUIRED)

set(APP_NAME fast_sampen)
set(LIB_NAME sampen)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(CPP_LIST utils.cpp sample_entropy_calculator.cpp random_sampler.cpp kdtree.cpp sampen_entropy_caculator_kd.cpp)
include_directories(${CMAKE_SOURCE_DIR}/include)

add_library(${LIB_NAME} ${CPP_LIST})
target_link_libraries(${LIB_NAME} GSL::gsl GSL::gslcblas)

add_executable(${APP_NAME} fast_sampen.cpp)
target_link_libraries(${APP_NAME} ${LIB_NAME})

add_executable(experiment_n0n1 experiment_n0n1.cpp kdtree.cpp)
target_link_libraries(experiment_n0n1 ${LIB_NAME})

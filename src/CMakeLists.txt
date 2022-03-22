find_package(GSL REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(CPP_LIST utils.cpp 
    sample_entropy_calculator.cpp
    random_sampler.cpp
    kdtree.cpp
    sampen_entropy_caculator_kd.cpp
    sample_entropy_calculator_direct.cpp)
include_directories(${CMAKE_SOURCE_DIR}/include)

add_library(${LIB_NAME} ${CPP_LIST})
target_link_libraries(${LIB_NAME} GSL::gsl GSL::gslcblas)

add_executable(${APP_NAME} fast_sampen.cpp)
target_link_libraries(${APP_NAME} ${LIB_NAME})

add_executable(experiment_n0n1 experiment_n0n1.cpp kdtree.cpp)
target_link_libraries(experiment_n0n1 ${LIB_NAME})

find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGICKPP REQUIRED Magick++)
add_executable(sampen2d sampen2d.cpp)
message(INFO "MAGICKPP_INCLUDE_DIRS: ${MAGICKPP_INCLUDE_DIRS}")
target_include_directories(sampen2d PUBLIC ${MAGICKPP_INCLUDE_DIRS})
target_link_libraries(sampen2d sampen ${MAGICKPP_LINK_LIBRARIES})

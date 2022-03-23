find_package(GSL REQUIRED)
include_directories(${CMAKE_SOURCE_DIR}/include)

set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(CPP_LIST utils.cpp 
    sample_entropy_calculator.cpp
    random_sampler.cpp
    kdtree.cpp
    sampen_entropy_caculator_kd.cpp
    sample_entropy_calculator_direct.cpp)

set(PUBLIC_HEADERS global_defs.h;utils.h;kdtree.h;kdpoint.h;sample_entropy_calculator.h;sample_entropy_calculator_kd.h;sample_entropy_calculator_direct.h;sample_entropy_calculator2d.h;random_sampler.h)
add_library(${LIB_NAME} SHARED ${CPP_LIST})
target_link_libraries(${LIB_NAME} GSL::gsl GSL::gslcblas)
target_include_directories(${LIB_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_executable(${APP_NAME} fast_sampen.cpp)
target_link_libraries(${APP_NAME} ${LIB_NAME})

add_executable(experiment_n0n1 experiment_n0n1.cpp)
target_link_libraries(experiment_n0n1 ${LIB_NAME})

find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGICKPP REQUIRED Magick++)
add_executable(sampen2d sampen2d.cpp)
message(INFO "MAGICKPP_INCLUDE_DIRS: ${MAGICKPP_INCLUDE_DIRS}")
target_include_directories(sampen2d PUBLIC ${MAGICKPP_INCLUDE_DIRS})
target_link_libraries(sampen2d sampen ${MAGICKPP_LINK_LIBRARIES})

install(TARGETS ${LIB_NAME} ${APP_NAME} sampen2d experiment_n0n1
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)

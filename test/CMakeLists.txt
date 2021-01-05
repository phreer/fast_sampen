find_package(GSL REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/include)
add_executable(test_swr test_swr.cpp ../src/random_sampler.cpp)
target_link_libraries(test_swr GSL::gsl GSL::gslcblas)

add_executable(test_grid test_grid.cpp ../src/random_sampler.cpp)
target_link_libraries(test_grid GSL::gsl GSL::gslcblas)
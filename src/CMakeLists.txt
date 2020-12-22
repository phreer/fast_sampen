find_package(GSL REQUIRED)

set(CPP_LIST utils.cpp sample_entropy_calculator.cpp random_sampler.cpp) 
include_directories(${CMAKE_SOURCE_DIR}/include)
add_executable(kdtree_sample kdtree_sample.cpp ${CPP_LIST})
target_link_libraries(kdtree_sample GSL::gsl GSL::gslcblas)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3")

set(CPP_LIST utils.cpp sample_entropy_calculator.cpp random_sampler.cpp) 
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3")

include_directories(${CMAKE_SOURCE_DIR}/include)
find_package(GSL REQUIRED)
add_executable(kdtree_sample kdtree_sample.cpp ${CPP_LIST})
target_link_libraries(kdtree_sample GSL::gsl GSL::gslcblas)

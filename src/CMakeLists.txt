set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3")

set(CPP_LIST utils.cpp sample_entropy_calculator.cpp) 
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3")

include_directories($ENV{HOME}/local/include)
link_directories($ENV{HOME}/local/lib)
add_executable(kdtree_mddc kdtree_mddc.cpp ${CPP_LIST})
add_executable(kdtree_sample kdtree_sample.cpp ${CPP_LIST} random_sampler.cpp)
target_link_libraries(kdtree_sample gsl gslcblas)
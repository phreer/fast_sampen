set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3")

set(CPP_LIST kdtree_mddc.cpp utils.cpp sample_entropy_calculator.cpp) 
set(HEADER_LIST kdpoint.h kdtree.h sample_entropy_calculator.h sample_entropy_calculator_kd.h sample_entropy_calculator_direct.h utils.h)
add_executable(kdtree_mddc ${CPP_LIST} ${HEADER_LIST})
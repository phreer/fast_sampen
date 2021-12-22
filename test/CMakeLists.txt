find_package(GSL REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/include)
add_executable(test_swr test_swr.cpp)
target_link_libraries(test_swr sampen)

add_executable(test_grid test_grid.cpp)
target_link_libraries(test_grid sampen)
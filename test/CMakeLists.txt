find_package(Boost 1.76 REQUIRED COMPONENTS unit_test_framework)


add_executable(test_ggq test_ggq.cpp)


target_link_libraries(test_ggq Boost::unit_test_framework)
#Required packages
find_package(Boost 1.76 REQUIRED COMPONENTS unit_test_framework)

#Excutables
add_executable(test_ggq test_ggq.cpp)
add_executable(test_discretizer test_discretizer.cpp)

#Link library
target_link_libraries(test_ggq ggq_library)
target_link_libraries(test_discretizer ggq_library)

#Link testing framework
target_link_libraries(test_ggq Boost::unit_test_framework)
target_link_libraries(test_discretizer Boost::unit_test_framework)
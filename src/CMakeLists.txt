set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(ggq ggq_python_wrapper.cpp)
target_link_libraries(ggq PRIVATE ggq_library)
target_include_directories(ggq PRIVATE ${ggq_library_INCLUDE_DIRS})
set_target_properties(ggq PROPERTIES PREFIX "")

install(TARGETS ggq DESTINATION .)
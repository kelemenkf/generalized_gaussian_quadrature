set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
add_definitions(-DPYBIND11_DEBUG)


pybind11_add_module(ggq MODULE ggq_python_wrapper.cpp)

set_target_properties(ggq PROPERTIES PREFIX "")
target_link_libraries(ggq PRIVATE ggq_library)
target_include_directories(ggq PRIVATE ${ggq_library_INCLUDE_DIRS})

install(TARGETS ggq DESTINATION .)
#include "discretizer.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h> 

namespace py = pybind11;

PYBIND11_MODULE(ggq, m) {
    using InputFunctionType = std::function<double(const double&)>;

    py::class_<Discretizer>(m, "Discretizer")
        .def(py::init<int, double, double, double, Discretizer::InputFunctionType>())
        .def("evaluate", &Discretizer::evaluateFunction)
        .def("get_nodes", &Discretizer::determineFinalNodes)
        .def("get_endpoints", &Discretizer::getFinalEndpoints)
        .def("determine_final_nodes", &Discretizer::determineFinalNodes)
        .def("get_lower_bound", &Discretizer::getLowerBound)
        .def("get_upper_bound", &Discretizer::getUpperBound)
    ;
}
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
    ;
}
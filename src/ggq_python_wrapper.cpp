#include "discretizer.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// namespace py = pybind11;

// template class Discretizer<double>;

// PYBIND11_MODULE(ggq, m) {
//     using InputMethodType = double (InputClass::*)(const double&);
//     using InputFunctionType = double(*)(const double&);
//     py::class_<Discretizer<InputClass>>(m, "Discretizer")
//         .def(py::init<int, double, double, double, double, double>())
//     ;
// }
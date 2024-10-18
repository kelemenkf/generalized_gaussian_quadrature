#include "discretizer.hpp"
#include "ggq.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h> 
#include <numbers>
using namespace std::numbers;

namespace py = pybind11;

using InputFunction = std::variant
<
    std::function<double(double)>,
    std::function<double(double, double)>,
    std::function<double(double, double, double)>,
    std::function<double(double, double, double, double)>
>;


template<typename... Parameters>
void declare_function_handler(py::module& m, const std::string& suffix) {
    using Class = FunctionHandler<Parameters...>;
    std::string class_name = "FunctionHandler" + suffix;

    py::class_<Class>(m, class_name.c_str())
        .def(py::init<py::function, Parameters...>())
        .def("get_params", &Class::getNumberOfParameters)
        .def("call_function", &Class::callFunction)
        ;
}


template<typename... Parameters>
void declare_discretizer(py::module& m, const std::string& suffix) {
    using FHClass = FunctionHandler<Parameters...>;
    using Discretizer = Discretizer<FHClass>;
    std::string class_name = "Discretizer" + suffix;

    py::class_<Discretizer>(m, class_name.c_str())
        .def(py::init<int, double, double, double, FHClass>())
        .def("get_nodes", &Discretizer::determineFinalNodes)
        .def("get_endpoints", &Discretizer::getFinalEndpoints)
        .def("determine_final_nodes", &Discretizer::determineFinalNodes)
        .def("get_lower_bound", &Discretizer::getLowerBound)
        .def("get_upper_bound", &Discretizer::getUpperBound)
    ;
}


template<typename... Parameters>
void declare_quadrature(py::module& m, const std::string& suffix) {
    using FHClass = FunctionHandler<Parameters...>;
    using Quadrature = QuadratureRule<FHClass>;
    std::string class_name = "Quadrature" + suffix;

    py::class_<Quadrature>(m, class_name.c_str())
        .def(py::init<double, double, FHClass>())
        .def("discretize", &Quadrature::discretizeFunctions)
    ;
}


PYBIND11_MODULE(ggq, m) {
    // Bind different instantiations
    declare_function_handler<>(m, "0Param");
    declare_function_handler<std::vector<double>>(m, "1Param");
    declare_function_handler<std::vector<double>, std::vector<double>>(m, "2Param");
    declare_function_handler<std::vector<double>, std::vector<double>, std::vector<double>>(m, "3Param");

    declare_discretizer<>(m, "0Param");
    declare_discretizer<std::vector<double>>(m, "1Param");
    declare_discretizer<std::vector<double>, std::vector<double>>(m, "2Param");
    declare_discretizer<std::vector<double>, std::vector<double>, std::vector<double>>(m, "3Param");

    declare_quadrature<>(m, "0Param");
    declare_quadrature<std::vector<double>>(m, "1Param");
    declare_quadrature<std::vector<double>, std::vector<double>>(m, "2Param");
    declare_quadrature<std::vector<double>, std::vector<double>, std::vector<double>>(m, "3Param");
}
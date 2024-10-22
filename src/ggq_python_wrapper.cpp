#include "ggq.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h> 
#include <numbers>
using namespace std::numbers;

using InputFunction = std::variant
<
    py::function,
    std::function<double(const double&)>,
    std::function<double(const double&, const double&)>,
    std::function<double(const double&, const double&, const double&)>,
    std::function<double(const double&, const double&, const double&, const double&)>
>;


template<typename... Parameters>
void declare_function_handler(py::module& m, const std::string& suffix) {
    using Class = FunctionHandler<Parameters...>;
    std::string class_name = "FunctionHandler" + suffix;

    py::class_<Class>(m, class_name.c_str())
        .def(py::init<InputFunction, Parameters...>())
        .def("param_space", &Class::getParameterCombinations)
        ;
}


template<typename... Parameters>
void declare_quadrature(py::module& m, const std::string& suffix) {
    using FHClass = FunctionHandler<Parameters...>;
    using Quadrature = QuadratureRule<FHClass>;
    std::string class_name = "Quadrature" + suffix;

    py::class_<Quadrature>(m, class_name.c_str())
        .def(py::init<double, double, FHClass>())
        .def("calculate_quadrature", &Quadrature::calculateQuadratureNodes)
        .def("get_endpoints", &Quadrature::getConsolidatedEndpoints)
        .def("get_nodes", &Quadrature::getNodes)
    ;
}


PYBIND11_MODULE(ggq, m) {
    declare_function_handler<>(m, "0Param");
    declare_function_handler<std::vector<double>>(m, "1Param");
    declare_function_handler<std::vector<double>, std::vector<double>>(m, "2Param");
    declare_function_handler<std::vector<double>, std::vector<double>, std::vector<double>>(m, "3Param");

    declare_quadrature<>(m, "0Param");
    declare_quadrature<std::vector<double>>(m, "1Param");
    declare_quadrature<std::vector<double>, std::vector<double>>(m, "2Param");
    declare_quadrature<std::vector<double>, std::vector<double>, std::vector<double>>(m, "3Param");

    m.def("test_function", [](){ return "Carey nem mariah"; });
}
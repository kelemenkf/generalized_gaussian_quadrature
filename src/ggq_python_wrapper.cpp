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
        .def(py::init<double, double, FHClass, double, double>())
        .def(py::init<double, double, FHClass, double, double, size_t>())
        .def("calculate_quadrature", &Quadrature::calculateQuadratureNodes)
        .def("compress_functions", &Quadrature::compressFunctionSpace)
        .def("interpolate_basis", &Quadrature::obtainBasisCoefficients)
        .def("get_endpoints", &Quadrature::getConsolidatedEndpoints)
        .def("get_nodes", &Quadrature::getNodes)
        .def("get_values", &Quadrature::getValues)
        .def("get_u", &Quadrature::getCompressedBasis)
        .def("get_split_nodes", &Quadrature::getSplitNodes)
        .def("get_alphas", &Quadrature::getBasisCoefficients)
        .def("get_split_u", &Quadrature::getSplitCompressedBasis)
    ;
}


template<typename... Parameters>
void declare_interpolator(py::module& m, const std::string& suffix) {
    using FHClass = FunctionHandler<Parameters...>;
    using Quadrature = QuadratureRule<FHClass>;
    using Interpolator = Interpolator<FHClass>;
    std::string class_name = "Interpolator" + suffix;

    py::class_<Interpolator>(m, class_name.c_str())
        .def(py::init<int, double, double, FHClass, std::vector<double>>())
        .def("interpolate", &Interpolator::interpolateFunction)
        .def("get_alphas", &Interpolator::getAlphaVector)
        .def("eval", &Interpolator::evaluate)
    ;
}

 
void declare_evaluator(py::module& m) {
    py::class_<Evaluator>(m, "Evaluator")
        .def(py::init<std::vector<double>, const std::vector<double>&, double, double>())
        .def("get_evaluated_nodes", &Evaluator::getOutput)
        .def("evaluate_input", &Evaluator::evaluateInput)
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
    
    declare_interpolator<>(m, "0Param");
    declare_interpolator<std::vector<double>>(m, "1Param");
    declare_interpolator<std::vector<double>, std::vector<double>>(m, "2Param");
    declare_interpolator<std::vector<double>, std::vector<double>, std::vector<double>>(m, "3Param");

    declare_evaluator(m);
}
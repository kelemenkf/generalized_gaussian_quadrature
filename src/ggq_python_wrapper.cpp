#include "discretizer.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h> 

namespace py = pybind11;

using InputFunction = std::variant
<
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
        // Bind other methods here
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
}
#include "function_handler.hpp"


template class FunctionHandler<>;
template class FunctionHandler<std::vector<double>>;
template class FunctionHandler<std::vector<double>, std::vector<double>>;
template class FunctionHandler<std::vector<double>, std::vector<double>, std::vector<double>>;

template<typename... Parameters>
size_t FunctionHandler<Parameters...>::combinationIndex = 0;
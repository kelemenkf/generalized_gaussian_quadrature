#include <chrono>
#include "discretizer.hpp"
using namespace std::chrono;


double testFunction(const double& x)
{
    return pow(x, 2);
}


void timeSingleFunctionDetermineFinalNodes()
{
    auto start = high_resolution_clock::now();

    const int k = 30; 
    const double lowerBound = 1; 
    const double upperBound = 2;
    const double precision = 1e-6;

    Discretizer discretizer(k, precision, lowerBound, upperBound, testFunction);

    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> result;
    result = discretizer.determineFinalNodes();

    auto end = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(end - start);

    std::cout << "Discretizes x^2 in " << duration.count() << std::endl;
}


int main()
{
    timeSingleFunctionDetermineFinalNodes();

    return 0;
}
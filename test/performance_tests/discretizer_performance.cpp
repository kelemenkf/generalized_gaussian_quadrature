#include <chrono>
#include <numbers>
#include "discretizer.hpp"
using namespace std::chrono;
using namespace std::numbers;

double zeta(const double& alpha, const double& beta)
{
    return - beta * tan((alpha*pi) / 2);
}


double h(const double& t, const double& x, const double& alpha, const double& beta)
{
    return (x-zeta(alpha, beta)) * t + zeta(alpha, beta)*pow(t,alpha);
}


double IFTReCh(const double& t, const double& x, const double& alpha, const double& beta)
{
    return cos(h(t, x, alpha, beta)) * pow(e, -alpha);
}


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

    std::cout << "Discretizes second order polynomial in " << duration.count() << std::endl;
}


int main()
{
    timeSingleFunctionDetermineFinalNodes();

    return 0;
}
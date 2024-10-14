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


double piecewiseSmoothFunction(const double& x)
{
    if (x <= 1)
    {
        return x * x;
    }
    else 
    {
        return 2 - x;
    }
}


double highlyOscillatoryFunction(const double& x)
{
    return sin(50 * x);
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


void timePiecewiseSmoothFunctionDetermineFinalNodes()
{
    auto start = high_resolution_clock::now();

    const int k = 30; 
    const double lowerBound = 1; 
    const double upperBound = 2;
    const double precision = 1e-6;

    Discretizer discretizer(k, precision, lowerBound, upperBound, piecewiseSmoothFunction);

    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> result;
    result = discretizer.determineFinalNodes();

    auto end = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(end - start);

    std::cout << "Discretizes second order polynomial in " << duration.count() << std::endl;
}


void timeHighlyOsicllatoryFunctionDetermineFinalNodes()
{
    auto start = high_resolution_clock::now();

    const int k = 30; 
    const double lowerBound = 1; 
    const double upperBound = 2;
    const double precision = 1e-6;

    Discretizer discretizer(k, precision, lowerBound, upperBound, highlyOscillatoryFunction);

    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> result;
    result = discretizer.determineFinalNodes();

    auto end = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(end - start);

    std::cout << "Discretizes second order polynomial in " << duration.count() << std::endl;
}


int main()
{
    timeSingleFunctionDetermineFinalNodes();

    timePiecewiseSmoothFunctionDetermineFinalNodes();

    timeHighlyOsicllatoryFunctionDetermineFinalNodes();

    return 0;
}
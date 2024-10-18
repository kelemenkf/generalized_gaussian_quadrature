#include <chrono>
#include <numbers>
#include "discretizer.hpp"
#include "function_handler.hpp"
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
    return cos(h(t, x, alpha, beta)) * exp(-pow(t, alpha));
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



FunctionHandler<> handlerPolynomial(testFunction);
FunctionHandler<> handlerPiecewiseSmooth(piecewiseSmoothFunction);
FunctionHandler<> handlerHighlyOscillatory(highlyOscillatoryFunction);
std::vector<double> x = {10};
std::vector<double> alpha = {0.5};
std::vector<double> beta = {0};
FunctionHandler<std::vector<double>, std::vector<double>, std::vector<double>> handlerStable(IFTReCh, x, alpha, beta);



void timeSingleFunctionDetermineFinalNodes()
{
    auto start = high_resolution_clock::now();

    const int k = 30; 
    const double lowerBound = 1; 
    const double upperBound = 2;
    const double precision = 1e-6;

    Discretizer discretizer(k, precision, lowerBound, upperBound, handlerPolynomial);

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

    Discretizer discretizer(k, precision, lowerBound, upperBound, handlerPiecewiseSmooth);

    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> result;
    result = discretizer.determineFinalNodes();

    auto end = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(end - start);

    std::cout << "Discretizes piecewise smooth in " << duration.count() << std::endl;
}


void timeHighlyOsicllatoryFunctionDetermineFinalNodes()
{
    auto start = high_resolution_clock::now();

    const int k = 30; 
    const double lowerBound = 1; 
    const double upperBound = 2;
    const double precision = 1e-6;

    Discretizer discretizer(k, precision, lowerBound, upperBound, handlerHighlyOscillatory);

    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> result;
    result = discretizer.determineFinalNodes();

    auto end = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(end - start);

    std::cout << "Discretizes highly oscillatory in " << duration.count() << std::endl;
}

void timeStableDetermineFinalNodes()
{
    auto start = high_resolution_clock::now();

    const int k = 30; 
    const double lowerBound = 0; 
    const double upperBound = 20;
    const double precision = 1e-5;

    Discretizer discretizer(k, precision, lowerBound, upperBound, handlerStable);

    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> result;
    result = discretizer.determineFinalNodes();

    auto end = high_resolution_clock::now();

    auto duration = duration_cast<seconds>(end - start);

    std::cout << "Discretizes stable in " << duration.count() << std::endl;
}


int main()
{
    timeSingleFunctionDetermineFinalNodes();

    timePiecewiseSmoothFunctionDetermineFinalNodes();

    // timeHighlyOsicllatoryFunctionDetermineFinalNodes();

    timeStableDetermineFinalNodes();

    return 0;
}
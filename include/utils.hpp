#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <algorithm>
#include <iostream>
#include <map>
#include <cmath>
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <Eigen/Core>

using namespace boost::math::tools;
using namespace boost::numeric::ublas;
using namespace Eigen; 


template<typename T>
void displayVector(const std::vector<T>& vector)
{
    for (auto element = vector.cbegin(); element != vector.cend(); 
    ++element)
    {
        std::cout << *element << " " << std::endl;
    }
    std::cout << std::endl;
}



template<typename T>
T innerProduct(const std::vector<T>& input1, const std::vector<T>& input2)
{   
    T result = 0.0;
    if (input1.size() == input2.size())
    {
        for (size_t i = 0; i < input1.size(); ++i)
        {
            result += (input1[i] * input2[i]);
        }
    }

    return result;
}


template<typename T>
std::vector<T> convertBoostVectorToStd(const vector<T>& input)
{
    std::vector<T> result(input.size()); 

    for (size_t i = 0; i < input.size(); ++i)
    {
        result[i] = input(i);
    }

    return result;
}



template<typename T>
vector<T> convertStdVectorToBoost(const std::vector<T>& input)
{
    vector<T> result(input.size());

    for (size_t i = 0; i < input.size(); ++i)
    {
        result(i) = input[i];
    }

    return result;
}


template<typename T>
VectorXd convertStdVectorToEigen(const std::vector<T>& input)
{
    VectorXd result(input.size());

    for (size_t i = 0; i < input.size(); ++i)
    {
        result[i] = input[i];
    }

    return result;
}


template<typename MatrixType>
MatrixType removeColumnFromEigenMatrix(MatrixType& input, int columnToRemove, int n)
{
    MatrixType result(input.rows(), 0);

    for (int col = 0; col < input.cols(); ++col)
    {
        if (col != columnToRemove && col != (columnToRemove + n))
        {
            result.conservativeResize(Eigen::NoChange, result.cols() + 1);
            result.col(result.cols() - 1) = input.col(col);
        }
    }

    return result;
}

#endif
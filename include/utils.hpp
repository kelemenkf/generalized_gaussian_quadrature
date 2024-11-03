#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <algorithm>
#include <iostream>
#include <map>
#include <cmath>

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


#endif
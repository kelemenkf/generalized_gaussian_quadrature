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



#endif
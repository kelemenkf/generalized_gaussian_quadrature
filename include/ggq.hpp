#include <vector>

template <typename InputClass = double>
class QuadratureRule
{
private:
    double lowerBound;
    double upperBound;

    using InputMethodType = double (InputClass::*)(const std::vector<double>&);
    InputMethodType methodPtr;
    InputClass& objectRef;

    using InputFunctionType = double ()(const std::vector<double>);
    InputFunctionType functionPtr;


public:
    QuadratureRule(double lowerBoundInput, double upperBoundInput, InputFunctionType function = NULL, 
    InputMethodType mehtod = NULL, InputClass& inputObject = NULL);

    ~QuadratureRule();


protected: 
    InputFunctionType* getFunctionPointer()
    {
        return functionPtr;
    }


    InputMethodType* getMethodPointer()
    {
        return methodPtr;
    }


private: 
    static double validateLowerBound(double input);

    static double validateUpperBound(double input);
};
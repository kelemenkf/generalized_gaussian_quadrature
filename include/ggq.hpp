#include <vector>

template <typename InputClass>
class QuadratureRule
{
private:
    double lowerBound;
    double upperBound;

    using InputFunctionType = double (InputClass::*)(const std::vector<double>&);
    InputFunctionType functionPtr;
    InputClass& objectRef;


public:
    QuadratureRule(InputFunctionType function, InputClass& inputClass, double lowerBoundInput, double upperBoundInput);

    ~QuadratureRule();


private: 
    static double validateLowerBound(double input);

    static double validateUpperBound(double input);
};
#ifndef DISCRETIZER_HPP
#define DISCRETIZER_HPP

#include "ggq.hpp"
#include "interval_divider.hpp"

template<typename InputClass>
class Discretizer: public QuadratureRule<InputClass>
{
private: 
    int k;
    double precision;
    std::vector<double> endpoints;
    std::vector<double> measureVector;

    using InputMethodType = double (InputClass::*)(const double&);
    using InputFunctionType = double(*)(const double&);


public:
    Discretizer(int kInput, double precisionInput, double lowerBoundInput, double upperBoundInput, InputFunctionType inputFunctionPtr = nullptr, 
    InputMethodType inputMethodPtr = nullptr, InputClass* objectPtr = nullptr) 
    : QuadratureRule<InputClass>(lowerBoundInput, upperBoundInput, inputFunctionPtr, inputMethodPtr, objectPtr), 
    k(validateK(kInput)), precision(validatePrecision(precisionInput)) {
        endpoints = {this->lowerBound, this->upperBound};
    };

    ~Discretizer() {};

    void discretizationRoutine()
    {
        while(!evaluateStoppingCondition())
        {
            //Loop over endpoints, deteremine measure with interval_divider, update measurevector
            measureVector.reserve(endpoints.size() - 1);
            for (size_t i = 0; i < endpoints.size() - 1; ++i)
            {
                IntervalDivider<InputClass> divider;
                if (this->functionPtr)
                {
                    divider = IntervalDivider<InputClass>(k, endpoints[i], endpoints[i+1], this->functionPtr, nullptr, nullptr);
                }
                else if (this->methodPtr)
                {
                    divider = IntervalDivider<InputClass>(k, endpoints[i], endpoints[i+1], nullptr, this->methodPtr, this->objectPtr);   
                }
                divider.processInterval();

                double measure;
                measure = divider.getMeasure();

                //Smart update so only previously unused intervals use measures
                measureVector.push_back(measure);
            }
        }
    }


    std::vector<double> getFinalEndpoints()
    {
        return endpoints;
    }



private: 
    bool evaluateStoppingCondition()
    {
        bool stop = true;
        std::vector<std::vector<double>::iterator> impreciseSubintervalIndeces;
        for (size_t i = 0; i < measureVector.size(); ++i)
        {
            if (measureVector[i] >= precision)
            {
                stop = false;
                auto it = std::find(endpoints.begin(), endpoints.end(), endpoints[i+1]);
                if (it != endpoints.end())
                    impreciseSubintervalIndeces.push_back(it);
            }
        }
        if (!stop) 
            determineNewEndpoints(impreciseSubintervalIndeces);

        return stop;
    }


    void determineNewEndpoints(const std::vector<std::vector<double>::iterator>& impreciseSubintervalIndeces)
    {
        for (size_t i = impreciseSubintervalIndeces.size() - 1; i >= 0; --i)
        {
            double newPoint = calculateNewEndpoint(*std::prev(impreciseSubintervalIndeces[i]), *impreciseSubintervalIndeces[i]);
            endpoints.insert(impreciseSubintervalIndeces[i], newPoint);
        } 

        measureVector.clear();
    }   


    double calculateNewEndpoint(double lower, double upper)
    {
        return (lower + upper) / 2;
    }


    static int validateK(int inputK)
    {
        if (inputK > 0)
        {
            return inputK;
        }
        else
        {
            throw std::invalid_argument("k has to be positive");
        }
    }

    static double validatePrecision(double inputPrecision)
    {
        if (inputPrecision > 0)
        {
            return inputPrecision;
        }
        else 
        {
            throw std::invalid_argument("Precision has to be positive");
        }
    }


protected: 
    double getPrecision() const
    {
        return precision;
    }
};

#endif
#ifndef DISCRETIZER_HPP
#define DISCRETIZER_HPP

#include "interpolator.hpp"
#include "utils.hpp"
#include <tuple>

template<typename T>
class Discretizer
{
private: 
    int k;
    double precision;
    double lowerBound;
    double upperBound;
    std::vector<double> endpoints;
    std::vector<double> measureVector;
    T handler;


public:
    Discretizer(int kInput, double precisionInput, double lowerBoundInput, double upperBoundInput, const T& handlerInput) 
    : k(validateK(kInput)), precision(validatePrecision(precisionInput)), lowerBound(lowerBoundInput), upperBound(upperBoundInput), handler(handlerInput)
    {
        endpoints = {this->lowerBound, this->upperBound};
        determineFinalEndpoints();
    };

    ~Discretizer() {};


    double getLowerBound() const
    {
        return lowerBound;
    }

    double getUpperBound() const 
    {
        return upperBound;
    }


    void determineFinalEndpoints()
    {
        calculateMeasures();

        while(!evaluateStoppingCondition())
        {
            calculateMeasures();
        }
    }


    void calculateMeasures()
    {
        measureVector.reserve(endpoints.size() - 1);
        for (size_t i = 0; i < endpoints.size() - 1; ++i)
        {
            Interpolator divider(k, endpoints[i], endpoints[i+1], handler);

            divider.processInterval();

            double measure;
            measure = divider.getMeasure();

            //Smart update so only previously unused intervals use measures
            measureVector.push_back(measure);
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
        for (size_t i = impreciseSubintervalIndeces.size(); i-- > 0;)
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


    std::vector<double> getMeasureVector() const
    {
        return measureVector;
    }
};

#endif
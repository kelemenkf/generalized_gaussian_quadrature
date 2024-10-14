#ifndef DISCRETIZER_HPP
#define DISCRETIZER_HPP

#include "interval_divider.hpp"
#include <tuple>

class Discretizer: private QuadratureRule
{
private: 
    int k;
    double precision;
    std::vector<double> endpoints;
    std::vector<double> measureVector;


public:
    using InputFunctionType = std::function<double(const double&)>;

    Discretizer(int kInput, double precisionInput, double lowerBoundInput, double upperBoundInput, InputFunctionType inputFunctionPtr) 
    : QuadratureRule(lowerBoundInput, upperBoundInput, inputFunctionPtr), 
    k(validateK(kInput)), precision(validatePrecision(precisionInput))
    {
        endpoints = {this->lowerBound, this->upperBound};
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
            IntervalDivider divider(k, endpoints[i], endpoints[i+1], this->function);

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


    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> determineFinalNodes()
    {
        std::vector<std::vector<double>> nodes;
        std::vector<std::vector<double>> values;
        for (size_t i = 0; i < endpoints.size() - 1; ++i)
        {
            IntervalDivider divider(k / 2, endpoints[i], endpoints[i+1], this->function);
            divider.calculateMesh();
            divider.transformMesh();
            divider.evaluateFunctionOnTransformedMesh();
            nodes.push_back(divider.getTransformedMesh());
            values.push_back(divider.getLagrangeVector());
        }

        return {nodes, values};
    }


    double evaluateFunction(const double& value)
    {
        return function(value);
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


    std::vector<double> getMeasureVector() const
    {
        return measureVector;
    }
};

#endif
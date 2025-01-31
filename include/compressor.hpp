#ifndef COMPRESSOR_HPP
#define COMPRESSOR_HPP

#include <Eigen/Dense>
#include <numeric>
#include "utils.hpp"
using namespace Eigen;


template<typename T>
class Compressor
{
private:
    const T* quadrature;
    double quadraturePrecision;
    MatrixXd A;
    MatrixXd U;
    MatrixXd scaledU;
    MatrixXd R;
    std::vector<double> nodes;
    std::vector<double> weights;
    std::vector<std::vector<double>> values;
    std::vector<double> normalizingFactors;
    std::vector<std::vector<double>> scaledDiscardedU;
    VectorXd rVector;
    MatrixXd B;
    MatrixXd R_11;
    MatrixXd Q;
    VectorXd z;
    std::vector<size_t> selectedK;
    std::vector<double> chebyshevNodes; 
    std::vector<double> chebyshevWeights; 


public:
    Compressor(const T* quadratureInput, double quadraturePrecisionInput = 1e-3) : quadrature(quadratureInput), quadraturePrecision(quadraturePrecisionInput)
    {
        nodes = quadrature->getNodes();
        weights = quadrature->getWeights();
        values = quadrature->getValues();
        std::cout << "Compressor got " << values.size() << " functions as input." << std::endl;
        std::cout << "Compressor got a quadrature input with " << nodes.size() << " nodes." << std::endl;
        constructA();
        decomposeIntoQR();
        scaleU();
        calculateNormalizingFactors();
        discardFunctions(); 
        calculateBasisFunctionIntegrals();
        constructB();
        doubleOrthogonalization();
        solveSystem();
        formNewQuadrature();
        std::cout << "Compressor outputs a function space with " << scaledU.cols() << " functions." << std::endl;
        std::cout << "Compressor outputs a quadrature with " << chebyshevNodes.size() << " nodes." << std::endl; 
    }


    ~Compressor() 
    {

    }


    MatrixXd getA() const
    {
        return A;
    }


    MatrixXd getU() const
    {
        return U;
    }


    MatrixXd getScaledU() const
    {
        return scaledU;
    }


    std::vector<double> getNormalizingFactors() const 
    {
        return normalizingFactors;
    }


    std::vector<std::vector<double>> getCompressedBasis() const
    {
        return scaledDiscardedU;
    }


    VectorXd getRVector() const
    {
        return rVector;
    }


    MatrixXd getB() const
    {
        return B;
    }


    MatrixXd getQ() const
    {
        return Q;
    }


    MatrixXd getR_11() const
    {
        return R_11;
    }


    std::vector<size_t> getSelectedK() const
    {
        return selectedK;
    }


    std::vector<double> getChebyshevNodes() const
    {
        return chebyshevNodes;
    }


    std::vector<double> getChebyshevWeights() const
    {
        return chebyshevWeights;
    }


private:
    void constructA()
    {
        A.resize(nodes.size(), values.size());
        for (size_t column = 0; column < values.size(); ++column)
        {
            for (size_t row = 0; row < nodes.size(); ++row)
            {
                A(row, column) = values[column][row] * sqrt(weights[row]);
            }
        }
    }


    void decomposeIntoQR()
    {
        HouseholderQR<MatrixXd> qr(A);

        U = qr.householderQ();

        std::cout << "Initial orthogonalization outputs matrix with size " << U.rows() << " " << U.cols() << std::endl;
        R = qr.matrixQR().triangularView<Eigen::Upper>();
    }


    void scaleU()
    {
        scaledU = U;
        for (size_t row = 0; row < weights.size(); ++row)
        {
            scaledU.row(row) /= sqrt(weights[row]);
        }

        std::cout << "Scaled U shape " << scaledU.rows() << " " << scaledU.cols() << std::endl;
    }


    void calculateNormalizingFactors()
    {  
        for (size_t i = 0; i < R.rows(); ++i)
        {
            for (size_t j = 0; j < R.cols(); ++j)
            {
                if (i == j)
                {
                    if (R(i,j) < 0)
                    {
                        //Correction for Householder QR producing negative factors
                        normalizingFactors.push_back(-R(i,j));
                    }
                    else 
                    {
                        normalizingFactors.push_back(R(i,j));
                    }
                }
            }
        }
    }


    void discardFunctions()
    {
        int numberOfDiscardedFunctions = 0;

        std::cout << "Normalizing factors" << std::endl;

        displayVector(normalizingFactors);

        std::cout << "Precision " << quadraturePrecision << std::endl;

        for (size_t i = 0; i < normalizingFactors.size(); ++i)
        {
            if (normalizingFactors[i] > quadraturePrecision)
            {
                ++numberOfDiscardedFunctions; 
                Eigen::VectorXd column = scaledU.col(i);
                std::vector<double> vec(column.data(), column.data() + column.size());
                scaledDiscardedU.push_back(vec);
            }
        }
        std::cout << "Number of functions kept " << numberOfDiscardedFunctions << std::endl;
    }


    void calculateBasisFunctionIntegrals()
    {
        rVector.resize(scaledDiscardedU.size());
        for (size_t i = 0; i < scaledDiscardedU.size(); ++i)
        {
            double res = calculateBasisIntegral(i);
            rVector[i] = res;
        }
    }


    double calculateBasisIntegral(size_t index)
    {
        double result = 0; 

        if (scaledDiscardedU[index].size() != weights.size()) {
            throw std::runtime_error("Size mismatch between vectors.");
        }

        for (size_t i = 0; i < scaledDiscardedU[index].size(); ++i)
        {
            result += (scaledDiscardedU[index][i] * weights[i]);
        }

        return result;
    }


    void constructB()
    {
        B.resize(scaledDiscardedU.size(), weights.size());
        for (size_t row = 0; row < scaledDiscardedU.size(); ++row)
        {
            for (size_t column = 0; column < weights.size(); ++column)
            {
                B(row, column) = scaledDiscardedU[row][column] * sqrt(weights[column]);
            }
        }
    }


    void doubleOrthogonalization()
    {
        std::tuple<MatrixXd, MatrixXd, std::vector<size_t>> result = doublePivotedGramSchmidt(B);
        Q = std::get<0>(result);
        MatrixXd tempR = std::get<1>(result);
        std::vector<size_t> perm = std::get<2>(result);

        size_t k = Q.cols(); 
        selectedK.resize(k); 

        for (size_t i = 0; i < k; ++i)
        {
            selectedK[i] = perm[i];
        }

        R_11.resize(Q.rows(), Q.cols());
        for (size_t i = 0; i < Q.rows(); ++i)
        {
            for (size_t j = 0; j < Q.cols(); ++j)
            {
                R_11(i, j) = tempR(i, j);
            }
        }
    }


    void solveSystem()
    {
        VectorXd b = Q.transpose() * rVector;
        z = VectorXd::Zero(b.size());

        for (int i = R_11.rows() - 1; i >= 0; --i)
        {
            double sum = 0.0;
            for (int j = i+1; j < R_11.rows(); ++j)
            {
                sum += R_11(i,j) * z(j);
            }
            z[i] = (b[i] - sum) / R_11(i,i);
        }
    }


    void formNewQuadrature()
    {
        chebyshevWeights.resize(selectedK.size());
        chebyshevNodes.resize(selectedK.size());

        std::cout << "SelectedK size " << selectedK.size() << std::endl;

        for (size_t i = 0; i < selectedK.size(); ++i)
        {   
            chebyshevWeights[i] = z(i) * sqrt(weights[selectedK[i]]);
            chebyshevNodes[i] = nodes[selectedK[i]];
        }
    }


protected:
    std::tuple<MatrixXd, MatrixXd, std::vector<size_t>> doublePivotedGramSchmidt(MatrixXd& inputMatrix) 
    {
        int n = inputMatrix.cols();
        int m = inputMatrix.rows();
        Q = MatrixXd::Zero(m, m);
        R = MatrixXd::Zero(m, n);
        MatrixXd V = inputMatrix;

        std::vector<size_t> perm(n);
        for (size_t i = 0; i < n; ++i)
        {
            perm[i] = i;
        }

        std::vector<double> norms(n);
        for (size_t j = 0; j < n; ++j)
        {
            norms[j] = V.col(j).norm();
        }
    

        for (size_t k = 0; k < m; ++k)
        {
            size_t maxNormId = k;
            double maxNorm = norms[k];
            for(size_t j = k+1; j < n; ++j)
            {
                if (norms[j] > maxNorm)
                {
                    maxNormId = j;
                    maxNorm = norms[j];
                }
            }


            if (maxNormId != k)
            {
                V.col(k).swap(V.col(maxNormId));
                std::swap(perm[k], perm[maxNormId]);

                for (size_t i = 0; i < k; ++i) 
                {
                    std::swap(R(i, k), R(i, maxNormId));
                }
            }

            VectorXd vk = V.col(k);

            double rkk = vk.norm();

            if (rkk < 1e-10) { 
                continue;
            }

            R(k,k) = rkk;
            VectorXd qk = vk / rkk;
            Q.col(k) = qk;       

            for (size_t j = k + 1; j < n; ++j)
            {
                double rkj = V.col(j).dot(Q.col(k));

                R(k, j) = rkj;
                V.col(j) -= Q.col(k) * rkj;
            }


            for (size_t j = k + 1; j < n; ++j)
            {
                double rkj = V.col(j).dot(Q.col(k));
                R(k,j) += rkj;
                V.col(j) -= Q.col(k) * rkj;
            }

            for (size_t j = k + 1; j < n; ++j)
            {
                norms[j] = V.col(j).norm();
            }
        } 


        return {Q, R, perm};
    }
};
 

#endif
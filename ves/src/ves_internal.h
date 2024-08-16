#pragma once
#ifndef VES_INTERNAL
#define VES_INTERNAL
#include <vector>
#include <Eigen/Eigen/Core>

namespace ves {
    
    double _pow(double base, int exp);

    // Returns vector with positions (0, 1) of edge node for VEM, follow the Gauss-Lobatto convention.
    const std::vector<double> EdgeNodePositions(const int order); // order is the order of the virtual element, not that of the quadrature
    // Returns vector with positions (0, 1) of points for Lobatto quadrature associated with sides of VEM element of input order.
    const std::vector<double> LobattoNodePositions(const int order); // order is the order of the virtual element, not that of the quadrature

    // Returns vector with pairs of positions (0,1) and weights (0,1) for the full Lobatto quadrature, ordered by position
    const std::vector<std::array<double, 2>> LobattoQuadratureOrdered(const int order); // order is the order of the virtual element, not that of the quadrature

    // Evaluates Lagrange polynomial based on the points being interpolated (NodePositions), the index of the corresponding function, and the value at which to evaluate (xsi)
    const double LagrangePolynomialEvaluation(const std::vector<double>& NodePositions, const size_t functionIndex, const double xsi);

    // Returns the external normal vector for a line given by start->end
    const Eigen::Vector2d WeightedNormalFromLine(const Eigen::Vector2d& start, const Eigen::Vector2d& end);
    const Eigen::Vector2d NormalFromLine(const Eigen::Vector2d& start, const Eigen::Vector2d& end);


}

#endif
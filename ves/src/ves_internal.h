#pragma once
#ifndef VES_INTERNAL
#define VES_INTERNAL
#include <vector>
#include <Eigen/Eigen/Core>

namespace ves {
    // Returns vector with positions (0,1) of edge node for VEM, follow the Gauss-Lobatto convention.
    const std::vector<double> EdgeNodePositions(const int order); // order is the order of the virtual element, not that of the quadrature

    // Returns vector with pairs of positions (0,1) and weights (0,1) for the full Lobatto quadrature, ordered by position
    const std::vector<std::array<double, 2>> LobattoQuadratureOrdered(const int order); // order is the order of the virtual element, not that of the quadrature

    // Returns the external normal vector for a line given by start->end
    const Eigen::Vector2d WeightedNormalFromLine(const Eigen::Vector2d& start, const Eigen::Vector2d& end);
    const Eigen::Vector2d NormalFromLine(const Eigen::Vector2d& start, const Eigen::Vector2d& end);
}

#endif
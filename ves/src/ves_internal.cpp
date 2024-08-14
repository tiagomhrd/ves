#include "ves_internal.h"
#include "mnl/include/glq.hpp"
#include <algorithm>
namespace ves {
    const std::vector<double> ves::EdgeNodePositions(const int order)
    {
        std::vector<double> out(size_t(order - 1));
        const auto gl = mnl::GaussLobattoR(2 * (order + 1) - 3);
        std::transform(gl.cbegin() + 2, gl.cend(), out.begin(), [](const auto& pair){ return pair[0]; });
        std::sort(out.begin(), out.end());
        return out;
    }
    const std::vector<std::array<double, 2>> LobattoQuadratureOrdered(const int order)
    {
        auto gl = mnl::GaussLobattoR(2 * (order + 1) - 3);
        std::sort(gl.begin(), gl.end(), [](const auto& p1, const auto& p2){ return p1[0] < p2[0];});
        return gl;
    }
    const Eigen::Vector2d WeightedNormalFromLine(const Eigen::Vector2d &start, const Eigen::Vector2d &end)
    {
        return Eigen::Vector2d(end(1) - start(1), start(0) - end(0));
    }
    const Eigen::Vector2d NormalFromLine(const Eigen::Vector2d &start, const Eigen::Vector2d &end)
    {
        return WeightedNormalFromLine(start, end).normalized();
    }
}
#include "V2D.h"

#include <Eigen/Eigen/Dense>
#include <ptp/ptp/src/ptp.h>
#include <mnl/include/mnl.hpp>

#include "ves_internal.h"

namespace ves {
    V2D::V2D(const std::vector<Eigen::Vector2d> &polygon, const int order, const int maxMonomialOrder)
        : m_Polygon{polygon}, m_Order{order}, m_InvDiameter{pow(ptp::Polygon2D::Diameter(polygon), -1.)}
    {
        // Centroid
        std::vector<double> integrals = ptp::Polygon2D::MonomialIntegrals(m_Polygon, 1);
        m_Centroid = Eigen::Vector2d(integrals[1], integrals[2]) / integrals[0];
        
        // Scaled Monomial Integrals
        m_SMIntegrals = ScaledMonomialIntegrals(std::max(maxMonomialOrder, 2 * order));

        // Method Matrices
        Init();
    }
    void V2D::Init()
    {
        const Eigen::MatrixXd GGrad = GGrad_Impl();
        const Eigen::MatrixXd BGrad = BGrad_Impl();
        m_PiGrad = GGrad.ldlt().solve(BGrad);

        const Eigen::MatrixXd G0 = G0_Impl();
        const Eigen::MatrixXd B0 = B0_Impl();
        m_Pi0 = G0.ldlt().solve(B0);
    }
    const double V2D::Area() const
    {
        return m_SMIntegrals[0];
    }
    const Eigen::Vector2d V2D::Centroid() const
    {
        return m_Centroid;
    }
    const Eigen::MatrixXd V2D::D() const
    {
        return D_Impl();
    }
    const Eigen::MatrixXd V2D::GGradTilde() const
    {
        return GGradTilde_Impl();
    }
    const Eigen::MatrixXd V2D::PiGrad() const
    {
        return m_PiGrad;
    }
    const Eigen::MatrixXd V2D::Pi0() const
    {
        return m_Pi0;
    }
    const Eigen::MatrixXd V2D::GGrad() const
    {
        return GGrad_Impl();
    }
    const Eigen::MatrixXd V2D::G0() const
    {
        return G0_Impl();
    }
    const Eigen::MatrixXd V2D::BGrad() const
    {
        return BGrad_Impl();
    }
    const double V2D::SM(const int alpha, const Eigen::Vector2d &pos) const
    {
        if (alpha == -1)
            return 0.0;
        if (alpha == 0)
            return 1.0;
        const Eigen::Vector2d scaledPos = ScaledCoord(pos);
        return _pow(scaledPos(0), mnl::PSpace2D::Exponent(alpha, 0)) * _pow(scaledPos(1), mnl::PSpace2D::Exponent(alpha, 1));
    }
    const double V2D::SMIntegral(const int alpha) const
    {
        if (alpha < 0)
            return 0.0;
        return m_SMIntegrals[(size_t) alpha];
    }
    const double *V2D::IntegralData() const
    {
        return m_SMIntegrals.data();
    }
    const Eigen::Vector2d V2D::ScaledCoord(const Eigen::Vector2d &pos) const
    {
        return (pos - m_Centroid) * m_InvDiameter;
    }
    const std::vector<double> V2D::ScaledMonomialIntegrals(const int maxOrder) const
    {
        // Get scaled version of polygon
        std::vector<Eigen::Vector2d> scaledPoints(m_Polygon.size());
        std::transform(m_Polygon.cbegin(), m_Polygon.cend(), scaledPoints.begin(), [this](const auto& pt) { return this->ScaledCoord(pt); }); // 
        
        // Compute integrals
        std::vector<double> Integrals = ptp::Polygon2D::MonomialIntegrals(scaledPoints, maxOrder);
        
        // Adjust jacobian
        const double jacobian = _pow(ptp::Polygon2D::Diameter(m_Polygon), 2);
        std::transform(Integrals.begin(), Integrals.end(), Integrals.begin(), [&jacobian](double integral){ return integral * jacobian; });
        
        return Integrals;
    }
    const Eigen::MatrixXd V2D::D_Impl() const
    {
        const int nV = (int) m_Polygon.size();
        const int nInner = mnl::PSpace2D::SpaceDim(m_Order - 2);
        const int nDOF = m_Order * nV + nInner;
        const int nk = mnl::PSpace2D::SpaceDim(m_Order);
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(nDOF, nk);
        // Boundary terms
        // Vertices
        D.block(0, 0, m_Order * nV, 1).setOnes();
        int v = 0;
        for (const auto& point : m_Polygon){
            for (int alpha = 1; alpha < nk; ++alpha)
                D(v, alpha) = SM(alpha, point);
            ++v;
        }
        if (m_Order == 1)
            return D;

        // Edge Nodes
        const auto edgePositions = EdgeNodePositions(m_Order);
        const int nEdgePoints = (int) edgePositions.size();
        for (int v = 0; v < nV; ++v){
            for (int e = 0; e < nEdgePoints; ++e){
                const double xi = edgePositions[e]; 
                const Eigen::Vector2d point = (1. - xi) * m_Polygon[v] + xi * m_Polygon[(v + 1) % nV];
                for (int alpha = 1; alpha < nk; ++alpha)
                    D(nV + v * nEdgePoints + e, alpha) = SM(alpha, point);
            }
        }

        // Internal DOFs
        for (int alpha = 0; alpha < nInner; ++alpha)
            for (int beta = 0; beta < nk; ++beta)
                D(m_Order * nV + alpha, beta) = SMIntegral(mnl::PSpace2D::Product(alpha, beta));
        D.block(m_Order * nV, 0, nInner, nk) /= Area();
        return D;
    }
    const Eigen::MatrixXd V2D::GGradTilde_Impl() const
    {
        const int nk = mnl::PSpace2D::SpaceDim(m_Order);
        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(nk, nk);
        for (int alpha = 1; alpha < nk; ++alpha) {
            const int   alphaX = mnl::PSpace2D::D(alpha, 0),
                        alphaY = mnl::PSpace2D::D(alpha, 1),
                        alphaEX = mnl::PSpace2D::Exponent(alpha, 0), 
                        alphaEY = mnl::PSpace2D::Exponent(alpha, 1);
            for (int beta = 1; beta < nk; ++beta){
                const int   betaX = mnl::PSpace2D::D(beta, 0),
                            betaY = mnl::PSpace2D::D(beta, 1),
                            betaEX = mnl::PSpace2D::Exponent(beta, 0), 
                            betaEY = mnl::PSpace2D::Exponent(beta, 1);
                
                if (alphaX * betaX >= 0)
                    G(alpha, beta) += alphaEX * betaEX * SMIntegral(mnl::PSpace2D::Product(alphaX, betaX));
                if (alphaY * betaY >= 0)
                    G(alpha, beta) += alphaEY * betaEY * SMIntegral(mnl::PSpace2D::Product(alphaY, betaY));
            }
        }
        G *= m_InvDiameter * m_InvDiameter;
        return G;
    }
    const Eigen::MatrixXd V2D::GGrad_Impl() const
    {
        const int nk = mnl::PSpace2D::SpaceDim(m_Order);
        Eigen::MatrixXd G = GGradTilde_Impl();

        // Add P0
        if (m_Order == 1){
            // First order P0(m_alpha) is the average value of m_alpha on the vertices
            // This can be obtained directly from the D matrix, which contains every monomial evaluated at every node.
            const int nV = (int) m_Polygon.size();
            G.row(0) = D_Impl().transpose() * Eigen::VectorXd::Ones(nV) / nV;
            return G;
        }

        // P0(m_alpha) for higher order is the average value over the domain.
        // Integrals for alpha = 1 and 2 is 0, because scaled monomials have origin on the centroid
        for (int alpha = 3; alpha < nk; ++alpha)
            G(0, alpha) = SMIntegral(alpha);
        G.row(0) /= Area();
        // Integral for alpha = 0 is left out to avoid unnecessary additional division
        G(0,0) = 1.;
        return G;
    }
    const Eigen::MatrixXd V2D::BGrad_Impl() const
    {
        const int nV = (int) m_Polygon.size();
        const int nInner = mnl::PSpace2D::SpaceDim(m_Order - 2);
        const int nDOF = m_Order * nV + nInner;
        const int nk = mnl::PSpace2D::SpaceDim(m_Order);
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nk, nDOF);

        // Boundary term for vertices
        const auto& lobattoQuad = LobattoQuadratureOrdered(m_Order);
        const double vertexWeight = lobattoQuad[0][1];
        for (int v = 1; v <= nV; ++v) {
            const int previous = v - 1,
                      current = v % nV,
                      next = (v + 1) % nV;
            const Eigen::Vector2d weightedNormalCombination = 
                WeightedNormalFromLine(m_Polygon[previous], m_Polygon[current]) + 
                WeightedNormalFromLine(m_Polygon[current], m_Polygon[next]);

            for (int alpha = 1; alpha < nk; ++alpha){
                const Eigen::Vector2d gradAlpha(
                    mnl::PSpace2D::Exponent(alpha, 0) * SM(mnl::PSpace2D::D(alpha, 0), m_Polygon[current]) * m_InvDiameter,
                    mnl::PSpace2D::Exponent(alpha, 1) * SM(mnl::PSpace2D::D(alpha, 1), m_Polygon[current]) * m_InvDiameter);
                B(alpha, current) += gradAlpha.dot(weightedNormalCombination) * vertexWeight;
            }
        }

        if (m_Order == 1) {
            // If first order, add P0 and return
            B.row(0) = Eigen::VectorXd::Ones(nV) / nV;
            return B;
        }

        // Boundary term for edges
        const int nEdgeNodes = (int)(lobattoQuad.size() - 2);
        for (int v = 1; v < nV + 1; ++v) {
            const int start = v - 1,
                      end = v % nV;
            const Eigen::Vector2d weightedNormal = WeightedNormalFromLine(m_Polygon[start], m_Polygon[end]);
            for (int e = 0; e < nEdgeNodes; ++e){
                const Eigen::Vector2d pos = (1 - lobattoQuad[1 + e][0]) * m_Polygon[start] + lobattoQuad[1 + e][0] * m_Polygon[end];
                for (int alpha = 1; alpha < nk; ++alpha){
                    const Eigen::Vector2d gradAlpha(
                        mnl::PSpace2D::Exponent(alpha, 0) * SM(mnl::PSpace2D::D(alpha, 0), pos) * m_InvDiameter,
                        mnl::PSpace2D::Exponent(alpha, 1) * SM(mnl::PSpace2D::D(alpha, 1), pos) * m_InvDiameter);
                    B(alpha, nV + nEdgeNodes * start + e) += gradAlpha.dot(weightedNormal) * lobattoQuad[1 + e][1];
                }
            }
        }

        // Laplacian term
        const double area = SMIntegral(0);
        for (int beta = 0; beta < nInner; ++beta){
            const int alphaX = mnl::PSpace2D::AD(mnl::PSpace2D::AD(beta, 0), 0),
                      alphaY = mnl::PSpace2D::AD(mnl::PSpace2D::AD(beta, 1), 1);
            const int expX = mnl::PSpace2D::Exponent(alphaX, 0),
                      expY = mnl::PSpace2D::Exponent(alphaY, 1);
            B(alphaX, m_Order * nV + beta) -= area * m_InvDiameter * m_InvDiameter * expX * (expX - 1);
            B(alphaY, m_Order * nV + beta) -= area * m_InvDiameter * m_InvDiameter * expY * (expY - 1);
        }

        // P0
        B(0, m_Order * nV) = 1.0;

        return B;
    }
    const Eigen::MatrixXd V2D::G0_Impl() const
    {
        const int nk = mnl::PSpace2D::SpaceDim(m_Order);
        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(nk, nk);
        for (int alpha = 0; alpha < nk; ++alpha) {
            G(alpha, alpha) = SMIntegral(mnl::PSpace2D::Product(alpha, alpha));
            for (int beta = alpha + 1; beta < nk; ++beta)
                G(alpha, beta) = G(beta, alpha) = SMIntegral(mnl::PSpace2D::Product(alpha, beta));
        }

        return G;
    }
    const Eigen::MatrixXd V2D::B0_Impl() const
    {
        const int nV = (int) m_Polygon.size();
        const int nInner = mnl::PSpace2D::SpaceDim(m_Order - 2);
        const int nDOF = m_Order * nV + nInner;
        const int nk = mnl::PSpace2D::SpaceDim(m_Order);
        Eigen::MatrixXd B = G0() * PiGrad();
        B.block(0, 0, nInner, nDOF).setZero();
        B.block(0, m_Order * nV, nInner, nInner) = Eigen::MatrixXd::Identity(nInner, nInner) * Area();
        return B;
    }
}
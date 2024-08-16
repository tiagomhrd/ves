#include "SV2D.h"

#include <mnl/include/mnl.hpp>
#include <mnl/include/pnl.hpp>
#include <mnl/include/glq.hpp>
#include <ptp/ptp/src/ptp.h>

#include "ves_internal.h"

namespace ves {
    SV2D::SV2D(const std::vector<Eigen::Vector2d> &polygon, const int order, const int maxMonomialOrder)
        : m_Polygon{polygon}, m_Order{order}, m_InvDiameter{pow(ptp::Polygon2D::Diameter(polygon), -1.0)}
    {
        Init();
    }
    void SV2D::Init()
    {
        // Centroid
        std::vector<double> integrals = ptp::Polygon2D::MonomialIntegrals(m_Polygon, 1);
        m_Centroid = Eigen::Vector2d(integrals[1], integrals[2]) / integrals[0];
        
        // Inner Order
        const int concavityOrder = CheckConcavity();
        m_InnerOrder = std::max(-1, m_Order - (int)ptp::Polygon2D::UniqueSides(m_Polygon).size() + concavityOrder);
        
        // Integrals
        m_SMIntegrals = ScaledMonomialIntegrals(std::max(2 * m_Order, (m_InnerOrder > -1 ? 0 : m_InnerOrder + concavityOrder) + m_Order));

        // Serendipity Projector
        const Eigen::MatrixXd D = D_Impl();
        const auto DT = D.transpose();
        m_PiS = (DT * D).ldlt().solve(DT);

        // Derivative Projectors
        const auto G0DSolver = G0_Impl().ldlt();
        const auto [B0Dx, B0Dy] = B0D_Impl();
        m_Pi0Dx = G0DSolver.solve(B0Dx);
        m_Pi0Dy = G0DSolver.solve(B0Dy);
    }
    const double SV2D::Area() const
    {
        return m_SMIntegrals[0];
    }
    const Eigen::Vector2d SV2D::Centroid() const
    {
        return m_Centroid;
    }
    const Eigen::MatrixXd SV2D::PiS() const
    {
        return m_PiS;
    }
    const Eigen::MatrixXd SV2D::Pi0() const
    {
        return m_PiS;
    }
    const Eigen::MatrixXd SV2D::Pi0Dx() const
    {
        return m_Pi0Dx;
    }
    const Eigen::MatrixXd SV2D::Pi0Dy() const
    {
        return m_Pi0Dy;
    }
    const double SV2D::SM(const int alpha, const Eigen::Vector2d &pos) const
    {
        if (alpha == -1)
            return 0.0;
        if (alpha == 0)
            return 1.0;
        const Eigen::Vector2d scaledPos = ScaledCoord(pos);
        return _pow(scaledPos(0), mnl::PSpace2D::Exponent(alpha, 0)) * _pow(scaledPos(1), mnl::PSpace2D::Exponent(alpha, 1));
    }
    const double SV2D::SMIntegral(const int alpha) const
    {
        if (alpha < 0)
            return 0.0;
        return m_SMIntegrals[(size_t) alpha];    
    }
    const double *SV2D::IntegralData() const
    {
        return m_SMIntegrals.data();
    }
    
    const std::vector<double> SV2D::ScaledMonomialIntegrals(const int maxOrder) const
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
    const Eigen::MatrixXd SV2D::D_Impl() const
    {
        const int nv = (int) m_Polygon.size();
        const bool convex = m_ConcavityFunction.empty();
        const int nInner = mnl::PSpace2D::SpaceDim(m_InnerOrder);
        const int nDof = m_Order * nv + (convex ? nInner : 2. * nInner);
        const int nk = mnl::PSpace2D::SpaceDim(m_Order);

        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(nDof, nk);

        const auto eNodePos = ves::EdgeNodePositions(m_Order);
        D.block(0, 0, 1, m_Order * nv) = Eigen::VectorXd::Ones(m_Order * nv); // Case for alpha = 0 treated separately 
        for (int v = 0; v < nv; ++v) {
            Eigen::Vector2d startScaledPos = ScaledCoord(m_Polygon[v]);
            Eigen::Vector2d endScaledPos = ScaledCoord(m_Polygon[(v + 1) % nv]);
            for (int alpha = 1; alpha < nk; ++alpha)
                D(v, alpha) = _pow(startScaledPos(0), mnl::PSpace2D::Exponent(alpha, 0)) * _pow(startScaledPos(1), mnl::PSpace2D::Exponent(alpha, 1));
            for (int i = 0; i < m_Order - 1; ++i) {
                Eigen::Vector2d edgeScaledPos = (1. - eNodePos[i]) * startScaledPos + eNodePos[i] * endScaledPos;
                for (int alpha = 1; alpha < nk; ++alpha)
                    D(nv + (m_Order - 1) * v + i, alpha) = _pow(edgeScaledPos(0), mnl::PSpace2D::Exponent(alpha, 0)) * _pow(edgeScaledPos(1), mnl::PSpace2D::Exponent(alpha, 1));
            }
        }

        if (m_InnerOrder == -1)
            return D;

        // Laplacian term if internal DOFs exist
        for (size_t beta{}; beta < nInner; ++beta) {
            for (size_t alpha{}; alpha < nk; ++alpha) {
                D(m_Order * nv + beta, alpha) += SMIntegral(mnl::PSpace2D::Product(alpha, beta));
            }
        }

        if (convex)
            return D;

        for (size_t beta{}; beta < nInner; ++beta) {
            for (size_t alpha{}; alpha < nk; ++alpha) {
                for (const auto& [index, val] : m_ConcavityFunction)
                    D(m_Order * nv + nInner + beta, alpha) += SMIntegral(mnl::PSpace2D::Product(mnl::PSpace2D::Product(alpha, beta), index)) * val;
            }
        }

        return D;
    }
    const Eigen::MatrixXd SV2D::G0_Impl() const
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
    const std::array<Eigen::MatrixXd, 2> SV2D::B0D_Impl() const
    {
        const int nkD = mnl::PSpace2D::SpaceDim(m_Order - 1);
        const int nv = (int) m_Polygon.size();
        const int nInner = mnl::PSpace2D::SpaceDim(m_InnerOrder);
        const int nDof = m_Order * nv + nInner;
        Eigen::MatrixXd B0Dx = Eigen::MatrixXd::Zero(nkD, nDof);
        Eigen::MatrixXd B0Dy = Eigen::MatrixXd::Zero(nkD, nDof);


        // Boundary term
        const auto lobattoPositions = LobattoNodePositions(m_Order);
        const size_t nEdgePoints = lobattoPositions.size() - 2;
        for (int n = 1; n <= m_Order; ++n) {
            const auto quadrature = mnl::GaussLegendreRN(n);
            const int startAlpha = mnl::PSpace2D::SpaceDim(2 * (n - 1) - 1),
                      endAlpha = mnl::PSpace2D::SpaceDim(2 * n - 1);
            for (size_t start{}; start < nv; ++start) {
                const size_t end = (start + 1) % nv;
                const Eigen::Vector2d weightedNormal = 
                    WeightedNormalFromLine(m_Polygon[start], m_Polygon[end]); // normal * length
                for (const auto&[xsi, weight] : quadrature) {
                    const Eigen::Vector2d pos = (1. - xsi) * m_Polygon[start] + xsi * m_Polygon[end];
                    const double startValue = LagrangePolynomialEvaluation(lobattoPositions, 0, xsi),
                                 endValue = LagrangePolynomialEvaluation(lobattoPositions, nEdgePoints + 1, xsi);
                    for (int alpha = startAlpha; alpha < endAlpha; ++alpha){
                        const double alphaValue = SM(alpha, pos);
                        B0Dx(alpha, start) += alphaValue * startValue * weightedNormal(0) * weight;
                        B0Dy(alpha, start) += alphaValue * startValue * weightedNormal(1) * weight;
                        B0Dx(alpha, end) += alphaValue * endValue * weightedNormal(0) * weight;
                        B0Dy(alpha, end) += alphaValue * endValue * weightedNormal(1) * weight;
                        for (size_t e = 0; e < nEdgePoints; ++e){
                            const size_t index = nv + (m_Order - 1) * start + e;
                            const double value = LagrangePolynomialEvaluation(lobattoPositions, e + 1, xsi);
                            B0Dx(alpha, index) += alphaValue * value * weightedNormal(0) * weight;
                            B0Dy(alpha, index) += alphaValue * value * weightedNormal(1) * weight;
                        }
                    }
                }
            }
        }

        // First order case has only vertex contributions to the boundary term
        if (m_Order == 1)
            return {B0Dx, B0Dy};

        // Interior term (Laplacian)

        // Contribution of internal DOFs (if existing)
        for (int beta = 0; beta < nInner; ++beta) {
            const int IdX = mnl::PSpace2D::AD(beta, 0),
                      IdY = mnl::PSpace2D::AD(beta, 1);

            const int expXIdX = mnl::PSpace2D::Exponent(IdX, 0),
                      expYIdY = mnl::PSpace2D::Exponent(IdY, 1);
            
            const int i = m_Order * nv + beta;
            B0Dx(IdX, i) -= m_InvDiameter * expXIdX;
            B0Dy(IdY, i) -= m_InvDiameter * expYIdY;
        }
        // Contribution of the Serendipity projection
        for (int beta = nInner, nLaplacianSpace = mnl::PSpace2D::SpaceDim(m_Order - 2); beta < nLaplacianSpace; ++beta) { // Works better this way
            const int IdX = mnl::PSpace2D::AD(beta, 0),
                      IdY = mnl::PSpace2D::AD(beta, 1);

            const int expXIdX = mnl::PSpace2D::Exponent(IdX, 0),
                      expYIdY = mnl::PSpace2D::Exponent(IdY, 1);

            for (int alpha{}, nk = mnl::PSpace2D::SpaceDim(m_Order); alpha < nk; ++alpha) {
                int pindex = mnl::PSpace2D::Product(beta, alpha);
                const Eigen::VectorXd aux = m_InvDiameter * SMIntegral(pindex) * m_PiS.row(alpha);
                B0Dx.row(IdX) -= aux * expXIdX;
                B0Dy.row(IdY) -= aux * expYIdY;
            }
        }

        return { B0Dx, B0Dy };
    }
    int SV2D::CheckConcavity()
    {
        auto RESides = ptp::Polygon2D::UniqueReentrantSides(m_Polygon);
        if (RESides.empty())
            return 0;

        mnl::pnl2D f;
        f.Terms[0] = 1.;
        f.Terms.reserve(mnl::PSpace2D::SpaceDim((int)RESides.size()));
        for (const auto& v : RESides) {
            mnl::pnl2D s;
            s.Terms.reserve(3);
            for (int i = 0; i < 3; ++i)
                s.Terms[i] = v[i];
            f *= s;
        }
        for (const auto& [key, val] : f.Terms)
            m_ConcavityFunction[key] = val;

        return (int)RESides.size();
    }
}
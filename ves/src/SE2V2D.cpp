#include "SE2V2D.h"

#include <mnl/include/mnl.hpp>
#include <mnl/include/pnl.hpp>
#include <mnl/include/glq.hpp>
#include <ptp/ptp/src/ptp.h>
#include <Eigen/Eigen/Dense>

#include "ves_internal.h"

namespace ves {
    SE2V2D::SE2V2D(const std::vector<Eigen::Vector2d> &polygon, const int order, const int gradOrder, const int maxMonomialOrder)
        : m_Polygon{polygon}, m_Order{order}, m_InvDiameter{pow(ptp::Polygon2D::Diameter(polygon), -1.0)}
    {
        // Centroid
        std::vector<double> integrals = ptp::Polygon2D::MonomialIntegrals(m_Polygon, 1);
        m_Centroid = Eigen::Vector2d(integrals[1], integrals[2]) / integrals[0];
        
        // Inner Order
        const int concavityOrder = CheckConcavity();
        m_InnerOrder = std::max(-1, m_GradOrder - 1 - (int)ptp::Polygon2D::UniqueSides(m_Polygon).size() + concavityOrder);

        // Grad Order
        m_GradOrder = (gradOrder > -1 ? gradOrder : StabFreeGradOrder());
        
        // Integrals
        const int maxOrder = std::max(maxMonomialOrder,                                             // Override
                             std::max(2 * m_GradOrder,                                              // G0D 
                             (m_InnerOrder > -1 ? 0 : m_InnerOrder + concavityOrder) + m_Order));
        m_SMIntegrals = ScaledMonomialIntegrals(maxOrder);
        
        Init();
    }
    void SE2V2D::Init()
    {
        const Eigen::MatrixXd D = D_Impl();
        {
            const Eigen::MatrixXd DT = D.transpose();
            m_PiS = (DT * D).ldlt().solve(DT);
        }
        const Eigen::MatrixXd G0D = G0D_Impl();
        const auto [B0Dx, B0Dy] = B0D_Impl();
        const auto G0DSolver = G0D.ldlt();
        m_Pi0Dx = G0DSolver.solve(B0Dx);
        m_Pi0Dy = G0DSolver.solve(B0Dy);
    }    
    const double SE2V2D::Area() const
    {
        return SMIntegral(0);
    }
    const Eigen::Vector2d SE2V2D::Centroid() const
    {
        return m_Centroid;
    }
    const Eigen::MatrixXd SE2V2D::D() const
    {
        return D_Impl();
    }
    const Eigen::MatrixXd SE2V2D::PiS() const
    {
        return m_PiS;
    }
    const Eigen::MatrixXd SE2V2D::Pi0() const
    {
        return m_PiS;
    }
    const Eigen::MatrixXd SE2V2D::Pi0Dx() const
    {
        return m_Pi0Dx;
    }
    const Eigen::MatrixXd SE2V2D::Pi0Dy() const
    {
        return m_Pi0Dy;
    }
    const double SE2V2D::SM(const int alpha, const Eigen::Vector2d &pos) const
    {
        if (alpha == -1)
            return 0.0;
        if (alpha == 0)
            return 1.0;
        const Eigen::Vector2d scaledPos = ScaledCoord(pos);
        return _pow(scaledPos(0), mnl::PSpace2D::Exponent(alpha, 0)) * _pow(scaledPos(1), mnl::PSpace2D::Exponent(alpha, 1));
    }
    const double SE2V2D::SMIntegral(const int alpha) const
    {
        if (alpha < 0)
            return 0.0;
        return m_SMIntegrals[(size_t) alpha];
    }
    const double *SE2V2D::IntegralData() const
    {
        return m_SMIntegrals.data();
    }
    const Eigen::Vector2d SE2V2D::ScaledCoord(const Eigen::Vector2d &pos) const
    {
        return (pos - m_Centroid) * m_InvDiameter;
    }
    const std::vector<double> SE2V2D::ScaledMonomialIntegrals(const int maxOrder) const
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
    const Eigen::MatrixXd SE2V2D::D_Impl() const
    {
        const int nv = (int) m_Polygon.size();
        const bool convex = m_ConcavityFunction.empty();
        const int nInner = mnl::PSpace2D::SpaceDim(m_InnerOrder);
        const int nDof = m_Order * nv + (convex ? nInner : 2. * nInner);
        // The Serendipity projector has to supply moments up to order m_GradOrder - 1
        const int nk = mnl::PSpace2D::SpaceDim(m_GradOrder - 1); // m_GradOrder - 1 >= m_Order by

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

        for (size_t beta{}; beta < nInner; ++beta) 
            for (size_t alpha{}; alpha < nk; ++alpha) 
                for (const auto& [index, val] : m_ConcavityFunction)
                    D(m_Order * nv + nInner + beta, alpha) += SMIntegral(mnl::PSpace2D::Product(mnl::PSpace2D::Product(alpha, beta), index)) * val;

        return D;
    }
    const Eigen::MatrixXd SE2V2D::G0D_Impl() const
    {
        const int nl = mnl::PSpace2D::SpaceDim(m_GradOrder);
        Eigen::MatrixXd G0 = Eigen::MatrixXd::Zero(nl, nl);
        for (int alpha{}; alpha < nl; ++alpha) {
            G0(alpha, alpha) = SMIntegral(mnl::PSpace2D::Product(alpha, alpha));
            for (int beta = alpha + 1; beta < nl; ++beta)
                G0(alpha, beta) = G0(beta, alpha) = SMIntegral(mnl::PSpace2D::Product(beta, alpha));   
        }
        return G0;
    }
    const std::array<Eigen::MatrixXd, 2> SE2V2D::B0D_Impl() const
    {
        const int nkD = mnl::PSpace2D::SpaceDim(m_GradOrder);
        const int nv = (int) m_Polygon.size();
        const int nInner = mnl::PSpace2D::SpaceDim(m_InnerOrder);
        const int nDof = m_Order * nv + nInner; // Concavity?
        Eigen::MatrixXd B0Dx = Eigen::MatrixXd::Zero(nkD, nDof);
        Eigen::MatrixXd B0Dy = Eigen::MatrixXd::Zero(nkD, nDof);


        // Boundary term
        const auto lobattoPositions = LobattoNodePositions(m_Order);
        const size_t nEdgePoints = lobattoPositions.size() - 2;
        const int maxN = ceil(0.5 * (m_Order + m_GradOrder + 1));
        for (int n = 1; n <= maxN; ++n) {
            const auto quadrature = mnl::GaussLegendreRN(n);
            const int startAlpha = mnl::PSpace2D::SpaceDim(2 * (n - 1) - 1),
                      endAlpha = std::min(mnl::PSpace2D::SpaceDim(2 * n - 1), nkD);
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

        // Contribution of the Serendipity projection
        for (int beta = nInner, nLaplacianSpace = mnl::PSpace2D::SpaceDim(m_GradOrder - 1); beta < nLaplacianSpace; ++beta) { // Works better this way
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
        
        if (m_InnerOrder < 0)
            return { B0Dx, B0Dy };

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

        

        return { B0Dx, B0Dy };
    }
    const Eigen::MatrixXd SE2V2D::Pi0_Impl() const
    {
        const Eigen::MatrixXd D = D_Impl().topRows(mnl::PSpace2D::SpaceDim(m_Order));
        const Eigen::MatrixXd DT = D.transpose();
        return (DT * D).ldlt().solve(DT);
    }
    int SE2V2D::CheckConcavity()
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
    int SE2V2D::StabFreeGradOrderEstimate(int nVertices, int order) const
    {
        return int(0.5 * ((int)m_Polygon.size() + 2 * m_Order - 5) + 1. - 1e-8);
    }
    int SE2V2D::StabFreeGradOrder() const
    {
        // m_Order + 1 is the one with no additional degrees of freedom
        // m_InnerOrder = std::max(-1, m_GradOrder - 1 - eta); - SE2V2D
        // m_InnerOrder = std::max(-1, m_Order - eta); - SV2D
        return std::max(m_Order + 1, StabFreeGradOrderEstimate((int) m_Polygon.size(), m_Order));
    }
}
#include "V3D.h"
#include "ptp.h"
#include "mnl/include/mnl.hpp"
#include "mnl/include/glq.hpp"
#include "Eigen/Cholesky"

namespace ves {
    V3D::V3D(const std::vector<Eigen::Vector3d>& vertices,
            const std::vector<std::vector<size_t>>& faces,
            const int order,
            const int maxMonomialOrder = -1)
        : m_Order{order}, m_Faces{faces}, m_Vertices{vertices}
        { 
            m_SMIntegrals = ScaledMonomialIntegrals(std::max(maxMonomialOrder, 2 * (m_Order - 1)));
            Init();
        }

        const double V3D::SM(const int alpha, const Eigen::Vector3d &pos) const
        {
            double SM = 1.0;
            if (alpha == 0)
                return SM;

            const auto scaledCoord = ScaledCoord(pos);
            int sumExp = 0;
            for (int x = 0; x < 3; ++x)
            {
                const int exp = mnl::PSpace3D::Exponent(alpha, x);
                if (exp == 0)
                    continue;
                SM *= pow(scaledCoord(x), exp);
                sumExp += exp;
            }
            SM *= pow(m_InvDiameter, sumExp);

            return SM;
        }

        const double V3D::SMIntegral(const int alpha) const
        {
            return m_SMIntegrals[alpha];
        }

        const double *V3D::IntegralData() const
        {
            return m_SMIntegrals.data();
        }

        void V3D::Init()
        {
            // Parse edges (give them a known order that can be recalled)
            if (m_Order > 1)
                ParseEdges();

            // L2-Projector
            const Eigen::MatrixXd G0 = G0_Impl();
            const Eigen::MatrixXd B0 = B0_Impl();
            m_Pi0 = G0.ldlt().solve(B0);

            // Grad-Projector
            const Eigen::MatrixXd GGrad = GGrad_Impl();
            const Eigen::MatrixXd BGrad = BGrad_Impl();
            m_PiGrad = GGrad.ldlt().solve(BGrad);
        }

    const double V3D::Volume() const
    {
        return m_SMIntegrals[0];
    }
    const Eigen::Vector3d V3D::Centroid() const
    {
        return m_Centroid;
    }

    const Eigen::MatrixXd V3D::D() const
    {
        return D_Impl();
    }

    const Eigen::MatrixXd V3D::GGradTilde() const
    {
        return GGradTilde_Impl();
    }

    const Eigen::MatrixXd V3D::PiGrad() const
    {
        return m_PiGrad;
    }

    const Eigen::MatrixXd V3D::Pi0() const
    {
        return m_Pi0;
    }

    const Eigen::MatrixXd V3D::GGrad() const
    {
        return GGrad_Impl();
    }

    const Eigen::MatrixXd V3D::G0() const
    {
        return G0_Impl();
    }

    void V3D::ParseEdges()
    {
        const int pnv = static_cast<int>(m_Vertices.size());
        const int pnf = static_cast<int>(m_Faces.size());
        const int pne = pnv + pnf - 2;
        m_EdgeNodes.reserve(pne * (m_Order - 1));

        std::vector<EdgeCode> unorderedCodes;
        unorderedCodes.reserve((m_Order - 1) * pne);
        for (const auto& face : m_Faces) {
            const int nv = face.size();
            for (int v = 0; v < nv; ++v) {
                for (int inner = 0; inner < m_Order - 1; ++inner) {
                    unorderedCodes.push_back(GetEdgeCode(face[v], face[(v + 1) % nv], inner));
                }
            }
        }
        std::ranges::sort(unorderedCodes);
        const auto& ret = std::ranges::unique(unorderedCodes);
        unorderedCodes.erase(ret.begin(), ret.end());

        m_EdgeNodes = unorderedCodes;
    }

    const Eigen::Vector3d V3D::ScaledCoord(const Eigen::Vector3d &pos) const
    {
        return m_InvDiameter * (pos - m_Centroid);
    }

    const std::vector<double> V3D::ScaledMonomialIntegrals(const int maxOrder) const
    {
        std::vector<Eigen::Vector3d> scaledVertices;
        scaledVertices.reserve(m_Vertices.size());
        std::transform( m_Vertices.cbegin(), 
                        m_Vertices.cend(),
                        std::back_inserter(scaledVertices),
                        [this](const auto& pos) { return ScaledCoord(pos); } );
        return ptp::Polyhedron::MonomialIntegrals(scaledVertices, m_Faces, maxOrder);
    }


    const Eigen::MatrixXd V3D::GGradTilde_Impl() const
    {
        const int nk = mnl::PSpace3D::SpaceDim(m_Order);
        Eigen::MatrixXd GGT = Eigen::MatrixXd::Zero(nk, nk);
        for (int r = 1; r < nk; ++r) {
            for (int c = 1; c < nk; ++r) {
                for (int x = 0; x < 3; ++x) {
                    const int rexpx = mnl::PSpace3D::Exponent(r, x);
                    const int cexpx = mnl::PSpace3D::Exponent(c, x);
                    
                    // The remaining code is skippable if any of the monomials does not include the current variable.
                    if (rexpx * cexpx == 0)
                        continue;

                    const int rdx = mnl::PSpace3D::D(r, x);
                    const int cdx = mnl::PSpace3D::D(c, x);
                    GGT(r,c) += rexpx * cexpx * m_SMIntegrals[mnl::PSpace3D::Product(rdx, cdx)];
                }
            }
        }
        // Apply 1/h^2 factor.
        return GGT * m_InvDiameter * m_InvDiameter;
    }

    const Eigen::MatrixXd V3D::GGrad_Impl() const
    {
        Eigen::MatrixXd GGrad = GGradTilde_Impl();

        const int nk = static_cast<int>(GGrad.rows());
        if (m_Order == 1){
            // P0 operator
            const double nv = static_cast<double>(m_Vertices.size());
            GGrad(0,0) = nv;
            for (const auto& vertex : m_Vertices)
                for (int alpha = 1; alpha < nk; ++alpha)
                    GGrad(0, alpha) += SM(alpha, vertex);
            GGrad.row(0) /= nv;
            return GGrad;
        }

        for (int alpha = 0; alpha < nk; ++alpha)
	        GGrad(0, alpha) += m_SMIntegrals[alpha];
        GGrad.row(0) /= m_SMIntegrals[0];
        
        return GGrad;
    }
    const Eigen::MatrixXd V3D::G0_Impl() const
    {
        const int nk = mnl::PSpace3D::SpaceDim(m_Order);
        Eigen::MatrixXd G0 = Eigen::MatrixXd::Zero(nk, nk);
        for (int r = 0; r < nk; ++r)
            for (int c = 0; c < nk; ++c)
                G0(r, c) = m_SMIntegrals[mnl::PSpace3D::Product(r, c)];
        return G0;
    }
    const Eigen::Vector3d V3D::EdgeNodePosition(EdgeCode code) const
    {
        const int key = static_cast<int>(m_Vertices.size());
        
        const auto qNats = mnl::GaussLobattoR(2 * (m_Order + 1) - 3);

        const int end = code % key;
        code /= key;
        const int start = code % key;
        code /= key;
        const int innerPos = code % key;

        const double xi = qNats[innerPos + 1][0];
        return m_Vertices[start] * (1. - xi) + m_Vertices[end] * xi;
    }
    const V3D::EdgeCode V3D::GetEdgeCode(int start, int end, int innerPos) const
    {
        const int key = static_cast<int>(m_Vertices.size());
        bool ordered = start < end;
        return (end > start ? innerPos : m_Order - 2 - innerPos) * pow(key, 2) + (end > start ? start : end) * key + (end > start ? end : start);
    }
};
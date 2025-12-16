#include "V3D.h"

#include <unordered_set>
#include <algorithm>
#include <numeric>

#include <Eigen/Cholesky>

#include "ptp.h"
#include "mnl/include/mnl.hpp"
#include "mnl/include/glq.hpp"

#include "ves_internal.h"


namespace ves {
    V3D::V3D(const std::vector<V3D_Face*>& faceElements,
            const int order,
            const std::vector<size_t>& invertedFaces,
            const int maxMonomialOrder)
        : m_Order{order}, m_FaceElements{faceElements}
    {
        // Compute vector of vertices of the polyhedron from face elements
        SetupVertices();

        // Compute face connectivity with respect to this vector of vertices
        SetupFaces(invertedFaces);

        // Compute and store data about edge nodes, if needed
        if (m_Order > 1)
            ParseEdges();

        // Compute and store scaled monomial integrals
        m_SMIntegrals = ScaledMonomialIntegrals(std::max(maxMonomialOrder, 2 * (m_Order - 1)));

        // Compute and store projectors
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
        // L2-Projector
        const Eigen::MatrixXd G0 = G0_Impl();
        const Eigen::MatrixXd B0 = B0_Impl();
        m_Pi0 = G0.ldlt().solve(B0);

        // Grad-Projector
        const Eigen::MatrixXd GGrad = GGrad_Impl();
        const Eigen::MatrixXd BGrad = BGrad_Impl();
        m_PiGrad = GGrad.ldlt().solve(BGrad);
    }

    void V3D::SetupVertices()
    {
        // Find Vertex list
        std::unordered_set<Eigen::Vector3d> verticesSet;
        const size_t upperBoundNVertices = std::transform_reduce(m_FaceElements.cbegin(),
                                                       m_FaceElements.cend(), 
                                                       size_t{}, 
                                                       std::plus<>{}, 
                                                       [](const auto facePtr){ return facePtr->Vertices().size(); });
        verticesSet.reserve(upperBoundNVertices);
        for (const auto facePtr : m_FaceElements)
            for (auto& point : facePtr->Vertices())
                verticesSet.insert(point);

        m_Vertices.reserve(verticesSet.size());
        std::copy(verticesSet.cbegin(), verticesSet.cend(), std::back_inserter(m_Vertices));
        
    }

    void V3D::SetupFaces(const std::vector<size_t> &invertedFaces)
    {
        m_Faces.reserve(m_FaceElements.size());
        for (const auto facePtr : m_FaceElements){
            auto& faceIndices = m_Faces.emplace_back();
            const auto& vertices = facePtr->Vertices();
            faceIndices.reserve(vertices.size());
            for (const auto& pos : vertices) {
                const auto& it = std::find(m_Vertices.cbegin(), m_Vertices.cend(), pos);
                // Try and strip this if not debug
                if (it == m_Vertices.cend())
                    return;
                faceIndices.emplace_back(size_t{std::distance(m_Vertices.cbegin(), it)});
            }
        }
        for (const auto i : invertedFaces)
            std::reverse(m_Faces[i].begin(), m_Faces[i].end());
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

    const Eigen::MatrixXd V3D::D_Impl() const
    {
        const int vertexDofs = static_cast<int>(m_Vertices.size());
        const int edgeDofs = static_cast<int>(m_EdgeNodes.size());
        const int nInnerFace = mnl::PSpace2D::SpaceDim(m_Order - 2);
        const int faceDofs = nInnerFace * static_cast<int>(m_Faces.size());
        const int volDofs = mnl::PSpace3D::SpaceDim(m_Order - 2);
        const int ndof = vertexDofs + edgeDofs + faceDofs + volDofs;
        
        const int nk = mnl::PSpace3D::SpaceDim(m_Order);

        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(ndof, nk);
        
        // Vertex DOFs
        int i = 0;
        for (const auto& vertex : m_Vertices){
            for (int alpha = 0; alpha < nk; ++alpha){
                D(i, alpha) = SM(alpha, vertex);
            }
            ++i;
        }
        if (m_Order == 1)
            return D;
        
        // Edge DOFs
        for (const auto& edgecode : m_EdgeNodes){
            const auto pos = EdgeNodePosition(edgecode);
            for (int alpha = 0; alpha < nk; ++alpha){
                D(i, alpha) = SM(alpha, pos);
            }
            ++i;
        }

        // Face DOFs
        for (const auto& f : m_FaceElements) {
            for (size_t beta{}; beta < nInnerFace; ++beta) {
                for (int alpha = 0; alpha < nk; ++alpha) {
                    D(i, alpha) = f->MonomialMoment(beta, alpha, m_Centroid, m_InvDiameter);
                }
                ++i;
            }
        }

        // Volume DOFs
        for (int beta{}; beta < volDofs; ++beta) {
            for (int alpha{}; alpha < nk; ++alpha) {
                D(i, alpha) = m_SMIntegrals[mnl::PSpace3D::Product(alpha, beta)];
            }
            ++i;
        }
        
        return D;
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
        
        const int end = code % key;
        code /= key;
        const int start = code % key;
        code /= key;
        const int innerPos = code % key;
        
        const auto nats = EdgeNodePositions(m_Order);
        const double& xi = nats[innerPos];
        return m_Vertices[start] * (1. - xi) + m_Vertices[end] * xi;
    }
    const V3D::EdgeCode V3D::GetEdgeCode(int start, int end, int innerPos) const
    {
        const int key = static_cast<int>(m_Vertices.size());
        bool ordered = start < end;
        return (end > start ? innerPos : m_Order - 2 - innerPos) * pow(key, 2) + (end > start ? start : end) * key + (end > start ? end : start);
    }
    const std::vector<Eigen::Vector3d> &V3D_Face::Vertices()
    {
        return m_Vertices;
    }
};
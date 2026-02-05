#include "V3D.h"

#include <numeric>

#include <Eigen/Eigen/Geometry>

#include "ptp/ptp/src/ptp.h"
#include "mnl/include/mnl.hpp"
#include "mnl/include/gtq.hpp"

#include "ves_internal.h"


namespace ves {
    V3D::V3D(const std::vector<V3D_Face*>& faceElements,
            const int order,
            const std::vector<size_t>& invertedFaces,
            const int maxMonomialOrder)
        : m_Order{order}, 
          m_FaceElements{faceElements}, 
          m_Vertices{ComputeVertices()}, 
          m_Faces{ComputeFaces(invertedFaces)}, 
          m_EdgeNodes{ComputeEdges()},
          m_InvDiameter{ComputeInvDiameter()},
          m_Centroid{ComputeCentroid()},
          m_SMIntegrals{ScaledMonomialIntegrals(std::max(maxMonomialOrder, 2 * m_Order))},
          m_PiGrad{ComputePiGrad()},
          m_Pi0{ComputePi0()}
    {
    }

    const double V3D::SM(const int alpha, const Eigen::Vector3d &pos) const
    {
        double SM = 1.0;
        if (alpha == 0)
            return SM;

        const auto scaledCoord = ScaledCoord(pos);
        for (int x = 0; x < 3; ++x)
        {
            const int exp = mnl::PSpace3D::Exponent(alpha, x);
            if (exp == 0)
                continue;
            SM *= pow(scaledCoord(x), exp);
        }

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

    const std::vector<Eigen::Vector3d> V3D::ComputeVertices() const
    {
        const size_t upperBoundNVertices = std::transform_reduce(m_FaceElements.cbegin(),
                                                       m_FaceElements.cend(), 
                                                       size_t{}, 
                                                       std::plus<>{}, 
                                                       [](const auto facePtr){ return facePtr->Vertices().size(); });

        std::vector<Eigen::Vector3d> vertices;
        vertices.reserve(upperBoundNVertices);

        // Populating all vertices
        for (const auto facePtr : m_FaceElements)
            for (auto& point : facePtr->Vertices())
                vertices.emplace_back(point);

        // Removing duplicates
        std::sort(vertices.begin(), vertices.end(),
            [](const auto& one, const auto& other){ return one[0] + 1e1 * one[1] + 1e2 * one[2] < other[0] + 1e1 * other[1] + 1e2 * other[2]; }
        );
        
        auto last = std::unique(vertices.begin(), vertices.end(), [](const auto& one, const auto& other){ return (one - other).norm() == 0.0; });
        vertices.erase(last, vertices.end());

        return vertices;
    }

    const std::vector<std::vector<size_t>> V3D::ComputeFaces(const std::vector<size_t> &invertedFaces) const
    {
        std::vector<std::vector<size_t>> faces;
        faces.reserve(m_FaceElements.size());
        for (const auto facePtr : m_FaceElements){
            auto& faceIndices = faces.emplace_back();
            const auto& vertices = facePtr->Vertices();
            faceIndices.reserve(vertices.size());
            for (const auto& pos : vertices) {
                const auto& it = std::find(m_Vertices.cbegin(), m_Vertices.cend(), pos);
                // // Try and strip this if not debug
                // if (it == m_Vertices.cend())
                //     return {};
                faceIndices.emplace_back(static_cast<size_t>(std::distance(m_Vertices.cbegin(), it)));
            }
        }
        for (const auto i : invertedFaces)
            std::reverse(faces[i].begin(), faces[i].end());

        return faces;
    }

    const double V3D::Volume() const
    {
        return m_SMIntegrals[0];
    }
    const Eigen::Vector3d V3D::Centroid() const
    {
        return m_Centroid;
    }

    const double V3D::InverseDiameter() const
    {
        return m_InvDiameter;
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

    const Eigen::MatrixXd V3D::BGrad() const
    {
        return GGrad() * m_PiGrad;
    }

    const Eigen::MatrixXd V3D::B0() const
    {
        return B0_Impl();
    }

    const std::vector<V3D::EdgeCode> V3D::ComputeEdges() const
    {
        if (m_Order == 1)
            return {};

        const int pnv = static_cast<int>(m_Vertices.size());
        const int pnf = static_cast<int>(m_Faces.size());
        const int pne = pnv + pnf - 2;

        std::vector<V3D::EdgeCode> unorderedCodes;
        unorderedCodes.reserve((m_Order - 1) * pne);
        for (const auto& face : m_Faces) {
            const int nv = face.size();
            for (int v = 0; v < nv; ++v) {
                for (int inner = 0; inner < m_Order - 1; ++inner) {
                    unorderedCodes.push_back(GetEdgeCode(face[v], face[(v + 1) % nv], inner));
                }
            }
        }
        std::sort(unorderedCodes.begin(), unorderedCodes.end());
        auto ret = std::unique(unorderedCodes.begin(), unorderedCodes.end());
        unorderedCodes.erase(ret, unorderedCodes.end());

        return unorderedCodes;
    }

    const Eigen::Vector3d V3D::ComputeCentroid() const
    {
        const auto ints = ptp::Polyhedron::MonomialIntegrals(m_Vertices, m_Faces, 1);
        return Eigen::Vector3d(ints[1], ints[2], ints[3]) / ints[0];
    }

    const double V3D::ComputeInvDiameter() const
    {
        return 1. / ptp::Polyhedron::Diameter(m_Vertices);
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
                        
        std::vector<double> Integrals = ptp::Polyhedron::MonomialIntegrals(scaledVertices, m_Faces, maxOrder);
        
        // Adjust jacobian
        const double jacobian = _pow(ptp::Polyhedron::Diameter(m_Vertices), 3);
        std::transform(Integrals.begin(), Integrals.end(), Integrals.begin(), [&jacobian](double integral){ return integral * jacobian; });
        return Integrals;
    }

    const Eigen::MatrixXd V3D::ComputePiGrad() const
    {
        const Eigen::MatrixXd GGrad = GGrad_Impl();
        const Eigen::MatrixXd BGrad = BGrad_Impl();
        
        return GGrad.fullPivLu().solve(BGrad);
    }

    const Eigen::MatrixXd V3D::ComputePi0() const
    {
        const Eigen::MatrixXd G0 = G0_Impl();
        const Eigen::MatrixXd B0 = B0_Impl();
        
        return G0.fullPivLu().solve(B0);
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
            for (int beta{}; beta < nInnerFace; ++beta) {
                for (int alpha = 0; alpha < nk; ++alpha) {
                    D(i, alpha) = f->MonomialMoment(beta, alpha, m_Centroid, m_InvDiameter);
                }
                ++i;
            }
        }

        // Volume DOFs
        const double invvolume = 1. / SMIntegral(0);
        for (int beta{}; beta < volDofs; ++beta) {
            for (int alpha{}; alpha < nk; ++alpha) {
                D(i, alpha) = SMIntegral(mnl::PSpace3D::Product(alpha, beta)) * invvolume;
            }
            ++i;
        }
        
        return D;
    }

    const Eigen::MatrixXd V3D::GGradTilde_Impl() const
    {
        const int nuk = mnl::PSpace3D::SpaceDim(m_Order);
        Eigen::MatrixXd GGT = Eigen::MatrixXd::Zero(nuk, nuk);
        for (int r = 1; r < nuk; ++r) {
            for (int c = 1; c < nuk; ++c) {
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

        const int nuk = static_cast<int>(GGrad.rows());
        if (m_Order == 1){
            // P0 operator
            const double nv = static_cast<double>(m_Vertices.size());
            GGrad(0,0) = nv;
            for (const auto& vertex : m_Vertices)
                for (int alpha = 1; alpha < nuk; ++alpha)
                    GGrad(0, alpha) += SM(alpha, vertex);
            GGrad.row(0) /= nv;
            return GGrad;
        }

        for (int alpha = 0; alpha < nuk; ++alpha)
	        GGrad(0, alpha) += m_SMIntegrals[alpha];
        GGrad.row(0) /= m_SMIntegrals[0];
        
        return GGrad;
    }
    const Eigen::MatrixXd V3D::BGrad_Impl() const
    {
        const int nuk = mnl::PSpace3D::SpaceDim(m_Order);
        const int nukgrad = mnl::PSpace3D::SpaceDim(m_Order - 1);
        const int nuki = mnl::PSpace3D::SpaceDim(m_Order - 2);
        const int nki = mnl::PSpace2D::SpaceDim(m_Order - 2);
        const int nv = (int)m_Vertices.size();
        const int ne = (int)m_EdgeNodes.size();
        const int nf = (int)m_Faces.size();
        const int nEdgeNodes = m_Order - 1;
        const int nVE = nv + ne; // Vertex and Edge DOFs, ne already is multiplied by m_Order - 1;
        const int nB = nVE + nki * nf;
        const int ndof = nB + nuki;
        
        // Setting up monomial data
        Eigen::MatrixXi muExp = Eigen::MatrixXi::Zero(nuk, 3);
        Eigen::MatrixXi muD = Eigen::MatrixXi::Zero(nuk, 3);
        for (int alpha{}; alpha < nuk; ++alpha) {
            for (int j{}; j < 3; ++j) {
                muExp(alpha, j) = mnl::PSpace3D::Exponent(alpha, j);
                muD(alpha, j) = mnl::PSpace3D::D(alpha, j);
            }
        }

        // Setting up edge data
        std::unordered_map<EdgeCode, int> edgeIndex;
        {
            int i = 0;
            for (const auto& ec : m_EdgeNodes)
                edgeIndex[ec] = i++;
        }
        
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nuk, ndof);
        // First row uses P0
        if (m_Order == 1)
            B.row(0) = Eigen::VectorXd::Ones(ndof) / nv;
        else
            B.block(0, nB, nuki, nuki) = Eigen::MatrixXd::Identity(nuki, nuki);
        
        // Boundary integral
        for (size_t f{}; f < (size_t)nf; ++f) {
            // BF(alpha, j) = \int_F{\mu_\alpha * phi_j * d\sigma}
            // j - being the local index (face)

            const auto& face = *m_FaceElements[f];
            const auto BF = face.B0(m_Centroid, m_InvDiameter);
            const Eigen::Vector3d normal = face.Normal();

            const auto& faceIndices = m_Faces[f];
            const size_t nvF = faceIndices.size();
            const int nBF = m_Order * (int)nvF;
            const int ndofF = nBF + nki;

            // Inverter loop de k e alpha, e primeiro verificar se abs(normal[k]) > tol.
            constexpr double tol = 1e-8;
            for (int k{}; k < 3; ++k) { // Direction
                if (abs(normal[k]) < tol) continue;
                for (int alpha = 1; alpha < nuk; ++alpha) {
                    if (muD(alpha, k) == -1) continue; // Derivative is 0
                    const double aux = normal[k] * muExp(alpha, k) * m_InvDiameter;

                    // Face-Internal DOFs
                    for (int j = 0; j < nki; ++j)
                        B(alpha, nVE + (int)f * nki + j) += BF(muD(alpha, k), nBF + j) * aux;

                    if (alpha < nukgrad) continue; // Boundary of Face DOFs are zero for monomial moments that constitute internal face DOFs.

                    // Vertex and Edge DOFs
                    for (int j{}; j < nvF; ++j) { // Vertices
                        const int i = (int)faceIndices[j]; // Corresponding polyhedron index
                        B(alpha, i) += BF(muD(alpha, k), j) * aux;
                        for (int e{}; e < m_Order - 1; ++e) { // Edge nodes
                            const auto edgeCode = GetEdgeCode(i, (int)faceIndices[(j + 1) % nvF], e);
                            const int faceEdgeIndex = (int)nvF + nEdgeNodes * j + e;
                            B(alpha, nv + edgeIndex[edgeCode]) += BF(muD(alpha, k), faceEdgeIndex) * aux;
                        }
                    }
                }
            }
        }

        // Domain integral
        const double cst = Volume() * m_InvDiameter * m_InvDiameter;
        for (int alpha{}; alpha < nuki; ++alpha) {
            for (int k{}; k < 3; ++k) {
                const int beta = mnl::PSpace3D::AD(mnl::PSpace3D::AD(alpha, k), k);
                const int exp = mnl::PSpace3D::Exponent(alpha, k);
                const double c = double((exp + 1) * (exp + 2));
                B(beta, nB + alpha) -= cst * c;
            }
        }

        return B;
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
    const Eigen::MatrixXd V3D::B0_Impl() const
    {
        const int nuk = mnl::PSpace3D::SpaceDim(m_Order);
        const int nuki = mnl::PSpace3D::SpaceDim(m_Order - 2);
        const int nukiq = nuk - nuki;
        const int nki = mnl::PSpace2D::SpaceDim(m_Order - 2);
        const int nv = (int)m_Vertices.size();
        const int ne = (int)m_EdgeNodes.size();
        const int nf = (int)m_Faces.size();
        const int nEdgeNodes = m_Order - 1;
        const int nB = nv + ne + nki * nf;
        const int ndof = nB + nuki;
        const double volume = Volume();
        
        Eigen::MatrixXd PI = Eigen::MatrixXd::Zero(nukiq, nuk);
        for (int alpha{}; alpha < nukiq; ++alpha)
            for (int beta{}; beta < nuk; ++beta)
                PI(alpha, beta) = SMIntegral(mnl::PSpace3D::Product(nuki + alpha, beta));

        Eigen::MatrixXd B0 = Eigen::MatrixXd::Zero(nuk, ndof);
        // Internal DOFs contribution
        B0.block(0, nB, nuki, nuki) = Eigen::MatrixXd::Identity(nuki, nuki) * volume;
        B0.block(nuki, 0, nukiq, ndof) = PI * m_PiGrad;

        return B0;
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

    V3D_Face::V3D_Face(const std::vector<Eigen::Vector3d> &vertices, const int order)
        : m_Order(order), m_Vertices(vertices), m_Centroid(ComputeCentroid()), m_ChangeBasis(ComputeChangeBasis()), m_LocalSpace(LocalVertices(), m_Order)
    {
    }

    const std::vector<Eigen::Vector3d> &V3D_Face::Vertices() const
    {
        return m_Vertices;
    }

    const Eigen::Vector3d V3D_Face::Normal() const
    {
        return ptp::Polygon3D::Normal(m_Vertices);
    }

    const double V3D_Face::MonomialMoment(const int beta2D, const int alpha3D, const Eigen::Vector3d &polyhedronCentroid, const double polyhedronInvDiameter) const
    {
        const auto triangulation = ptp::Polygon2D::Triangulation(LocalVertices());
        const double area = m_LocalSpace.SMIntegral(0);
        double moment = 0.0;
        for (const auto& triangleIndices : triangulation) {
            const std::vector<Eigen::Vector3d> triangle{
                m_Vertices[triangleIndices[0]],
                m_Vertices[triangleIndices[1]],
                m_Vertices[triangleIndices[2]]
            };
            const double triArea = ptp::Polygon3D::MonomialIntegrals(triangle, 0)[0];
            const int k2D = mnl::PSpace2D::MonOrder(beta2D);
            const int k3D = mnl::PSpace3D::MonOrder(alpha3D);
            const auto quadrature = mnl::GaussLegendreTriangle(k2D + k3D);
            for (const auto& qData : quadrature) {
                const double& xi0 = qData[0],
                            & xi1 = qData[1],
                              xi2 = 1. - xi0 - xi1,
                            & w = qData[2];
                const Eigen::Vector3d qPos = xi0 * triangle[0] + xi1 * triangle[1] + xi2 * triangle[2];
                const double muAlpha = SM3D(alpha3D, qPos, polyhedronCentroid, polyhedronInvDiameter);
                const double mBeta = m_LocalSpace.SM(beta2D, LocalCoordinate(qPos));
                moment += muAlpha * mBeta * w * triArea;
            }
        }
        moment /= area;
        return moment;
    }

    const double V3D_Face::SM3D(int alpha, const Eigen::Vector3d &pos, const Eigen::Vector3d &polyhedronCentroid, const double polyhedronInvDiameter) const
    {
        if (alpha == -1)
            return 0.0;

        double SM = 1.0;
        if (alpha == 0)
            return SM;

        const auto scaledCoord = (pos - polyhedronCentroid) * polyhedronInvDiameter;
        for (int x = 0; x < 3; ++x)
            SM *= pow(scaledCoord(x), mnl::PSpace3D::Exponent(alpha, x));
        return SM;
    }

    const Eigen::MatrixXd V3D_Face::DM(const Eigen::Vector3d &polyhedronCentroid, const double polyhedronInvDiameter) const
    {
        const int nv = m_Vertices.size();
        const int nB = nv * m_Order;
        const int nki = mnl::PSpace2D::SpaceDim(m_Order - 2);
        const int ndof = nB + nki;
        const int nukgrad = mnl::PSpace3D::SpaceDim(m_Order - 1);
        const int nuki = mnl::PSpace3D::SpaceDim(m_Order - 2);
        const int nukiq = nukgrad - nuki;
        Eigen::MatrixXd DF = Eigen::MatrixXd::Zero(ndof, nukgrad);

        // Boundary DOFs and additional setup
        for (int v = 0; v < nv; ++v) {
            const Eigen::Vector3d& start = m_Vertices[(size_t) v];
            // Vertex DOFs
            for (int alpha = 0; alpha < nukgrad; ++alpha) {
                DF(v, alpha) = SM3D(alpha, start, polyhedronCentroid, polyhedronInvDiameter);
            }
        }
        
        if (m_Order == 1)
            return DF;
        
        const auto edgePoints = EdgeNodePositions(m_Order);
        const int nEdgePoints = (int) edgePoints.size();
        const auto localVertices = LocalVertices();

        Eigen::Vector3d edgePoint = Eigen::Vector3d::Zero();
        for (int v = 0; v < nv; ++v) {
            const int next = (v + 1) % nv;
            const Eigen::Vector3d& start = m_Vertices[(size_t) v];
            const Eigen::Vector3d& end = m_Vertices[(size_t)next];
            
            // Edge DOFs
            for (size_t e = 0; e < nEdgePoints; ++e) {
                const int index = nv + nEdgePoints * v + e;
                const double& xi = edgePoints[e];
                edgePoint = (1. - xi) * start + xi * end;
                for (int alpha = 0; alpha < nukgrad; ++alpha) {
                    DF(index, alpha) = SM3D(alpha, edgePoint, polyhedronCentroid, polyhedronInvDiameter);
                }
            }
        }
        
        // Internal DOFs, need to integrate.
        const auto triangulation = ptp::Polygon2D::Triangulation(localVertices);
        const double area = m_LocalSpace.SMIntegral(0);
        for (const auto& triangleIndices : triangulation){
            const std::vector<Eigen::Vector3d> triangle{
                m_Vertices[triangleIndices[0]],
                m_Vertices[triangleIndices[1]],
                m_Vertices[triangleIndices[2]]
            };
            const double triArea = ptp::Polygon3D::MonomialIntegrals(triangle, 0)[0];
            for (int k2D = 0; k2D <= m_Order - 2; ++k2D) {
                const int startBeta = mnl::PSpace2D::SpaceDim(k2D - 1);
                const int maxBeta = mnl::PSpace2D::SpaceDim(k2D);
                for (int k3D = 0; k3D < m_Order; ++k3D) {
                    const int startAlpha = mnl::PSpace3D::SpaceDim(k3D - 1);
                    const int maxAlpha = mnl::PSpace3D::SpaceDim(k3D);
                    const auto quadrature = mnl::GaussLegendreTriangle(k2D + k3D);
                    for (const auto& qData : quadrature){
                        const double& xi0 = qData[0],
                                    & xi1 = qData[1],
                                      xi2 = 1. - xi0 - xi1,
                                    & w   = qData[2];
                        const Eigen::Vector3d qPos = xi0 * triangle[0] + xi1 * triangle[1] + xi2 * triangle[2];
                        for (int alpha = startAlpha; alpha < maxAlpha; ++alpha) {
                            const double muAlpha = SM3D(alpha, qPos, polyhedronCentroid, polyhedronInvDiameter);
                            for (int beta = startBeta; beta < maxBeta; ++beta) {
                                const double mBeta = m_LocalSpace.SM(beta, LocalCoordinate(qPos));
                                DF(nB + beta, alpha) += muAlpha * mBeta * w * triArea;
                            }
                        }
                    }
                }
            }
        }
        DF.block(nB, 0, nki, nukgrad) /= area;
        return DF;
    }

    const Eigen::MatrixXd V3D_Face::B0(const Eigen::Vector3d &polyhedronCentroid, const double polyhedronInvDiameter) const
    {
        const int nB = m_Order * m_Vertices.size();
        const int nk = mnl::PSpace2D::SpaceDim(m_Order);
        const int nkgrad = mnl::PSpace2D::SpaceDim(m_Order - 1);
        const int nki = mnl::PSpace2D::SpaceDim(m_Order - 2);
        const int nkiq = nkgrad - nki; // Dimension of the quotient space P_k-1/P_ki
        const int nukgrad = mnl::PSpace3D::SpaceDim(m_Order - 1);
        const int nuki = mnl::PSpace3D::SpaceDim(m_Order - 2);
        const int nukiq = nukgrad - nuki;
        const int ndof = nB + nki;

        Eigen::MatrixXd BF = Eigen::MatrixXd::Zero(nukgrad, ndof);

        // The projection of a polynomial is itself
        const Eigen::MatrixXd Pi0 = m_LocalSpace.Pi0();
        // DM(i, \alpha) = DOF_i(\mu_\alpha).
        // MF(\beta,\alpha) - m_\beta component of \mu_\alpha.
        const Eigen::MatrixXd MF = Pi0.block(0, 0, nkgrad, ndof) * DM(polyhedronCentroid, polyhedronInvDiameter);

        // Contribution of monomials associated with internal DOFs
        const double area = m_LocalSpace.SMIntegral(0);
        BF.block(0, nB, nukgrad, nki) += area * MF.block(0, 0, nki, nukgrad).transpose();

        Eigen::MatrixXd ProductIntegrals = Eigen::MatrixXd::Zero(nkiq, nk);
        for (int beta = nki; beta < nkgrad; ++beta)
            for (int gamma = 0; gamma < nk; ++gamma)
                ProductIntegrals(beta - nki, gamma) = m_LocalSpace.SMIntegral(mnl::PSpace2D::Product(beta, gamma));

        BF.block(nuki, 0, nukiq, ndof) += MF.block(nki, nuki, nkiq, nukiq).transpose() * ProductIntegrals * Pi0;

        return BF;
    }

    const std::vector<Eigen::Vector2d> V3D_Face::LocalVertices() const
    {
        std::vector<Eigen::Vector2d> out;
        out.reserve(m_Vertices.size());
        for (const auto& vertex : m_Vertices)
            out.emplace_back(LocalCoordinate(vertex));
        return out;
    }

    const Eigen::Vector2d V3D_Face::LocalCoordinate(const Eigen::Vector3d &globalCoord) const
    {
        return m_ChangeBasis * (globalCoord - m_Centroid);
    }

    const Eigen::Matrix<double, 2, 3> V3D_Face::ComputeChangeBasis() const
    {
        Eigen::Matrix<double, 2, 3> out = Eigen::Matrix<double, 2, 3>::Zero();
        const Eigen::Vector3d e0 = (m_Vertices[0] - m_Centroid).normalized(),
        e2 = ptp::Polygon3D::Normal(m_Vertices);
        const Eigen::Vector3d e1 = e2.cross(e0);

        out.row(0) = e0;
        out.row(1) = e1;
        return out;
    }

    const Eigen::Vector3d V3D_Face::ComputeCentroid() const
    {
        const auto integrals = ptp::Polygon3D::MonomialIntegrals(m_Vertices, 1);
        return {
            integrals[1] / integrals[0],
            integrals[2] / integrals[0],
            integrals[3] / integrals[0]
        };
    }
};
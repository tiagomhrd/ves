#pragma once
#ifndef VES_V3D
#define VES_V3D 
#include <Eigen/Eigen/Core>
#include "V2D.h"

namespace ves {
    class V3D_Face;

    /*  V3D - 3D Virtual Element in Modified Formulation
    
    This class provides basic structures for the Modified VEM formulation 
    (see https://www.sciencedirect.com/science/article/pii/S0898122113003179#s000015)
    REVIEW
    The constructor takes:
    polyhedron                  - Vector of vertices and faces defining the polyhedron
    order                       - Order of the element
    (optional) maxMonomialOrder - Order up to which integrals of Scaled Monomials are computed
    */
    class V3D {
    public:
        V3D(const std::vector<V3D_Face*>& faceElements,
            const int order,
            const std::vector<size_t>& invertedFaces = {},
            const int maxMonomialOrder = -1);

        // Auxiliar geometry functions
        const double            Volume() const;
        const Eigen::Vector3d   Centroid() const;
        const double            InverseDiameter() const;

        // Main VEM structures for this formulation
        const Eigen::MatrixXd   D() const;
        const Eigen::MatrixXd   GGradTilde() const;
        const Eigen::MatrixXd   PiGrad() const;
        const Eigen::MatrixXd   Pi0() const;

        // Auxiliar VEM structures that may come in hand
        const Eigen::MatrixXd   GGrad() const;
        const Eigen::MatrixXd   G0() const;
        const Eigen::MatrixXd   BGrad() const;
        const Eigen::MatrixXd   B0() const;

        // Scaled Monomial functions
        const double            SM(const int alpha, const Eigen::Vector3d& pos) const; // Scaled monomial of index alpha evaluated at position pos
        const double            SMIntegral(const int alpha) const;                     // Integral of scaled monomial of index alpha over the polyhedron
        const double*           IntegralData() const;

    protected:
        const std::vector<Eigen::Vector3d> ComputeVertices() const;
        const std::vector<std::vector<size_t>> ComputeFaces(const std::vector<size_t>& invertedFaces) const;
        
        using EdgeCode = int;
        const Eigen::Vector3d EdgeNodePosition(EdgeCode code) const;
        const EdgeCode GetEdgeCode(int start, int end, int innerPos) const;
        const std::vector<EdgeCode> ComputeEdges() const;

        const Eigen::Vector3d ComputeCentroid() const;
        const double ComputeInvDiameter() const;
        const std::vector<double> ScaledMonomialIntegrals(const int maxOrder) const;

        const Eigen::MatrixXd ComputePiGrad() const;
        const Eigen::MatrixXd ComputePi0() const;

        const Eigen::Vector3d ScaledCoord(const Eigen::Vector3d& pos) const;

        const Eigen::MatrixXd D_Impl() const;

        const Eigen::MatrixXd GGradTilde_Impl() const;
        const Eigen::MatrixXd GGrad_Impl() const;
        const Eigen::MatrixXd BGrad_Impl() const;

        const Eigen::MatrixXd G0_Impl() const;
        const Eigen::MatrixXd B0_Impl() const;

    protected:
        const int m_Order;
        // Polyhedron representation
        const std::vector<V3D_Face*> m_FaceElements;
        const std::vector<Eigen::Vector3d> m_Vertices;
        const std::vector<std::vector<size_t>> m_Faces; // Each face is defined by a vector of vertex indices
        
        // Internal representation of edges
        const std::vector<EdgeCode> m_EdgeNodes; // Each entry encodes start and end point and innerPosition
        
        // Scaled monomial related storage
        const double m_InvDiameter;
        const Eigen::Vector3d m_Centroid;
        const std::vector<double> m_SMIntegrals;

        // Projector storage
        const Eigen::MatrixXd m_PiGrad;
        const Eigen::MatrixXd m_Pi0;
    };

    /*  V3D_Face - Face of a 3D Virtual Element
     
     This class provides basic structures for a face of a 3D Virtual Element
     to allow computation of face integrals and other face-related operations.
     
     The constructor takes:
     vertices    - Vector of counter-clockwise ordered vertices defining the face
     order       - Order of the face
    */
    class V3D_Face {
    public:
        V3D_Face(const std::vector<Eigen::Vector3d>& vertices, const int order);

        const std::vector<Eigen::Vector3d>& Vertices() const;
        const Eigen::Vector3d Normal() const;

        // Required for face DOFs in D
        const double MonomialMoment(const int beta2D, const int alpha3D, const Eigen::Vector3d& polyhedronCentroid, const double polyhedronInvDiameter) const;

        /* 
            B0(alpha, i) = \int_F{\mu_\alpha\phi_i d\sigma}
            * \mu_\alpha    Scaled monomial of index alpha of the polyhedron the current object is a face of.
            * \phi_i        Ansatz function of index i of the face.
            * polyhedronCentroid and polyhedronInvDiameter are required to compute \mu_\alpha
            
            This assumes \alpha < dim P_{m_Order - 1}.
        */
        const Eigen::MatrixXd B0(const Eigen::Vector3d& polyhedronCentroid, const double polyhedronInvDiameter) const;
        const Eigen::MatrixXd DM(const Eigen::Vector3d& polyhedronCentroid, const double polyhedronInvDiameter) const;

    protected:
        // Initialization
        const Eigen::Vector3d ComputeCentroid() const;
        const Eigen::Matrix<double, 2, 3> ComputeChangeBasis() const;
        const Eigen::VectorXd GaussLobattoWeightVector() const;
        const std::vector<Eigen::Vector2d> LocalVertices() const;
        const Eigen::Vector2d LocalCoordinate(const Eigen::Vector3d& globalCoord) const;

        // Computation
        const double SM3D(int alpha, const Eigen::Vector3d& pos, const Eigen::Vector3d& polyhedronCentroid, const double polyhedronInvDiameter) const;

    protected:
        int m_Order;
        std::vector<Eigen::Vector3d> m_Vertices;
        const Eigen::Vector3d m_Centroid;
        const Eigen::Matrix<double, 2, 3> m_ChangeBasis;
        const Eigen::VectorXd m_BoundaryIntegrationWeights;
        const V2D m_LocalSpace;
    };
}

#endif // VES_V3D
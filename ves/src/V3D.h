#pragma once
#ifndef VES_V3D
#define VES_V3D 
#include <Eigen/Eigen/Core>

namespace ves {
    /*  V3D - 3D Virtual Element in Modified Formulation
    
    This class provides basic structures for the Modified VEM formulation 
    (see https://www.sciencedirect.com/science/article/pii/S0898122113003179#s000015)
    
    The constructor takes:
    polyhedron                  - Vector of vertices and faces defining the polyhedron
    order                       - Order of the element
    (optional) maxMonomialOrder - Order up to which integrals of Scaled Monomials are computed
    */
    class V3D {
    public:
        V3D(const std::vector<Eigen::Vector3d>& vertices,
            const std::vector<std::vector<size_t>>& faces,
            const int order,
            const int maxMonomialOrder = -1);

        // Auxiliar geometry functions
        const double            Volume() const;
        const Eigen::Vector3d   Centroid() const;

        // Main VEM structures for this formulation
        const Eigen::MatrixXd   D() const;
        const Eigen::MatrixXd   GGradTilde() const;
        const Eigen::MatrixXd   PiGrad() const;
        const Eigen::MatrixXd   Pi0() const;

        // Auxiliar VEM structures that may come in hand
        const Eigen::MatrixXd   GGrad() const;
        const Eigen::MatrixXd   G0() const;
        const Eigen::MatrixXd   BGrad() const;

        // Scaled Monomial functions
        const double            SM(const int alpha, const Eigen::Vector3d& pos) const; // Scaled monomial of index alpha evaluated at position pos
        const double            SMIntegral(const int alpha) const;                     // Integral of scaled monomial of index alpha over the polyhedron
        const double*           IntegralData() const;

    protected:
        void Init();
        void ParseEdges();

        const Eigen::Vector3d ScaledCoord(const Eigen::Vector3d& pos) const;
        const std::vector<double> ScaledMonomialIntegrals(const int maxOrder) const;

        const Eigen::MatrixXd D_Impl() const;

        const Eigen::MatrixXd GGradTilde_Impl() const;
        const Eigen::MatrixXd GGrad_Impl() const;
        const Eigen::MatrixXd BGrad_Impl() const;

        const Eigen::MatrixXd G0_Impl() const;
        const Eigen::MatrixXd B0_Impl() const;

    protected:
        using EdgeCode = int;
        const Eigen::Vector3d EdgeNodePosition(EdgeCode code) const;
        const EdgeCode GetEdgeCode(int start, int end, int innerPos) const;

    protected:
        int m_Order;
        // Polyhedron representation
        std::vector<Eigen::Vector3d> m_Vertices;
        std::vector<std::vector<size_t>> m_Faces; // Each face is defined by a vector of vertex indices
        
        // Internal representation of edges
        std::vector<EdgeCode> m_EdgeNodes; // Each entry encodes start and end point and innerPosition
        
        // Scaled monomial related storage
        std::vector<double> m_SMIntegrals;
        Eigen::Vector3d m_Centroid;
        double m_InvDiameter;

        // Projector storage
        Eigen::MatrixXd m_PiGrad, m_Pi0;
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

    protected:
        int m_Order;
        std::vector<Eigen::Vector3d> m_Vertices;

    };
}

#endif // VES_V3D
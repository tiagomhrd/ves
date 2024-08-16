#pragma once
#ifndef VES_V2D
#define VES_V2D
#include <Eigen/Eigen/Core>

namespace ves {
    class V2D {
    public:
        /*  V2D - 2D Virtual Element in Modified Formulation

            This class provides basic structures for the Modified VEM formulation 
            (see https://www.sciencedirect.com/science/article/pii/S0898122113003179#s000015)

            The constructor takes:
            polygon                     - Vector of counter-clockwise ordered vertices
            order                       - Order of the element
            (optional) maxMonomialOrder - Order up to which integrals of Scaled Monomials are computed
        */
        V2D(const std::vector<Eigen::Vector2d>& polygon, 
            const int order, 
            const int maxMonomialOrder = -1);

        // Auxiliar geometry functions
        const double            Area() const;
        const Eigen::Vector2d   Centroid() const;

        // Main VEM structures for this formulation
        const Eigen::MatrixXd   D() const { return D_Impl(); }
        const Eigen::MatrixXd   GGradTilde() const { return GGradTilde_Impl(); }
        const Eigen::MatrixXd   PiGrad() const;
        const Eigen::MatrixXd   Pi0() const;

        // Auxiliar VEM structures that may come in hand
        const Eigen::MatrixXd   GGrad() const { return GGrad_Impl(); }
        const Eigen::MatrixXd   G0() const { return G0_Impl(); }
        const Eigen::MatrixXd   BGrad() const { return BGrad_Impl(); }
        
        // Scaled Monomial functions
        const double            SM(const int alpha, const Eigen::Vector2d& pos) const; // Scaled monomial of index alpha evaluated at position pos
        const double            SMIntegral(const int alpha) const;                     // Integral of scaled monomial of index alpha over the polygon
        const double*           IntegralData() const;

    protected:
        void Init();

        const Eigen::Vector2d ScaledCoord(const Eigen::Vector2d& pos) const;
        const std::vector<double> ScaledMonomialIntegrals(const int maxOrder) const;

        const Eigen::MatrixXd D_Impl() const;

        const Eigen::MatrixXd GGradTilde_Impl() const;
        const Eigen::MatrixXd GGrad_Impl() const;
        const Eigen::MatrixXd BGrad_Impl() const;

        const Eigen::MatrixXd G0_Impl() const;
        const Eigen::MatrixXd B0_Impl() const;

    protected:
        int m_Order;
        std::vector<Eigen::Vector2d> m_Polygon;

        // Scaled monomial related storage
        std::vector<double> m_SMIntegrals;
        Eigen::Vector2d m_Centroid;
        double m_InvDiameter;

        // Projector storage
        Eigen::MatrixXd m_PiGrad, m_Pi0;
    };
}
#endif
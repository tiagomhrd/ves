#pragma once
#ifndef VES_SV2D
#define VES_SV2D
#include <Eigen/Eigen/Core>
#include <unordered_map>

namespace ves {
    class SV2D {
    public:
        SV2D(const std::vector<Eigen::Vector2d>& polygon, 
             const int order, 
             const int maxMonomialOrder = -1);

        // Auxiliar geometry functions
        const double            Area() const;
        const Eigen::Vector2d   Centroid() const;

        // Main VEM structures for this formulation
        const int InnerOrder() const;

        const Eigen::MatrixXd   D() const;
        
        const Eigen::MatrixXd   PiS() const;
        const Eigen::MatrixXd   Pi0() const;
        const Eigen::MatrixXd   Pi0Dx() const;
        const Eigen::MatrixXd   Pi0Dy() const;

        // Scaled Monomial functions
        const double            SM(const int alpha, const Eigen::Vector2d& pos) const; // Scaled monomial of index alpha evaluated at position pos
        const double            SMIntegral(const int alpha) const;                     // Integral of scaled monomial of index alpha over the polygon
        const double*           IntegralData() const;

    protected:
        void Init();

        const Eigen::Vector2d ScaledCoord(const Eigen::Vector2d& pos) const;
        const std::vector<double> ScaledMonomialIntegrals(const int maxOrder) const;

        const Eigen::MatrixXd D_Impl() const;
        const Eigen::MatrixXd G0_Impl() const;
        const std::array<Eigen::MatrixXd, 2> B0D_Impl() const;

        int CheckConcavity();

    protected:
        int m_Order;
        int m_InnerOrder; // Order up to which internal degrees of freedom have to exist.
        std::vector<Eigen::Vector2d> m_Polygon;

        // Scaled monomial related storage
        std::vector<double> m_SMIntegrals;
        Eigen::Vector2d m_Centroid;
        double m_InvDiameter;

        // Polynomial resulting from the product of the unique boundary support lines that define the concavities
        std::unordered_map<int, double> m_ConcavityFunction;

        // Projector storage
        Eigen::MatrixXd m_Pi0Dx, m_Pi0Dy, m_PiS;
    };
}

#endif;
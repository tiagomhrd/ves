#include "V3D.h"
#include "ptp.h"

namespace ves {
    const double V3D::Volume() const
    {
        return m_SMIntegrals[0];
    }
    const Eigen::Vector3d V3D::Centroid() const
    {
        return m_Centroid;
    }
};
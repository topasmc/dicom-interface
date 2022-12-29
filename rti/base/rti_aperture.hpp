#ifndef RTI_APERTURE_H
#define RTI_APERTURE_H

/// \file
///
/// Geometry model for an aperture

#include <rti/base/rti_geometry.hpp>

namespace rti{

/// \class aperture
/// \brief A aperture class 

class aperture : public geometry{
public:
    /// whehter aperture shape is box or cylinder type
    bool is_rectangle ;

    /// aperture dimension, (Lx, Ly, Lz) for box or (R, H, dummy) for cylinder
    const rti::vec3<float> volume;
    
    /// X-Y opening points. Divergence is not considered.
    const std::vector<std::vector<std::array<float,2>>> block_data;

public:

    /// Creates an aperture
    /// \param xypts a set of (x,y) points. 
    /// \param v a volume, e.g., (Lx,Ly,Lz) or (R, thickness, ignored). 
    /// \param p a center position of the aperture. 
    /// \param is_rect tells whether aperture is box shape or cylinder shape.
    aperture(
        std::vector<std::vector<std::array<float,2>>> xypts,
        rti::vec3<float>& v,
        rti::vec3<float>& p,
        rti::mat3x3<float>& r,
        bool is_rect = true) 
    : geometry (p, r, rti::geometry_type::BLOCK),
      is_rectangle(is_rect),
      volume(v),
      block_data(xypts)
    {;}

    /// Creates a copy from an existing aperture
    aperture(
        const rti::aperture& rhs)
    : geometry(rhs.pos, rhs.rot, rhs.geotype),
      is_rectangle(rhs.is_rectangle),
      volume(rhs.volume),
      block_data(rhs.block_data)
    {;}

    /// Destructor
    ~aperture(){;}

    /// Assignment operator
    const aperture&
    operator= 
    (const rti::aperture& rhs)
    {
        return rhs;
    }

};


}
#endif

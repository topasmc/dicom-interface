#ifndef RTI_RANGESHIFTER_H
#define RTI_RANGESHIFTER_H

/// \file
///
/// RT-Ion geometry for rangeshifter

#include <rti/base/rti_geometry.hpp>

namespace rti{

/// \class rangeshifter
/// supports box and cylinder type rangeshifter
/// \note only one material
class rangeshifter : public geometry{

public:

    const bool        is_rectangle; ///< type of volume
    const rti::vec3<float> volume;  ///< volume dimension

    /// Constructor
    rangeshifter(
        rti::vec3<float>& v,    ///< x,y,z or r, theta, thickness
        rti::vec3<float>& p,    ///< position
        rti::mat3x3<float>& r,  ///< rotation matrix
        bool is_rect=true) 
    : geometry(p, r, rti::geometry_type::RANGESHIFTER),
      is_rectangle(is_rect),
      volume(v)
    {;}

    /// Copy constructor
    rangeshifter(const rti::rangeshifter& rhs)
    : geometry(rhs.pos, rhs.rot, rhs.geotype),
      is_rectangle(rhs.is_rectangle),
      volume(rhs.volume)
    {;}

    /// Destructor
    ~rangeshifter(){;}

    /// Assignment operator
    const rangeshifter&
    operator= 
    (const rti::rangeshifter& rhs)
    {
        return rhs;
    }

};

}

#endif

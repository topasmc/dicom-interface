#ifndef RTI_BEAMLINE_H
#define RTI_BEAMLINE_H

/// \file
///
/// A beamline is a collection of physical components in the beamline

#include <rti/base/rti_coordinate_transform.hpp>
#include <rti/base/rti_geometry.hpp>

namespace rti{

/// \class beamline
/// 
/// \tparam T type of units
/// \note we may not need this class if this is a just container for geometries.
/// Will be determined later.
template <typename T>
class beamline {
public:

protected:

    /// A container for components
    std::vector<rti::geometry*> geometries_;

    /// A coordinate transform. 
    /// \note not sure we need this in future.
    rti::coordinate_transform<T> p_coord_;

public:

    /// An empty constructor
    beamline(){;}

    /// Clear all geometries in the vector
    ~beamline(){
        
        geometries_.clear();
    }

    /// Add new geometry to the container.
    void
    append_geometry(
        rti::geometry* geo)
    {
        geometries_.push_back(geo);
    }

    /// Returns geometry container (const reference)
    const std::vector<rti::geometry*>&
    get_geometries(){
        return geometries_;
    }

    /*
    void dump(){
        for(auto i : rotation_){
            std::cout<<i.first<<" : " << 
            i.second.angle <<" deg, " << 
            i.second.direction << std::endl;
        }    
        for(auto i : position_){
            std::cout<<i.first<<" : " 
                     << i.second.x << ", " 
                     << i.second.y << ", "
                     << i.second.z << std::endl;
        }
    }
    */
};

}
#endif
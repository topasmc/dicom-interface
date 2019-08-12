#ifndef RTI_STOPPING_POWER_H
#define RTI_STOPPING_POWER_H

/// \file
///
/// Stopping power interface

#include <random>
#include <functional>
#include <map>
#include <array>
#include <string>

#include <rti/base/rti_utils.hpp>

namespace rti{

    #include "rti_stopping_power.icc"

    /// Returns stopping power for given material and energy
    inline float 
    stopping_power(
       float energy, 
       const char* material="water")
    {
        if( !std::strcmp(material, "water") ){
            return rti::interp_linear<float,6>(icru_stopping_power_water, energy, 2);
        }else if( !std::strcmp(material, "air") ){
            return rti::interp_linear<float, 6>(icru_stopping_power_air, energy, 2);
        }   
        assert("Can't find given stopping power table");
        return -1.0f;
    }

    /// Returns water equivalent thickness of air of pos_z (thickness)
    inline float 
    wet_air(
       float energy,  //MeV
       float pos_z)   //cm
    {  
       //density: 0.00120479
       
       return pos_z*(0.00120479)* stopping_power(energy, "air")/stopping_power(energy, "water");  
    }

} //namespace

#endif

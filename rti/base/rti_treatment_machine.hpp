#ifndef RTI_TREATMENT_MACHINE_H
#define RTI_TREATMENT_MACHINE_H

/// \file
///
/// Abstraction for treatment machine 

#include <rti/base/rti_beamline.hpp>
#include <rti/base/rti_beamsource.hpp>
#include <rti/base/rti_beam_module.hpp>
#include <rti/base/rti_beamlet.hpp>

namespace rti{

/// \class treatment_machine
/// 
/// Describes an abstraction of all types (RT and ION) treatment machine
/// Typically it consists of geometries (beam limiting devices) and sources (particle fluence)
/// This class has three user-implemented methods, create_beamline, create_beamsource, and create_coordinate_transform
/// \tparam T type of phase-space variables, e.g., float or double
/// \note treatment machine is distinguished in its name
template<typename T>
class treatment_machine {
protected:

    ///< Machine name in string (site:system:mc_code)
    const std::string name_             ;

    ///< Source to Axis distance, 
    ///< neccessary to calculate beam divergence
    std::array<float, 2>  SAD_  ;

    ///< Distance from phase-space plan to isocenter
    float source_to_isocenter_mm_ ; 

public: 

    /// Default constructor
    treatment_machine(){;}

    /*
    /// Default constructor to set infinite SAD
    /// \note not sure yet to include.
    treatment_machine(): 
        SAD_{std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity()}
    {;}
    */
    
    /// Default destructor
    ~treatment_machine(){;}

    /// Returns beamline model 
    /// \param ds : pointer of dataset
    /// \param m  : modality type such as RTIP, US, PASSIVE, etc 
    /// \return beamline
    virtual 
    rti::beamline<T> 
    create_beamline(
        const rti::dataset* ds, 
        rti::modality_type m
    ) = 0 ;

    
    
    /// Returns beamsource model 
    /// \param ds : pointer of dataset
    /// \param m  : modality type such as IMPT, US, PASSIVE, etc
    /// \param pcoord: coordinate system of beam geometry including gantry, patient support, etc
    /// \param scalefactor : total number of histories can be calculated using this number   
    virtual 
    rti::beamsource<T>
    create_beamsource(
        const rti::dataset* ds,
        const rti::modality_type  m,
        const rti::coordinate_transform<T> pcoord,
        const float scalefactor = -1,
        const float source_to_isocenter_mm = 390.0  
    ) = 0 ;

    
    /// Returns coordinate transform information 
    /// \param ds : pointer of dataset
    /// \param m  : modality type such as IMPT, US, PASSIVE, etc
    /// \return pcoord: coordinate system of beam geometry including gantry, patient support, etc
    virtual 
    rti::coordinate_transform<T>
    create_coordinate_transform(
        const rti::dataset* ds,
        const rti::modality_type m
    ) = 0 ;
    

};

}

#endif
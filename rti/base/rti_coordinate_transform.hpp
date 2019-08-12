#ifndef RTI_COORDINATE_TRANSFORM_HPP
#define RTI_COORDINATE_TRANSFORM_HPP

/// \file
///
/// A coordinate transform to map a history from beamsource to map a IEC coordinate or DICOM coordinate.

#include <iostream>
#include <fstream>
#include <random>
#include <functional>
#include <array>
#include <tuple>

#include <rti/base/rti_vec.hpp>
#include <rti/base/rti_matrix.hpp>


namespace rti{

/// \class coordinate_transform
/// default units are mm for length and degree for angle.
/// rotation direction is CCW
/// order of rotation is collimator->gantry->couch->iec2dicom
/// \tparam T for types for return values by the distributions
template <typename T>
class coordinate_transform{
public:

    /// A constant to convert degree to radian.
    const float deg2rad = M_PI/180.0;

    /// Final rotation matrix and translation vector
    rti::mat3x3<T>  rotation    ;
    rti::vec3<T>    translation ;

    /// rotation matrices 
    rti::mat3x3<T>  collimator      ; ///< Rotation due to collimator angle
    rti::mat3x3<T>  gantry          ; ///< Rotation due to gantry angle
    rti::mat3x3<T>  patient_support ; ///< Rotation due to couch
    rti::mat3x3<T>  iec2dicom       ; ///-90 deg (iec2dicom), 90 deg (dicom2iec)

    
    /// Constructor with 4-angles and position
    /// \param angles[4] angles of collimator, gantry, couch, and iec
    /// \param pos move origin, i.e, isocenter
    CUDA_HOST_DEVICE
    coordinate_transform(
        std::array<T, 4>& angles, 
        rti::vec3<T>& pos) 
    : translation(pos)
    {
        collimator = rti::mat3x3<T>(0, 0, angles[0]);
        gantry     = rti::mat3x3<T>(0, angles[1], 0);
        patient_support      = rti::mat3x3<T>(0, 0, angles[2]);
        iec2dicom  = rti::mat3x3<T>(angles[3], 0, 0); 
        rotation = iec2dicom * patient_support *  gantry * collimator    ; 
    }

    /// Constructor with constant 4-angles and position
    /// \param const angles[4] angles of collimator, gantry, couch, and iec
    /// \param const pos
    CUDA_HOST_DEVICE
    coordinate_transform(
        const std::array<T, 4>& angles, 
        const rti::vec3<T>& pos) 
    : translation(pos)
    {
        collimator = rti::mat3x3<T>(0,0, angles[0]);
        gantry     = rti::mat3x3<T>(0, angles[1], 0);
        patient_support      = rti::mat3x3<T>(0, 0, angles[2]);
        iec2dicom  = rti::mat3x3<T>(angles[3], 0, 0);
        rotation = iec2dicom * patient_support *  gantry * collimator    ; 
    }

    /// A copy constructor
    CUDA_HOST_DEVICE
    coordinate_transform(const coordinate_transform<T>& ref)
    {
        collimator = ref.collimator;
        gantry     = ref.gantry;
        patient_support      = ref.patient_support;
        iec2dicom  = ref.iec2dicom; 
        rotation = iec2dicom * patient_support *  gantry * collimator    ; 
        translation  = ref.translation;
    }

    /// Destructor
    CUDA_HOST_DEVICE
    coordinate_transform(){;}

    /// Assignment operator
    CUDA_HOST_DEVICE
    coordinate_transform<T>&
    operator=(const coordinate_transform<T>& ref)
    {
        collimator = ref.collimator;
        gantry     = ref.gantry;
        patient_support      = ref.patient_support;
        iec2dicom  = ref.iec2dicom; 
        rotation = iec2dicom * patient_support *  gantry * collimator    ; 
        translation  = ref.translation;
        return *this;
    }

    /// Prints out matrix
    CUDA_HOST
    void
    dump()
    {
        std::cout<<"--- coordinate transform---"<< std::endl;
        std::cout<<"    translation ---"<< std::endl;
        translation.dump();
        std::cout<<"    rotation ---"<< std::endl;
        rotation.dump();
        std::cout<<"    gantry ---"<< std::endl;
        gantry.dump();
        std::cout<<"    patient_support ---"<< std::endl;
        patient_support.dump();
        std::cout<<"    collimator ---"<< std::endl;
        collimator.dump();
        std::cout<<"    IEC2DICOM ---"<< std::endl;
        iec2dicom.dump();
    }

};

}

#endif 

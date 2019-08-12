#ifndef RTI_BEAM_MODULE_H
#define RTI_BEAM_MODULE_H

/// \file
///
/// Top abstraction to interpret DICOM beam modules, RT and RT-Ion.

#include <rti/base/rti_dataset.hpp>

namespace rti {

/// \class beam_module
///
/// A class that interprets IonControlPoints and converts fluence_map to be actually delivered.
/// beam_module        : provides overall interface to MC engine
/// beam_module_rt     : photon/electron  -> DYNAMIC, STATIC
/// beam_module_rtion  : proton/particles -> UNIFORM, MODULATED, MODULATED_SPEC
class beam_module { 
public:    
protected: 
    const rti::dataset* ds_; ///< Item of IonBeamSequence
    const rti::modality_type modality_ ;
    
public:
    beam_module(
        const rti::dataset* d,
        rti::modality_type m)
        :ds_(d), modality_(m)
    {;}
    ~beam_module(){;}
    virtual void dump() const {;}


};

}

#endif

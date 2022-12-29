//
// ********************************************************************
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * TOPAS collaboration.                                             *
// * Use or redistribution of this code is not permitted without the  *
// * explicit approval of the TOPAS collaboration.                    *
// * Contact: Joseph Perl, perl@slac.stanford.edu                     *
// *                                                                  *
// ********************************************************************
// author: Jungwook Shin @ mgh (jshin13@mgh.harvard.edu)

#ifndef TSRTION_COMPONENT_HH
#define TSRTION_COMPONENT_HH

/// \file
///
/// TOPAS geometry extiontion for RT-Ion component


#include <rti/base/rti_treatment_session.hpp>
#include <rti/base/rti_rtdose.hpp>

#include "TsVGeometryComponent.hh"

/// \class TsRTIonComponents
class TsRTIonComponents : public TsVGeometryComponent
{
protected: 

	std::unique_ptr<rti::treatment_session<float>> rti_session_ {nullptr} ; 

	rti::beamline<float> beamline_;

	///< Get corrdinate system first
    rti::coordinate_transform<float> rt_coordinate_dicom_;
	rti::coordinate_transform<float> rt_coordinate_topas_;

	///< Local origin w.r.t World origin
	const G4ThreeVector* pos_global_ ; 
	
	///< shared by all components
	const G4RotationMatrix* rot_local_ ;  
	const G4RotationMatrix* rot_global_  ; 
	const G4RotationMatrix* rot_posture_ ; 

	G4ThreeVector*
	ConvertRTIPosition2TOPASPosition(
		const rti::vec3<float>& p
	);

    void
    ExportDICOMCoordinateToParameters(
        rti::coordinate_transform<float>& p
    );

    void ExportDoseGridToParameters(
		rti::rect3d<float, float>& dg, 
		rti::vec3<float>& ct_p,
		rti::vec3<float>& trans_p);

	void 
	Wrap3Vector(
		rti::vec3<float> vec, 
		const char* name, 
		const char* unit=" mm"
	);

	void 
	Wrap1Vector(
		float value, 
		const char* name, 
		const char* unit=" deg"
	);

	void 
	Wrap3IntVector(
		rti::vec3<int> vec,
		const char* name
	);
	

	bool 
	ConstructRangeShifter(rti::rangeshifter* rs);

	bool 
	ConstructBlock(rti::aperture* apt);

	

public:
	TsRTIonComponents(
		TsParameterManager* pM, 
		TsExtensionManager* eM, 
		TsMaterialManager*  mM, 
		TsGeometryManager*  gM,
		TsVGeometryComponent* parentComponent, 
		G4VPhysicalVolume* parentVolume, 
		G4String& name);
	~TsRTIonComponents();

	G4VPhysicalVolume* Construct();

	/*
	bool
	ConstructDoseGrid(
		//std::shared_ptr<const rti::dosegrid> dg
	);
	*/


};

#endif

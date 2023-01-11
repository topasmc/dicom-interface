// Component for TsRTIonComponents
// This code implementation is the intellectual property of the MGH.
// by Jungwook Shin jshin13@mgh.harvard.edu

/// \file
///
/// This is an user extension for RT-Ion plan (RTIP) and beam treatment record (RTIBTR)


#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include <iostream>
#include <map>
#include <vector>
#include <array>
#include <tuple>
#include <algorithm>
#include <fstream>

#include "G4LogicalVolume.hh"
#include "G4UIcommand.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"
#include "G4ExtrudedSolid.hh"
#include "G4PVPlacement.hh"

//Headers for TOPAS
#include "TsParameterManager.hh"

//Headers for extension
#include "TsRTIonComponents.hh"

TsRTIonComponents::TsRTIonComponents(
	TsParameterManager *pM,
	TsExtensionManager *eM,
	TsMaterialManager *mM,
	TsGeometryManager *gM,
	TsVGeometryComponent *parentComponent,
	G4VPhysicalVolume *parentVolume,
	G4String &name)
	: TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
	//Important: 'false' caused error
	fIsGroup = true;

	G4String rt_ion_file = fPm->GetStringParameter(GetFullParmName("File"));

	G4String machine_name = "";
	if (fPm->ParameterExists(GetFullParmName("machinename")))
	{
		machine_name = fPm->GetStringParameter(GetFullParmName("machinename"));
	}

	rti_session_ =
		std::unique_ptr<rti::treatment_session<float>>(new rti::treatment_session<float>(rt_ion_file, machine_name));
}

TsRTIonComponents::~TsRTIonComponents() { ; }

G4VPhysicalVolume *
TsRTIonComponents::Construct()
{

	rti::beam_id_type beam_id = {rti::beam_id_type::NUM, 1};
	if (fPm->ParameterExists(GetFullParmName("beamnumber")))
	{
		G4int bnb = fPm->GetIntegerParameter(GetFullParmName("beamnumber"));
		beam_id.number = bnb;
	}
	else
	{
		if (fPm->ParameterExists(GetFullParmName("beamname")))
		{
			G4String bname = fPm->GetStringParameter(GetFullParmName("beamname"));
			beam_id.type = rti::beam_id_type::STR;
			beam_id.name = bname.c_str();
		}
		else
		{
			std::cerr << " Neither BeamNumber nor BeamName was given." << std::endl;
			exit(1);
		}
	}

	/**
    * Export DICOM coordinate information and export as TOPAS parameters
    */
	//rti::coordinate_transform<float> rt_coordinate_dicom;
	switch (beam_id.type)
	{
	case rti::beam_id_type::NUM:
		std::cout<<"Retrieving coordinate information from beam number " << beam_id.number << std::endl;
		rt_coordinate_dicom_ = rti_session_->get_coordinate(beam_id.number);
		break;
	case rti::beam_id_type::STR:
		std::cout<<"Retrieving coordinate information from beam name " << beam_id.name << std::endl;;
		rt_coordinate_dicom_ = rti_session_->get_coordinate(beam_id.name);
		break;
	}

	this->ExportDICOMCoordinateToParameters(rt_coordinate_dicom_);
	std::cout << "RTIonComponent isocenter " << std::endl;
	rt_coordinate_dicom_.translation.dump();

	/**
    * Export Patient CT's center position as transient parameter
	* We need center position of patient and size of patient (for default dosegrid)
    */
	bool is_patient_created = false;
	rti::vec3<float> ct_center(0.0,0.0,0.0);
	rti::vec3<float> ct_size(0,0,0);

	if (fPm->ParameterExists(GetFullParmName("imgdirectory")))
	{

		G4String image_dir = fPm->GetStringParameter(GetFullParmName("imgdirectory"));

        struct stat ss;
        if (stat(image_dir.c_str(), &ss) == 0) //exist
		{
			//check imgdirectory is DIR 
            (ss.st_mode & S_IFDIR) ? true : throw std::runtime_error("Error: ImgDirectory should be a directory");
			
			//check # of files in the directory
			DIR* dir = opendir(image_dir.c_str());
			size_t files = 0;
			struct dirent* ent ;
            while((ent = readdir(dir))) files++;
			closedir(dir);
			if(files>2){ // . and .. are counted.
				rti::ct<float> patient_ct(image_dir);
				ct_center = patient_ct.get_center();
				ct_size = patient_ct.get_size();
				is_patient_created = true;
			}

		}

		this->Wrap3Vector(ct_center, "ImgCenter", "mm");
		
	}

	//2. TOPAS operation for group
	//BeginConstruction must get called
	//after ExportDICOMCoordinateToParameters()
	//after Export ImgCenter
	BeginConstruction();
	fEnvelopePhys = fParentVolume;

	// Coordinate system handling
	//Angle order: 1,2,3,4
	std::array<float, 4> rotations;
	rotations[0] = static_cast<float>(fPm->GetDoubleParameter(GetFullParmName("RotCollimator"), "Angle") / rad);
	rotations[1] = static_cast<float>(fPm->GetDoubleParameter(GetFullParmName("RotGantry"), "Angle") / rad);
	rotations[2] = static_cast<float>(fPm->GetDoubleParameter(GetFullParmName("RotPatientSupport"), "Angle") / rad);
	rotations[3] = static_cast<float>(fPm->GetDoubleParameter(GetFullParmName("RotIEC2DICOM"), "Angle") / rad);

	std::cout << "RTION geometry coordinate position (x,y,z): " << std::endl;
	std::cout << "RTION geometry Rotation (colli, gantry, patient support, iec2dicom): "
			  << rotations[0] / deg << ", "
			  << rotations[1] / deg << ", "
			  << rotations[2] / deg << ", "
			  << rotations[3] / deg << std::endl;

	rti::vec3<float> p_xyz;
	rt_coordinate_topas_ = rti::coordinate_transform<float>(rotations, p_xyz);

	//Needs three column vector
	rot_local_ = new G4RotationMatrix(
		G4ThreeVector(rt_coordinate_topas_.rotation.xx, rt_coordinate_topas_.rotation.yx, rt_coordinate_topas_.rotation.zx),
		G4ThreeVector(rt_coordinate_topas_.rotation.xy, rt_coordinate_topas_.rotation.yy, rt_coordinate_topas_.rotation.zy),
		G4ThreeVector(rt_coordinate_topas_.rotation.xz, rt_coordinate_topas_.rotation.yz, rt_coordinate_topas_.rotation.zz));

	pos_global_ = GetTransForPlacement();
	rot_global_ = GetRotForPlacement();
	rot_posture_ = new G4RotationMatrix((*rot_global_ * *rot_local_).inverse());

	//
	switch (beam_id.type)
	{
	case rti::beam_id_type::NUM:
		beamline_ = rti_session_->get_beamline(
			beam_id.number);

		break;
	case rti::beam_id_type::STR:
		beamline_ = rti_session_->get_beamline(
			beam_id.name);
		break;
	}

	//5. TOPAS routine to create machine's geometry in TOPAS context
	for (auto geometry : beamline_.get_geometries())
	{
		switch (geometry->geotype)
		{
		case rti::geometry_type::BLOCK:
		{
			std::cout << "Creating Aperture" << std::endl;
			this->ConstructBlock(dynamic_cast<rti::aperture *>(geometry));
			break;
		}
		case rti::geometry_type::RANGESHIFTER:
		{
			std::cout << "Creating Rangeshifter" << std::endl;
			this->ConstructRangeShifter(dynamic_cast<rti::rangeshifter *>(geometry));
			break;
		}
		default:
		{
			std::cout << "Method for this geometry is not implemented." << std::endl;
			break;
		}

		} //switch
	}

	//6. Construct RT-Dose as a part of RTIonComponents
	//   Altneratively in higher TOPAS version ( > 3.0.1, as this is commissioned at MGH)
	//	 We can utilize CloneRTDose as a subcomponent of Patient
	/**
		 * DoseGrid handling
		 * -> DoseGrid is created only when a patient geometry is created
		 * Output dose-grid parameters.
		 */
	if (fPm->ParameterExists(GetFullParmName("includedosegridifexist")) && is_patient_created)
	{

		if (fPm->GetBooleanParameter(GetFullParmName("includedosegridifexist")))
		{

			if (fPm->ParameterExists(GetFullParmName("rtdosefile")))
			{
				G4String dose_file = fPm->GetStringParameter(GetFullParmName("rtdosefile"));
                struct stat ss;
				bool is_default_grid = true;
                if (stat(dose_file.c_str(), &ss) == 0)
				{
                    if (ss.st_mode & S_IFLNK || ss.st_mode & S_IFREG)
					{
						std::cout << "Creating DoseGrid from RTDOSE" << std::endl;
						rti::rtdose<float> dose_grid(dose_file);
						rti::vec3<float> trans_p(pos_global_->x(), pos_global_->y(), pos_global_->z());
						this->ExportDoseGridToParameters(dose_grid, ct_center, trans_p);
						is_default_grid = false;
					}
				}
				if (is_default_grid)
				{
					std::cout << "Creating DoseGrid from default" << std::endl;
					//create default dosegrid
					float nx = std::ceil(ct_size.x / 2.0);
					float ny = std::ceil(ct_size.y / 2.0);
					float nz = std::ceil(ct_size.z / 2.0);
					std::vector<float> x(nx);
					std::vector<float> y(ny);
					std::vector<float> z(nz);
					for (size_t n = 0; n < nx; ++n)
						x[n] = ct_center.x - (nx - 1) + n * 2.0;
					for (size_t n = 0; n < ny; ++n)
						y[n] = ct_center.y - (ny - 1) + n * 2.0;
					for (size_t n = 0; n < nz; ++n)
						z[n] = ct_center.z - (nz - 1) + n * 2.0;
					rti::rect3d<float, float> dose_grid(x, y, z);
					std::cout << "ct---" << std::endl;
					ct_center.dump();
					ct_size.dump();
					std::cout << "dg---" << std::endl;
					dose_grid.get_center().dump();
					dose_grid.get_size().dump();
					rti::vec3<float> trans_p(pos_global_->x(), pos_global_->y(), pos_global_->z());
					this->ExportDoseGridToParameters(dose_grid, ct_center, trans_p);
				} //default grid
			}	 //rtdosefile

#if (TOPAS_VERSION_MAJOR >= 3  && TOPAS_VERSION_MINOR >= 2) || TOPAS_VERSION_MAJOR >= 4
    InstantiateChild(fName + "/DoseGrid")->Construct();
#else
	InstantiateChild(fName + "/DoseGrid", fParentVolume)->Construct();
#endif
		} //include_dosegrid

	} //is_patient_created

	return fEnvelopePhys;
}

void TsRTIonComponents::ExportDICOMCoordinateToParameters(
	rti::coordinate_transform<float> &p)
{
	//1. isocenter
	this->Wrap3Vector(p.translation, "IsoCenter");

	//Angle 0. beam limiting devide pitch angle
	this->Wrap1Vector(p.collimator.z, "CollimatorAngle");

	//Angle 1. gantry angle
	this->Wrap1Vector(p.gantry.y, "GantryAngle");

	//Angle 2. patient support angle
	this->Wrap1Vector(p.patient_support.z, "PatientSupportAngle");

}

void TsRTIonComponents::Wrap3Vector(
	rti::vec3<float> vec,
	const char *name,
	const char *unit)
{
	G4String param_name;
	G4String trans_value;

	param_name = "dc:Ge/" + fName + "/" + name + "X";
	trans_value = G4UIcommand::ConvertToString(vec.x) + " " + unit;
	fPm->AddParameter(param_name, trans_value);

	param_name = "dc:Ge/" + fName + "/" + name + "Y";
	trans_value = G4UIcommand::ConvertToString(vec.y) + " " + unit;
	fPm->AddParameter(param_name, trans_value);

	param_name = "dc:Ge/" + fName + "/" + name + "Z";
	trans_value = G4UIcommand::ConvertToString(vec.z) + " " + unit;
	fPm->AddParameter(param_name, trans_value);
}

void TsRTIonComponents::Wrap3IntVector(
	rti::vec3<int> vec,
	const char *name)
{
	G4String param_name;
	G4String trans_value;

	param_name = "ic:Ge/" + fName + "/" + name + "X";
	trans_value = G4UIcommand::ConvertToString(vec.x);
	fPm->AddParameter(param_name, trans_value);

	param_name = "ic:Ge/" + fName + "/" + name + "Y";
	trans_value = G4UIcommand::ConvertToString(vec.y);
	fPm->AddParameter(param_name, trans_value);

	param_name = "ic:Ge/" + fName + "/" + name + "Z";
	trans_value = G4UIcommand::ConvertToString(vec.z);
	fPm->AddParameter(param_name, trans_value);
}

void TsRTIonComponents::Wrap1Vector(
	float value,
	const char *name,
	const char *unit)
{
	G4String param_name;
	G4String trans_value;

	param_name = "dc:Ge/" + fName + "/" + name;
	trans_value = G4UIcommand::ConvertToString(value) + " " + unit;
	fPm->AddParameter(param_name, trans_value);
}

bool TsRTIonComponents::ConstructRangeShifter(
	rti::rangeshifter *rs)
{
	G4bool include_rs = true;
	if (fPm->ParameterExists(GetFullParmName("includerangeshifterifexist")))
	{
		include_rs = fPm->GetBooleanParameter(GetFullParmName("includerangeshifterifexist"));
	}

	if (!include_rs)
		return false;

	G4String rs_name = "rangeshifter";
	G4String rs_mat = fPm->GetStringParameter(GetFullParmName(G4String(rs_name + "/material")));



	G4VSolid *rs_solid = nullptr;
	const rti::vec3<float> &vol = rs->volume;

	if (rs->is_rectangle){
		rs_solid = new G4Box(rs_name, 0.5 * vol.x * mm, 0.5 * vol.y * mm, 0.5 * vol.z * mm);
	}
	else{
		rs_solid = new G4Tubs(rs_name, 0.0, vol.x * mm, 0.5*vol.z * mm, 0.0 * deg, 360.0 * deg);
	}
	
	//G4Box *rs_solid = new G4Box(rs_name, 0.5 * vol.x * mm, 0.5 * vol.y * mm, 0.5 * vol.z * mm);

	G4LogicalVolume *rs_log = CreateLogicalVolume(rs_name, rs_mat, rs_solid);

	CreatePhysicalVolume(rs_name,
						 rs_log,
						 (G4RotationMatrix *)rot_posture_,
						 this->ConvertRTIPosition2TOPASPosition(rs->pos),
						 fEnvelopePhys);

	return true;
}

G4ThreeVector *
TsRTIonComponents::ConvertRTIPosition2TOPASPosition(
	const rti::vec3<float> &p)
{
	//this method utilizes members of pos_global_, rot_global_
	rti::vec3<float> pos = rt_coordinate_topas_.rotation * p;

	G4ThreeVector pos_local(pos.x, pos.y, pos.z);

	G4ThreeVector *topas_position = new G4ThreeVector(*pos_global_);
	(*topas_position) += (*rot_global_ * pos_local);

	return topas_position;
}

bool TsRTIonComponents::ConstructBlock(
	rti::aperture *apt)
{
	G4bool include_blk = true;
	if (fPm->ParameterExists(GetFullParmName("includeblockifexist")))
	{
		include_blk = fPm->GetBooleanParameter(GetFullParmName("includeblockifexist"));
	}

	if (!include_blk)
		return false;

	//block disk
	G4String blk_name = "block";
	G4String blk_mat = fPm->GetStringParameter(GetFullParmName(G4String(blk_name + "/material")));

	G4VSolid *gblock = nullptr;

	const rti::vec3<float> &vol = apt->volume;

	if (apt->is_rectangle){
		gblock = new G4Box(blk_name, 0.5 * vol.x * mm, 0.5 * vol.y * mm, 0.5 * vol.z * mm);
	}
	else{
		gblock = new G4Tubs(blk_name, 0.0, vol.x * mm, 0.5*vol.z * mm, 0.0 * deg, 360.0 * deg);
	}

	G4LogicalVolume *lblock = CreateLogicalVolume(blk_name, blk_mat, gblock);

	G4VPhysicalVolume *pblock = CreatePhysicalVolume(blk_name,
													 lblock,
													 (G4RotationMatrix *)rot_posture_,
													 this->ConvertRTIPosition2TOPASPosition(apt->pos),
													 fEnvelopePhys);

	//holes
    //int start = 0;

    for (size_t i = 0; i < apt->block_data.size(); ++i)
	{
		auto xypts = apt->block_data[i];

		std::cout << "Hole:" << i << " nb_pts:" << xypts.size() << std::endl;

		std::vector<G4TwoVector> xy_pts;

		for (const auto xy : xypts)
		{
			xy_pts.push_back(G4TwoVector(xy[0], xy[1]));
		}

		G4String vname = blk_name + std::string("_void") + std::to_string(i);
		G4String voidMaterial = fParentComponent->GetResolvedMaterialName();

		G4VSolid *voidSolid = new G4ExtrudedSolid(vname, xy_pts, 0.5 * vol.z * mm, G4TwoVector(0, 0), 1.0, G4TwoVector(0, 0), 1.0);
		G4LogicalVolume *voidLog = CreateLogicalVolume(vname, voidMaterial, voidSolid);
		G4RotationMatrix *rot = new G4RotationMatrix(0.0, 0.0, 0.0); //rotation
		G4ThreeVector *position = new G4ThreeVector(0.0, 0.0, 0.0);
		CreatePhysicalVolume(vname, voidLog, rot, position, pblock);
	}

	return true;
}

void TsRTIonComponents::ExportDoseGridToParameters(
	rti::rect3d<float, float> &dg,
	rti::vec3<float> &ct_p,
	rti::vec3<float> &trans_p)
{

	G4String gridName = fName + "/DoseGrid";
	G4String parameterName;
	G4String transValue;

	//Size
	auto lxyz = dg.get_size();
	parameterName = "d:Ge/" + gridName + "/HLX";
	transValue = G4UIcommand::ConvertToString(0.5 * lxyz.x) + " mm";
	fPm->AddParameter(parameterName, transValue);

	parameterName = "d:Ge/" + gridName + "/HLY";
	transValue = G4UIcommand::ConvertToString(0.5 * lxyz.y) + " mm";
	fPm->AddParameter(parameterName, transValue);

	parameterName = "d:Ge/" + gridName + "/HLZ";
	transValue = G4UIcommand::ConvertToString(0.5 * lxyz.z) + " mm";
	fPm->AddParameter(parameterName, transValue);

	//Position & rotation
	/*
	 * As DoseGrid should not move along with RTIonComponents as it is same coordinate system with CT
	 * but the vector for compensation for the  
	 */
	auto pxyz = dg.get_center() - ct_p - trans_p;
	parameterName = "d:Ge/" + gridName + "/TransX";
	transValue = G4UIcommand::ConvertToString(pxyz.x) + " mm ";
	fPm->AddParameter(parameterName, transValue);

	parameterName = "d:Ge/" + gridName + "/TransY";
	transValue = G4UIcommand::ConvertToString(pxyz.y) + " mm ";

	fPm->AddParameter(parameterName, transValue);

	parameterName = "d:Ge/" + gridName + "/TransZ";
	transValue = G4UIcommand::ConvertToString(pxyz.z) + " mm ";
	fPm->AddParameter(parameterName, transValue);

	parameterName = "d:Ge/" + gridName + "/RotX";
	transValue = "0. deg";
	fPm->AddParameter(parameterName, transValue);

	parameterName = "d:Ge/" + gridName + "/RotY";
	transValue = "0. deg";
	fPm->AddParameter(parameterName, transValue);

	parameterName = "d:Ge/" + gridName + "/RotZ";
	transValue = "0. deg";
	fPm->AddParameter(parameterName, transValue);

	//Bins
	rti::vec3<size_t> nxyz = dg.get_nxyz();
	parameterName = "i:Ge/" + gridName + "/XBins";
	transValue = G4UIcommand::ConvertToString(G4int(nxyz.x));
	fPm->AddParameter(parameterName, transValue);

	parameterName = "i:Ge/" + gridName + "/YBins";
	transValue = G4UIcommand::ConvertToString(G4int(nxyz.y));
	fPm->AddParameter(parameterName, transValue);

	parameterName = "i:Ge/" + gridName + "/ZBins";
	transValue = G4UIcommand::ConvertToString(G4int(nxyz.z));
	fPm->AddParameter(parameterName, transValue);

	parameterName = "s:Ge/" + gridName + "/Type";
	transValue = "\"TsBox\"";
	fPm->AddParameter(parameterName, transValue);

	parameterName = "b:Ge/" + gridName + "/IsParallel";
	transValue = "\"True\"";
	fPm->AddParameter(parameterName, transValue);
}

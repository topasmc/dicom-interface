// Particle Source for TsRTIonSource
// This code implementation is the intellectual property of the MGH.
// this is an user extension for RT-Ion plan (RTIP) and beam treatment record (RTIBTR)
// by Jungwook Shin jshin13@mgh.harvard.edu

//Headers for C library
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

//Headers for C++ std library
#include <iostream>
#include <map>
#include <vector>
#include <array>
#include <tuple>
#include <algorithm>

//Headers for Geant4
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4RotationMatrix.hh"
#include "G4UIcommand.hh"
#include "G4SystemOfUnits.hh"
#include "G4ExtrudedSolid.hh"
#include "G4PVPlacement.hh"

//Headers for TOPAS
#include "TsParameterManager.hh"

//Headers for extension
#include "TsRTIonSource.hh"

#ifdef TOPAS_MT
#include "G4AutoLock.hh"
namespace {G4Mutex move_to_nextspot = G4MUTEX_INITIALIZER;}
#endif

TsRTIonSource::TsRTIonSource(TsParameterManager* pM, 
                             TsSourceManager* sM, 
                             G4String name)
                             :TsSource(pM, sM, name)
{

    std::cout<<"Creating RTIon Source."<<std::endl;

    G4String rt_ion_file   = fPm->GetStringParameter(GetFullParmName("File"));
	
    G4String machine_name = "";
    
    if ( fPm->ParameterExists(GetFullParmName("machinename")) ){
        machine_name = fPm->GetStringParameter(GetFullParmName("machinename")) ;	
    }
    
    rti_session_ =  
        std::unique_ptr<rti::treatment_session<float>>
        (new rti::treatment_session<float>(rt_ion_file, machine_name)); //topas version

    
    //ReadHistories will be accessed from RTIonSourceGenerator
    ReadHistories = std::bind(
                              &TsRTIonSource::ReadHistoriesPerBeamlet, 
                              this, 
                              std::placeholders::_1);

    this->ResolveParameters();
}

TsRTIonSource::~TsRTIonSource() {;}


void TsRTIonSource::ResolveParameters(){
    
    /**
     * Particle type: proton (default)
     */
    if(fPm->ParameterExists(GetFullParmName("beamparticle"))){
        TsParticleDefinition resolvedDef = 
            fPm->GetParticleDefinition(fPm->GetStringParameter(GetFullParmName("beamparticle")));
        particle_definition_ = resolvedDef.particleDefinition;
        if (!particle_definition_)
            std::cerr << "Unknown particle type read from parameter = " << GetFullParmName("beamparticle") << std::endl;
    }else{
        particle_definition_ = G4Proton::ProtonDefinition();
    }
    
    /**
     * Beam id
     * \brief beam number (required) or beam name (optional)
     * If beam number is given, we ignore beam name.
     * When, beam name is used primarily in your clinic, it is important not to use beam number.
     * Because in paramter file chain, beam number might be defined elsewhere.
     */
    rti::beam_id_type beam_id = {rti::beam_id_type::NUM, 1};
    if ( fPm->ParameterExists(GetFullParmName("beamnumber")) ){
        G4int bnb = fPm->GetIntegerParameter(GetFullParmName("beamnumber"));
        beam_id.number = bnb;
    }else{
        if ( fPm->ParameterExists(GetFullParmName("beamname")) ){
            G4String bname = fPm->GetStringParameter(GetFullParmName("beamname")) ;
            beam_id.type = rti::beam_id_type::STR;
            beam_id.name = bname.c_str();
        }else{
            std::cerr << " Neither BeamNumber nor BeamName was given." << std::endl;
            exit(1);
        }
    }

    
    /**
     * Export DICOM coordinate information and export as TOPAS parameters
     */
    rti::coordinate_transform<float> rt_coordinate_dicom;
    switch(beam_id.type){
    case rti::beam_id_type::NUM: 
        rt_coordinate_dicom = rti_session_->get_coordinate(beam_id.number);
        break;
    case rti::beam_id_type::STR: 
        rt_coordinate_dicom = rti_session_->get_coordinate(beam_id.name);
        break;
    }
    this->ExportDICOMCoordinateToParameters(rt_coordinate_dicom);
    std::cout<<"isocenter: " << std::endl;
    rt_coordinate_dicom.translation.dump();

    /**
     * SID: source to isocenter distance
     */
    float sid_mm = 0.0*mm ;
    if(fPm->ParameterExists(GetFullParmName("sid"))){
        sid_mm = fPm->GetDoubleParameter(GetFullParmName("sid"), "Length")/mm;
    }

    /**
     * patient's position information & dosegrid
     */
    G4String transValue    ;
    G4String parameterName ;

    if ( fPm->ParameterExists(GetFullParmName("imgdirectory"))){
        //Writing img center should come before readign TransX/Y/Z
        G4String image_dir = fPm->GetStringParameter(GetFullParmName("imgdirectory"));

        struct stat s;
        if (stat(image_dir.c_str(), &s) == 0) //exist
            {
                //check imgdirectory is DIR 
                (s.st_mode & S_IFDIR) ? true : throw std::runtime_error("Error: ImgDirectory should be a directory");
			
                //check # of files in the directory
                DIR* dir = opendir(image_dir.c_str());
                size_t files = 0;
                struct dirent* ent ;
                while(ent = readdir(dir)) files++;
                closedir(dir);
                if(files>2){ // . and .. are counted.
                    rti::ct<float> patient_ct(image_dir);
                    auto ct_center = patient_ct.get_center();
                    this->Wrap3Vector(ct_center, "ImgCenter", "mm");
                }
            }
    }
    
    /**
     * Updated coordinate from TOPAS UI
     * Rotation of Gantry, Coli, Patient support, and iso-center shift
     */
    rti::vec3<float> p_xyz;
    p_xyz.x = static_cast<float> (fPm->GetDoubleParameter(GetFullParmName("ShiftX"), "Length")/mm);
    p_xyz.y = static_cast<float> (fPm->GetDoubleParameter(GetFullParmName("ShiftY"), "Length")/mm);
    p_xyz.z = static_cast<float> (fPm->GetDoubleParameter(GetFullParmName("ShiftZ"), "Length")/mm);

    //Geant4's rotation direction is same with IEC
    //i.e.: +angle -> CCW rotation viewing from the plus of rotating axis
    //TOPAS: CW viewing from +axis
    //RTI (DICOM): CW viewing from -axis

    //Angle order: 1,2,3,4
    std::array<float, 4> rotations;
    rotations[0] = static_cast<float> (fPm->GetDoubleParameter(GetFullParmName("RotCollimator"), "Angle")/rad);
    rotations[1] = static_cast<float> (fPm->GetDoubleParameter(GetFullParmName("RotGantry"), "Angle")/rad);
    rotations[2] = static_cast<float> (fPm->GetDoubleParameter(GetFullParmName("RotPatientSupport"), "Angle")/rad);
    rotations[3] = static_cast<float> (fPm->GetDoubleParameter(GetFullParmName("RotIEC2DICOM"), "Angle")/rad);
    
    std::cout<<"RTION source coordinate position (x,y,z): " << std::endl; 
    std::cout<<"RTION source Rotation (colli, gantry, patient support, iec2dicom): " 
             << rotations[0]/deg <<", "
             << rotations[1]/deg <<", "
             << rotations[2]/deg <<", "
             << rotations[3]/deg << std::endl;

    rti::coordinate_transform<float> rt_coordinate_topas(rotations, p_xyz);
    rt_coordinate_topas.rotation.dump();

    /**
     * Particles per history
     * 
     */    
    G4double particles_per_history = fPm->GetUnitlessParameter(GetFullParmName("particlesperhistory"));
    assert( (particles_per_history==0) || (particles_per_history==-1) || (particles_per_history>0) );

    switch(beam_id.type){
        case rti::beam_id_type::NUM: 
            beam_source_ = rti_session_->get_beamsource(
                                                    beam_id.number,
                                                    rt_coordinate_topas,
                                                    particles_per_history,
                                                    sid_mm
                                                    );

            break;
        case rti::beam_id_type::STR: 
            beam_source_ = rti_session_->get_beamsource(
                                                    beam_id.name,
                                                    rt_coordinate_topas,
                                                    particles_per_history,
                                                    sid_mm
                                                    );
            break;
    }

    /**
     * Set Beamlet range to be simulated. 
     * By default, starts with first beamlet (id=1) and ends with last beamlets
     */    
    if(fPm->ParameterExists(GetFullParmName("beamletrange"))){ 

        if( fPm->GetVectorLength(GetFullParmName("beamletrange")) != 2){
            throw std::runtime_error("Beamlet start_id is lower than end_id.");
        }

        int* v = fPm->GetIntegerVector(GetFullParmName("beamletrange"));
        
        (v[0] <= 0)    ? throw std::runtime_error("Beamlet range error. start id should be from 1") : false ;
        (v[0] >= v[1]) ? throw std::runtime_error("Beamlet range error. [start=1, stop=end]") : false ;
        
        //due to index start 0 in beam source
        counter_.start = v[0] - 1 ; 
        counter_.stop  = v[1] - 1 ;

    }else{
        counter_.start   = 0;
        counter_.stop    = beam_source_.total_beamlets()-1;
    }

    counter_.current = counter_.start;
    fNumberOfHistoriesInRun = std::get<2>(beam_source_[counter_.stop]) 
        - std::get<2>(beam_source_[counter_.start])
        + std::get<1>(beam_source_[counter_.start]);
    std::cout<<"Total number of histories to be simulated from RTIonSource: " << fNumberOfHistoriesInRun 
             <<", Total number of beamlets: "<< beam_source_.total_beamlets() << std::endl;
    
}

G4bool TsRTIonSource::ReadHistoriesPerBeamlet(std::queue<TsPrimaryParticle>* particleBuffer){
    #ifdef TOPAS_MT
	    G4AutoLock l(&move_to_nextspot);
    #endif
    
    if( counter_.current > counter_.stop) return false;
    
    auto current_beamlet        = beam_source_[counter_.current++];
    auto beamlet_distribution   = std::get<0>(current_beamlet);
    size_t histories_of_beamlet = std::get<1>(current_beamlet);

    while(histories_of_beamlet--){
        
        auto h   = beamlet_distribution();
        TsPrimaryParticle p;
        p.particleDefinition = particle_definition_; 

        p.posX = h.pos.x; 
        p.posY = h.pos.y;
        p.posZ = h.pos.z;
        
        p.dCos1 = h.dir.x ;
        p.dCos2 = h.dir.y ;
        p.dCos3 = h.dir.z ;
        p.kEnergy         = h.ke;
        p.weight          = 1;
        p.isNewHistory    = true;
        p.isOpticalPhoton = false;
        p.isGenericIon    = false;
        p.ionCharge       = 1 ;

        particleBuffer->push(p);
    }
    
    return true;
}


bool
TsRTIonSource::ExportDICOMCoordinateToParameters(
        rti::coordinate_transform<float>& p
){
    //1. isocenter
	this->Wrap3Vector(p.translation, "IsoCenter");

    //Angle 0. beam limiting devide pitch angle
	this->Wrap1Vector(p.collimator.z, "CollimatorAngle");

	//Angle 1. gantry angle
	this->Wrap1Vector(p.gantry.y , "GantryAngle");

	//Angle 2. patient support angle
	this->Wrap1Vector(p.patient_support.z, "PatientSupportAngle");

    //Angle 3. IEC2DICOM angle
    this->Wrap1Vector(p.iec2dicom.x, "Iec2DicomAngle");

}

void 
TsRTIonSource::Wrap3Vector(
	rti::vec3<float> vec, 
    const char* name, 
    const char* unit)
{
	G4String param_name; 
	G4String trans_value;

	param_name = "dc:So/" + fSourceName + "/" + name + "X";
	trans_value = G4UIcommand::ConvertToString(vec.x) + " " + unit;
	fPm->AddParameter(param_name, trans_value);

	param_name = "dc:So/" + fSourceName + "/" + name + "Y";
	trans_value = G4UIcommand::ConvertToString(vec.y)+ " " + unit;
	fPm->AddParameter(param_name, trans_value);

	param_name = "dc:So/" + fSourceName + "/" + name + "Z";
	trans_value = G4UIcommand::ConvertToString(vec.z) + " " + unit;
	fPm->AddParameter(param_name, trans_value);
}

void 
TsRTIonSource::Wrap1Vector(
    float value, 
    const char* name, 
    const char* unit)
{
	G4String param_name; 
	G4String trans_value;

	param_name = "dc:So/" + fSourceName + "/" + name ;
	trans_value = G4UIcommand::ConvertToString(value) + " " + unit;
	fPm->AddParameter(param_name, trans_value);
}

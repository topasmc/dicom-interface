#include <queue>
#include <thread>
#include <functional>

#include <rti/base/rti_treatment_session.hpp>
#include "cli_rti.hpp"


using namespace std;

//typedef double phase_space_type;
typedef float phase_space_type;


void describe(const rti::dataset* beam);


int main(int argc, char** argv){

    const std::string file_name(argv[1]);
    const int         beam_idx = std::stoi(argv[2]) ;
    
    rti::modality_type mtype_;
    rti::dataset* rti_ds_ = nullptr;
    const rti::dataset* beam_ds ;
    
    ///< Sequence tag dictionary for modality specific
    const std::map<
        const std::string,
        const gdcm::Tag
        >* seq_tags_ ;

    
    gdcm::Reader reader;

    reader.SetFileName(file_name.c_str());
        
    const bool is_file_valid = reader.Read();
    if(!is_file_valid) throw std::runtime_error("Invalid DICOM file is given to treatment_session.");

    gdcm::MediaStorage ms;
    ms.SetFromFile(reader.GetFile());

    switch(ms){
    case gdcm::MediaStorage::RTPlanStorage :
        mtype_ = rti::RTPLAN;
        break;
        //case gdcm::MediaStorage::RTBeamsTreatmentRecordStorage :
        //gdcm doesn't have definition of RTBeamsTreatmentRecordStorge
        //case gdcm::MediaStorage::RTBeamsTreatmentRecordStorage : //1.2.840.10008.5.1.4.1.1.481.4
        //    mtype_ = RTRECORD;
        //    break;
    case gdcm::MediaStorage::RTIonPlanStorage :
        mtype_ = rti::IONPLAN;
        break;
    case gdcm::MediaStorage::RTIonBeamsTreatmentRecordStorage :
        mtype_ = rti::IONRECORD;
        break;
    default:
        throw std::runtime_error("treatment_session does not supports given RTMODALITY");
    }

    rti_ds_   = new rti::dataset(reader.GetFile().GetDataSet(), true);
    seq_tags_ = &rti::seqtags_per_modality.at(mtype_);

    beam_ds = (*rti_ds_)( seq_tags_->at("beam"))[beam_idx];

    describe(beam_ds);

    //rti::beam_module_ion ion_beam(beam_ds, mtype_);
    std::string tune_id_;
    std::vector<int>         nb_pts(1)  ;
    std::vector<float>       energy(1)  ;
    std::vector<float>       fwhm_xy(2) ;
    std::vector<std::string> tune_id(1) ;
    std::vector<float> xy     ; 
    std::vector<float> weight ; 

    
    int layer_nb = 0;
    auto seq_tags = &rti::seqtags_per_modality.at(mtype_);
    auto ictrl    = (*beam_ds)(seq_tags->at("ctrl"));
    for(auto b : ictrl){
            
        if ( (layer_nb++)%2) continue;
        //b->get_values("ScanSpotTuneID", tune_id);
        b->get_values("NominalBeamEnergy", energy);
        std::cout <<"Energy: " << energy[0] << "\n";
        //b->get_values("NumberOfScanSpotPositions", nb_pts) ;
        //b->get_values("ScanningSpotSize", fwhm_xy);
        //b->get_values("ScanSpotPositionMap", xy);

        //b->get_values(str_weight.c_str(), weight); //ScanSpotmetersetweights?
        //PrivateElementWithEmptyPrivateCreator: (300b,1094), VR: ??, length: 160
        // tag (300b,1094)
        // vl = 160

        //Used for Sumitomo Line Scanning.
        //>Line Spot Tune ID (300B,1090) SH 3 
        //>Number of Line Scan Spot  Positions (300B,1092) IS 3 
        //>Line Scan Position Map (300B,1094) FL 3 
        //>Line Scan Meterset Weights (300B,1096) FL 3 
        //>Line Scanning Spot Size (300B,1098) FL 3 
        //>Number of Line Scan Spot Paintings (300B,109A) IS 3
        b->dump();

        b->get_values( (*b)[gdcm::Tag(0x300b,0x1090)], tune_id);
        b->get_values( (*b)[gdcm::Tag(0x300b,0x1092)], nb_pts);
        b->get_values(gdcm::Tag(0x300b,0x1094), xy);
        b->get_values(gdcm::Tag(0x300b,0x1096), weight);
        b->get_values(gdcm::Tag(0x300b,0x1098), fwhm_xy);

        //std::cout << "Tune ID size: "<< tune_id.size() << ", value: " << tune_id[0] << "\n";
        
        //(*b)[gdcm::Tag(0x300b,0x1094)]; //-> gdcm::DataElement ?
        //b->print_dataelement( (*b)[gdcm::Tag(0x300b,0x1094)] );
        //std::cout<< "\n";

        for(int j = 0 ; j < nb_pts[0] ; ++j){
            //tune_id_ = tune_id[0];
            std::cout << energy[0] <<", "
                      << xy[j*2] << ", " << xy[j*2+1] << ", "
                      << fwhm_xy[0] << ", " << fwhm_xy[1] << ", "
                      << weight[j] << "\n";
            //sequence_.push_back({energy[0], xy[j*2], xy[j*2+1], fwhm_xy[0], fwhm_xy[1], weight[j]});
        }//per spot


        //nb_spots_per_layer_.push_back(nb_pts[0]);

    }//per layer


    
    return 0;
}

void describe(const rti::dataset* beam)
{
    std::cout << "======== BEAM info =======" << std::endl;
    std::vector<int> data_int;
    std::vector<int> data_int_tag;
    std::vector<float> data_float;
    std::vector<std::string> data_str;
    
    beam->get_values("BeamNumber", data_int);
    std::cout<< "BeamNumber: " << data_int[0] << std::endl;

    beam->get_values("BeamName", data_str);
    std::cout<< "BeamName: " << data_str[0] << std::endl;

    beam->get_values("BeamType", data_str);
    std::cout<< "BeamType: " << data_str[0] << std::endl;

    beam->get_values("TreatmentDeliveryType", data_str);
    std::cout<<"TreatmentDeliveryType : " << data_str[0] << std::endl;

    beam->get_values("NumberOfWedges", data_int);
    std::cout <<"NumberOfWedges: " << data_int[0] << std::endl;

    beam->get_values("NumberOfRangeShifters", data_int);
    std::cout <<"NumberOfRangeShifters: " << data_int[0] << std::endl;

    beam->get_values("NumberOfCompensators", data_int);
    std::cout <<"NumberOfCompensators : " << data_int[0] << std::endl;

    beam->get_values("NumberOfBoli", data_int);
    std::cout <<"Number Of Boli : " << data_int[0] << std::endl;

    beam->get_values("NumberOfBlocks", data_int);
    std::cout <<"NumberOfBlocks : " << data_int[0] << std::endl;

    beam->get_values("NumberOfControlPoints", data_int);
    std::cout <<"NumberOfControlPoints : " << data_int[0] << std::endl;
    
    beam->get_values("VirtualSourceAxisDistances", data_float);
    std::cout <<"VirtualSourceAxisDistances : " << data_float[0] <<", " << data_float[1] << std::endl;

}

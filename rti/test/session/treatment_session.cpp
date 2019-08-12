
#include <queue>
#include <thread>
#include <functional>

#include <rti/base/rti_treatment_session.hpp>
#include "cli_rti.hpp"


using namespace std;

//typedef double phase_space_type;
typedef float phase_space_type;


void describe(const rti::dataset* beam);
void describe_rtibtr(const rti::dataset* beam);

template<class T>
void generator(rti::beamlet<T>& bl, size_t thread_id, size_t n_histories);


int main(int argc, char** argv){
    
    rti::cli_beam_read cl_opts;
    cl_opts.read(argc, argv);

    std::string machine_name = cl_opts["--machine"].size() == 0 ? "" : cl_opts["--machine"][0];
    
    std::string mc_code = cl_opts["--mc_code"].size() == 0 ? "topas:3.2.0" : cl_opts["--mc_code"][0];
    
    std::unique_ptr<rti::treatment_session<phase_space_type>> 
    treatment_session_(new rti::treatment_session<phase_space_type>(cl_opts["--rti"][0], machine_name, mc_code));
    
    treatment_session_->summary();

    rti::vec3<phase_space_type>(1.0,1.0,0.1);

    /*
      position displacement vector
    */
    rti::vec3<float>   p_xyz_(0,0,0);  
    if( cl_opts["--pxyz"].size() >=1 ){
        p_xyz_.x = std::stof(cl_opts["--pxyz"][0]);
        p_xyz_.y = std::stof(cl_opts["--pxyz"][1]);
        p_xyz_.z = std::stof(cl_opts["--pxyz"][2]);
    }

    //local rotation to account gantry rotation, couch rotation, couch pitch (better check).
    rti::mat3x3<float> r_xyz_(0,0,0);  
    if( cl_opts["--rxyz"].size() >=1 ){
        rti::vec3<float> rad ;
        rad.x = M_PI*std::stof(cl_opts["--rxyz"][0])/180.0;
        rad.y = M_PI*std::stof(cl_opts["--rxyz"][1])/180.0;
        rad.z = M_PI*std::stof(cl_opts["--rxyz"][2])/180.0;
        r_xyz_ = rti::mat3x3<float>(rad.x, rad.y, rad.z);
    }

    //particles per history
    double particles_per_history = -1;  //generate 1 history per beamlet
    if (cl_opts["--pph"].size()>=1) particles_per_history = std::stod(cl_opts["--pph"][0]);

    //source to isocenter distance (mm)
    float sid = 0.0 ;
    if (cl_opts["--sid"].size()>=1) sid = std::stof(cl_opts["--sid"][0]);

    //coordinate system: translational & rotational
    rti::coordinate_transform<phase_space_type> p_coord;

    //Objects set by either beam name or beam number
    rti::beamsource<phase_space_type> beam_src;
    rti::beamline<phase_space_type>   beam_line;
    
    if( cl_opts["--bname"].size() >=1 ){
    //if( cl_opts["--bname"].size() >=1 ){
        std::string bname = cl_opts["--bname"][0];
        p_coord   = treatment_session_->get_coordinate(bname);
        beam_src  = treatment_session_->get_beamsource(bname, p_coord, particles_per_history, sid);
        beam_line = treatment_session_->get_beamline(bname);
    }

    if( cl_opts["--bnumber"].size() >=1 ){
    //if( cl_opts["--bnumber"].size() >=1 ){
        int bnb   = std::stoi(cl_opts["--bnumber"][0]);
        p_coord   = treatment_session_->get_coordinate(bnb);
        beam_src  = treatment_session_->get_beamsource(bnb, p_coord, particles_per_history, sid);
        beam_line = treatment_session_->get_beamline(bnb);
    }

    p_coord.dump();
    std::cout<<"gantry rotation: "<<p_coord.gantry.y<< std::endl;
    std::cout<<"couch rotation: "<<p_coord.patient_support.z<< std::endl;
    std::cout<<"colli rotation: "<<p_coord.collimator.z<< std::endl;

    //History generation
    //default spot-generation : no generation
    size_t start_id = 0; 
    size_t end_id   = beam_src.total_beamlets(); 

    const unsigned int beamlet_opt_nb = cl_opts["--beamlets"].size();
    //const unsigned int beamlet_opt_nb = cl_opts["--beamlets"].size();
    assert(beamlet_opt_nb<3);

    switch(beamlet_opt_nb){
    case 0:
        break;
    case 1:
    {
        start_id = std::stoi(cl_opts["--beamlets"][0]); 
        end_id = start_id + 1 ;
    };
    break;
    case 2:
    {
        start_id = std::stoi(cl_opts["--beamlets"][0]);
        end_id = std::stoi(cl_opts["--beamlets"][1]);    
    };
    break;
    default :
        break;
    }
    assert(end_id > start_id);


    /*
    * To generate 2.6M history, it took about 2 sec in CPU & memory usage was < 7 MB
    * The memory usage doesn't change with number of histories.
    */
    std::cout<< "The size of source object: " 
             << sizeof(beam_src) 
             << ", total history:"
             << beam_src.total_histories()<< std::endl;

    std::queue< std::tuple<phase_space_type, rti::vec3<phase_space_type>, rti::vec3<phase_space_type>> > histories;

    
    for(size_t current_id = start_id ; current_id < end_id ; ++current_id){
        auto beamlet = beam_src[current_id];
	    auto beamlet_distribution = std::get<0>(beamlet);
        size_t histories_of_beamlet = std::get<1>(beamlet);

        while(histories_of_beamlet--){
            auto history = beamlet_distribution();
            histories.push(history);
        }
        
	
	

    }

    std::cout<<"Number of of histories in queue: "<< histories.size() << std::endl;
    const auto size_of_history = (sizeof(phase_space_type)+2*sizeof(rti::vec3<phase_space_type>));
    std::cout<<"Byte size of the queue: "<< histories.size()*size_of_history << std::endl;

    std::ofstream output;

    if( cl_opts["--output"].size() >=1 ){
        
        if( cl_opts["--output"][0].compare("")) output.open(cl_opts["--output"][0]);

        std::ostream& out = cl_opts["--output"][0].compare("") ? output : std::cout;

        while(!histories.empty()){
            auto history   = histories.front();
            auto energy    = std::get<0>(history);
            auto pos       = std::get<1>(history);
            auto momentum  = std::get<2>(history);
            out<< pos.x <<" "<< pos.y <<" "<< pos.z<< " " 
               << momentum.x << " " << momentum.y << " " << momentum.z << " "
               << energy << std::endl;
            histories.pop();
        }
    }

    


    return 0;
}
template<typename T>
void generator(
	rti::beamlet<T>& bl, 
	size_t thread_id, 
	size_t n_histories)
{
    //std::cout<<"thread id: " << thread_id <<", n_histories: " << n_histories << std::endl;
    //std::queue< std::tuple<phase_space_type, rti::vec3<phase_space_type>, rti::vec3<phase_space_type>> > histories;
	while(n_histories--){
	    auto history = bl();
        //histories.push(history);
	}	
        
    /* 
    //how we return this to thread outside?
    while(!histories.empty()){
        auto history   = histories.front();
        auto energy    = std::get<0>(history);
        auto pos       = std::get<1>(history);
        auto momentum  = std::get<2>(history);
	    std::cout<< pos.x <<" "<< pos.y <<" "<< pos.z<< " " 
            << momentum.x << " " << momentum.y << " " << momentum.z << " "
            << energy << std::endl;
            histories.pop();
    }
    */

}

void describe(const rti::dataset* beam){
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


void describe_rtibtr(const rti::dataset* beam){
    std::cout << "======== BEAM info RTIBTR =======" << std::endl;
    std::vector<int> data_int;
    std::vector<int> data_int_tag;
    std::vector<float> data_float;
    std::vector<std::string> data_str;
    
    //beam->get_values("BeamNumber", data_int);
    //std::cout<< "BeamNumber: " << data_int[0] << std::endl;

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
    
    //beam->get_values("VirtualSourceAxisDistances", data_float);
    //std::cout <<"VirtualSourceAxisDistances : " << data_float[0] <<", " << data_float[1] << std::endl;

}

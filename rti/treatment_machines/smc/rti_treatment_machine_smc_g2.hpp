#ifndef RTI_TREATMENT_MACHINE_SMC_GTR1_H
#define RTI_TREATMENT_MACHINE_SMC_GTR1_H

#include <rti/base/rti_treatment_machine_ion.hpp>

namespace rti{

namespace smc{
    
/*
template<typename R>
class gtr1_material_t: public patient_material_t<R> {
public:
    CUDA_HOST_DEVICE
    gtr1_material_t(): patient_material_t<R>(){;}

    CUDA_HOST_DEVICE
    gtr1_material_t(int16_t hu): patient_material_t<R>(hu){;}

    CUDA_HOST_DEVICE
    ~g2_material_t(){;}
};
*/
    
/**
 * The <code>g2</code> class represents beam model for Sumitomo IMPT machine in SMC
 */
template <typename T>
class g2 : public treatment_machine_ion<T> {
    
protected:
    const T cm2mm = 10.0;
    const T mm2cm = 0.1 ;

    const T nA  = 1e-9    ;
    const T sec = 1.0     ;
    const T ms  = 1e-3*sec;
    const T us  = 1e-6*sec;
    const T ns  = 1e-9*sec;
    const T mm  = 1.0     ;
    const T  m  = 1000.0*mm;
    const T const_e = 1.6021766208e-19;

    
    ///< spot, i.e., line-segment holder
    std::map<float,
             std::vector< rti::beam_module_ion::spot >
             > sequence_;

    ///<
    std::vector<int> nb_spots_per_layer_; //maybe we don't need

    ///< tunning model but not used
    //std::string tune_id_;
    
public:
    g2()
    {
        treatment_machine_ion<T>::SAD_[0] = 2597.0 ;
        treatment_machine_ion<T>::SAD_[1] = 2223.0 ;
    }

    ~g2(){;}


    ///
    virtual
    void
    build_line_segment_sequence(const rti::dataset* ds,
                                rti::modality_type  mod)
    {
        /// Initializes containers (temporal)
        std::vector<std::string>  nb_pts(1)    ;
        std::vector<float>        energy(1)    ;
        std::vector<float>        fwhm_xy(2)   ;
        std::vector<std::string>  tune_id(1)   ;
        std::vector<float>        xy           ; 
        std::vector<float>        weight       ; 
        int                       layer_nb = 0 ;

        auto seq_tags = &rti::seqtags_per_modality.at(mod);
        auto ictrl    = (*ds)(seq_tags->at("ctrl"));

        for(auto b : ictrl){
            
            if ( (layer_nb++)%2) continue;

            b->get_values("NominalBeamEnergy", energy);

            //Used for Sumitomo Line Scanning.
            //>Line Spot Tune ID (300B,1090) SH 3  (short string)
            //>Number of Line Scan Spot  Positions (300B,1092) IS 3   (integer string)
            //>Line Scan Position Map (300B,1094) FL 3  (single precision, 4 bytes)
            //   np.frombuffer(c0[0x300b,0x1094].value, dtype=np.float32)
            //>Line Scan Meterset Weights (300B,1096) FL 3 
            //>Line Scanning Spot Size (300B,1098) FL 3 
            //>Number of Line Scan Spot Paintings (300B,109A) IS 3
            //b->print_dataelement( (*b)[gdcm::Tag(0x300b,0x1098)] );

            b->get_values_from_element( (*b)[gdcm::Tag(0x300b,0x1090)], tune_id);
            b->get_values_from_element( (*b)[gdcm::Tag(0x300b,0x1092)], nb_pts);
            b->get_values_from_element( (*b)[gdcm::Tag(0x300b,0x1094)], xy);
            b->get_values_from_element( (*b)[gdcm::Tag(0x300b,0x1096)], weight);
            b->get_values_from_element( (*b)[gdcm::Tag(0x300b,0x1098)], fwhm_xy);

            //std::cout << "Tune ID size: "<< tune_id.size() << ", value: " << tune_id[0] << "\n";
            std::vector< rti::beam_module_ion::spot > layer ;
            for(int j = 0 ; j < std::stoi(nb_pts[0]) ; ++j){
                /*
                std::cout << energy[0] <<", "
                          << xy[j*2] << ", " << xy[j*2+1] << ", "
                          << fwhm_xy[0] << ", " << fwhm_xy[1] << ", "
                          << weight[j] << "\n";
                */
                layer.push_back({energy[0], xy[j*2], xy[j*2+1], fwhm_xy[0], fwhm_xy[1], weight[j]});
            }//per line-segment
            sequence_.insert( std::make_pair(energy[0], layer));
            nb_spots_per_layer_.push_back( std::stoi(nb_pts[0]));
        }//per layer
    
    }

    
    ///
    virtual
    rti::beamsource<T>
    create_beamsource( const rti::dataset* ds,
                       const rti::modality_type  m,
                       const rti::coordinate_transform<T> pcoord,
                       const float scalefactor = -1,
                       const float source_to_isocenter_mm = 390.0)
    {

        std::cout << "SMC beam source\n";
        
        treatment_machine<T>::source_to_isocenter_mm_ = source_to_isocenter_mm;


        std::vector<std::string> scan_mode(1);
        ds->get_values("ScanMode",scan_mode);

        if( !scan_mode[0].compare("LINE") ) std::runtime_error("Only LINE scan is supported!");
        
        ///Let's create a beam sequence
        // -> will create one of member, i.e., "sequence_"
        build_line_segment_sequence(ds, m);

        //< beamsource with line-segment
        rti::beamsource<T> beamsource;

        ///< scan from higher energy to lower energy
        //for (auto iter = sequence_.rbegin(); iter != std::next(sequence_.rbegin()); ++iter) {
        for (auto iter = sequence_.rbegin(); iter != sequence_.rend(); ++iter) {
            //std::cout << iter->first << ": " << iter->second << std::endl;
            const std::vector< rti::beam_module_ion::spot >& positions  = iter->second ;
            
            //std::cout << iter-> first << "\n";
            
            ///< position [i, i+1) or (i, i+1] => will determine layer.size() -1 or not
            ///  meterset from i+1. meterset is not 
            for(size_t i = 0 ; i < positions.size() -1 ; ++i){
                /*
                std::cout << "From (" << positions[i].x  << ", " << positions[i].y << ") "
                          << "To  (" << positions[i+1].x  << ", " << positions[i+1].y << ") "
                          << "MU  "  << positions[i+1].meterset << "\n";
                */
                
                std::array<T, 2> time_on_off = {1.0, 0.0}; //temporal
                
                size_t nb_histories = (scalefactor == -1) ? 1 : this->characterize_history(positions[i+1], scalefactor);

                beamsource.append_beamlet(this->characterize_beamlet(positions[i], positions[i+1]),
                                          nb_histories,
                                          pcoord,
                                          time_on_off[0],
                                          time_on_off[1]);
            }//positions
            
        }//layer

        return beamsource;
    }
    
    virtual
    rti::beamlet<T>
    characterize_beamlet( const rti::beam_module_ion::spot& s0,
                          const rti::beam_module_ion::spot& s1)
    {
        
        auto energy = new rti::norm_1d<T>({s0.e},{0.0});
        rti::vec3<T> dir0(std::atan(s0.x/treatment_machine_ion<T>::SAD_[0]),
                          std::atan(s0.y/treatment_machine_ion<T>::SAD_[1]),
                          -1.0);
        rti::vec3<T> dir1(std::atan(s1.x/treatment_machine_ion<T>::SAD_[0]),
                          std::atan(s1.y/treatment_machine_ion<T>::SAD_[1]),
                          -1.0);

        rti::vec3<T> pos0(0, 0, rti::treatment_machine_ion<T>::source_to_isocenter_mm_);
        pos0.x = (treatment_machine_ion<T>::SAD_[0] - pos0.z) * dir0.x ;
        pos0.y = (treatment_machine_ion<T>::SAD_[1] - pos0.z) * dir0.y ;
	
        rti::vec3<T> pos1(0, 0, rti::treatment_machine_ion<T>::source_to_isocenter_mm_);
        pos1.x = (treatment_machine_ion<T>::SAD_[0] - pos1.z) * dir1.x ;
        pos1.y = (treatment_machine_ion<T>::SAD_[1] - pos1.z) * dir1.y ;

        std::array<T,6> spot_pos_range = {pos0.x, pos1.x, pos0.y, pos1.y, pos0.z, pos1.z};
        std::array<T,6> spot_sigma     = {s0.fwhm_x/float(SIGMA2FWHM),  s0.fwhm_y/float(SIGMA2FWHM), 0.0, 0.0, 0.0, 0.0};
        std::array<T,2> corr           = {0.0, 0.0};

        //spot direction mean is calculated from the sampled position between x0 and x1
        //so passing spot-range to a specialized distribution for handling this case is prefered option.
        auto spot_= new rti::phsp_6d_fanbeam<T>(spot_pos_range, spot_sigma, corr, treatment_machine_ion<T>::SAD_);

        return rti::beamlet<T>(energy, spot_);
    }

    virtual
    size_t
    characterize_history( const rti::beam_module_ion::spot& s0,
                          const rti::beam_module_ion::spot& s1,
                          float scale)
    {
        return 0;
    }


    virtual
    size_t
    characterize_history(const rti::beam_module_ion::spot& s,
                         float scale)
    {

        //This is a madeup.
        /*
        double db_IC2_gain_coeff1  = 3.0*1.23547*153.69; //      #unitless
        double db_IC2_gain_coeff2  = -0.6891;             //    #Exponent of kinetic energy
        double db_ChargePerMU_IC2  = 3.0e-9;              //
        double db_MU_per_IC_charge = 1.0/db_ChargePerMU_IC2 ; //
        //#number of PBS MU's per Columnb in signal from IC2(MU2/IC2_C)
        double db_C_per_Gp = 1.0e+9/6.241e+18;
        double MU2GP = db_C_per_Gp*db_IC2_gain_coeff1*std::pow(s.e, db_IC2_gain_coeff2)*db_MU_per_IC_charge;
        */
        ///TODO
        return 1e9*s.meterset/scale;
    }


    virtual
    rti::beamlet<T>
    characterize_beamlet(const rti::beam_module_ion::spot& s)
    {

        auto energy = new rti::const_1d<T>({s.e},{0});

        rti::vec3<T> dir(std::atan(s.x/treatment_machine_ion<T>::SAD_[0]),
                         std::atan(s.y/treatment_machine_ion<T>::SAD_[1]),
                         -1.0);

        rti::vec3<T> pos(0, 0, rti::treatment_machine_ion<T>::source_to_isocenter_mm_);
        pos.x = (treatment_machine_ion<T>::SAD_[0] - pos.z) * dir.x ;
        pos.y = (treatment_machine_ion<T>::SAD_[1] - pos.z) * dir.y ;

        if ( rti::treatment_machine_ion<T>::source_to_isocenter_mm_ > 0){
            //For energy correction
            //due to the air amount between Isocenter and your beam starting position
        }

        //Define phsp distribution
        std::array<T,6> spot_mean = {pos.x, pos.y, pos.z, dir.x, dir.y, 0};
        std::array<T,6> spot_sigm = {s.fwhm_x/T(SIGMA2FWHM) , s.fwhm_y/T(SIGMA2FWHM), 0, 0, 0, 0};
        std::array<T,2> corr = {0.0, 0.0};
        auto spot_= new rti::phsp_6d<T>(spot_mean, spot_sigm, corr);

        return rti::beamlet<T>(energy, spot_);
    }

    rti::rangeshifter*
    characterize_rangeshifter(
        const rti::dataset* ds,
        rti::modality_type m)
    {
        auto seq_tags = &rti::seqtags_per_modality.at(m);

        //1. rangeshifter sequence
        auto  rs_ds = (*ds)( seq_tags->at("rs")) ;
        assert(rs_ds.size() >=1);

        //2. Snout position from control point 0
        //layer0 for snout position
        std::vector<float> ftmp;
        auto layer0    = (*ds)(seq_tags->at("ctrl"))[0];

        ///< TODO: potential bug for SnoutPosition. it got zero.
        ///  check!!! Snoutposition is FL. get_values may not work properly for FL type value
        layer0->get_values( "SnoutPosition", ftmp);
        std::cout << "SnoutPosition: " << ftmp[0] << "\n";
        
        rti::vec3<float>    lxyz(300.0, 300.0, 0.0) ;
        
        /// TODO !!!!! : temporarily set to 350.0 for all rangeshifter
        //rti::vec3<float>    pxyz(0.0, 0.0, ftmp[0]) ;
        rti::vec3<float>    pxyz(0.0, 0.0, 350.0) ;
        
        rti::mat3x3<float>  rxyz(0.0, 0.0, 0.0) ;

        //3. There must be at least one range shifter sequence
        for(auto r : rs_ds){
            std::vector<std::string> rs_id(0);
            r->get_values("RangeShifterID" , rs_id);

            std::cout<< "RangeShifterID value: " << rs_id[0] << std::endl;
            if ( rs_id[0].compare("SNOUT_DEG_B")){
                lxyz.z = 10.0;
            }else{ //from RTIBTR, needs a space
                lxyz.z = 0;
            }
            assert(lxyz.z>0);
        }

        std::cout<<"Range shifter thickness: " << lxyz.z
                 <<" (mm) and position: " << pxyz.z <<" (mm)" << std::endl;

        //!!!! I assumed no-gap between snout and range shifter
        //Shift down to
        pxyz.z -= lxyz.z ;
        pxyz.dump();

        return new rti::rangeshifter(lxyz, pxyz, rxyz);

    }

    rti::aperture*
    characterize_aperture(const rti::dataset* ds,
                          rti::modality_type m)
    {
        auto xypts = this->characterize_aperture_opening(ds,m);

        auto seq_tags = &rti::seqtags_per_modality.at(m);

        //1. block sequence
        auto  blk_ds = (*ds)( seq_tags->at("blk")) ;
        assert(blk_ds.size() >=1);

        std::vector<float> ftmp;

        //note:
        // 1. Assumed opening points are physical points of geometry.
        //    (there must be tag for this) => this should be taken cared by users
        //apt_lunder -> set_dimension(150.0, 67.5);
        //400 is temporal
        blk_ds[0]->get_values( "BlockThickness", ftmp);
        rti::vec3<float> lxyz(400.0, 400.0, ftmp[0]);

        //I omitted reading patient_side or snout_side
        blk_ds[0]->get_values("IsocenterToBlockTrayDistance", ftmp);
        rti::vec3<float> pxyz(0.0, 0.0, ftmp[0]);

        rti::mat3x3<float> rxyz(0.0, 0.0, 0.0);

        return new rti::aperture(xypts, lxyz, pxyz, rxyz);
    }
    
};

}
}
#endif

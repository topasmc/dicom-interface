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

    ///< Time factor
    const T Vx_max = 30.0 * m/sec  ; //scan_x : maximum scan speed X
    const T Vy_max =  3.0 * m/sec  ; //scan_y : maximum scan speed y
    const T Te     =  2.0 * sec    ; //Time to energy change
    
    ///< spot, i.e., line-segment holder
    std::map<float,
             std::vector< rti::beam_module_ion::spot >
             > sequence_;

    ///<
    std::vector<int> nb_spots_per_layer_; //maybe we don't need


    ///< Beam data
    
    // MUtoParticle = protonPerDoseCorrection * dosePerMUCountCorrection
    struct mu_correction{

        T proton_per_dose_correction ;
        T dose_per_mu_count_correction;
    };
    // SigmaX, SigmaY, SigmaXPrime, SigmaYPrime
    struct spot_correction{
        T cross_line ;
        T in_line    ;
        T angular_cross_line ;
        T angular_in_line    ;
    };
    
    const std::map<T, mu_correction> mu_to_particle_   = {
        {70.0,  { 1.0, 1.0}},
        {80.0,  { 1.12573609032495, 0.989255716854649}},
        {90.0,  { 1.25147616113001, 0.973421729297953 }},
        {100.0, { 1.36888442326936, 0.967281770613755 }},
        {110.0, { 1.48668286253201, 0.958215625815887 }},
        {120.0, { 1.60497205195899, 0.946937840980162 }},
        {130.0, { 1.71741194754422, 0.942685675037711 }},
        {140.0, { 1.82898327045955, 0.940168906626851 }}, 
        {150.0, { 1.94071715123743, 0.931161417057087 }},
        {160.0, { 2.04829230739643, 0.918762676945622 }},
        {170.0, { 2.16168786761159, 0.904569498824145 }},
        {180.0, { 2.27629228444253, 0.888164591949398 }}, 
        {190.0, { 2.39246901674031, 0.876689052268837 }},
        {200.0, { 2.50561983301185, 0.872826195199581 }},
        {210.0, { 2.63593473689952, 0.871540965585644 }},
        {220.0, { 2.75663921459094, 0.859481169160383 }}, 
        {230.0, { 2.89392497566575, 0.8524232713089 }}
    };

    const std::map<T, spot_correction>  spot_to_emittance_  = {
        {70.0 , {7.0,  7.6,  0.0043,  0.0043 }},
        {100.0, {5.1,  6.0,  0.0032,  0.0032 }},
        {150.0, {4.1,  4.4,  0.0024,  0.0024 }},
        {162.0, {3.8,  3.9,  0.0023,  0.0023 }},
        {174.0, {1.4,  1.4,  0.0015,  0.0015 }},
        {182.0, {1.0,  1.0,  0.0014,  0.0014 }},
        {190.0, {1.0,  0.9,  0.0014,  0.0014 }},
        {206.0, {0.94, 1.05, 0.0013,  0.0013 }},
        {230.0, {2.0 , 2.0,  0.00125, 0.00125 }}
    };

    ///< Interpolation function
    std::function<T(T, T, T, T, T)> intpl = 
        [](T x, T x0, T x1, T y0, T y1){return ( x1 == x0)? y0 : y0 + (x-x0)*(y1-y0)/(x1 - x0);};

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


                std::array<T, 2> time_on_off = {1.0, 0.0}; 
                const T Tx = std::abs(positions[i+1].x - positions[i].x) / Vx_max ;
                const T Ty = std::abs(positions[i+1].y - positions[i].y) / Vy_max ;
                
                //Beam_on time > 0 always due to continous raster scan
                //Beam_off time = 0 or Te when changing energy
                time_on_off[0]  = (Tx > Ty) ? Tx : Ty;
                time_on_off[1]  = ( i == (positions.size() - 2) ) ? Te : 0.0 ;
                
                size_t nb_histories = (scalefactor == -1) ? 1 : this->characterize_history(positions[i+1], scalefactor);
                
                std::cout << "From (" << positions[i].x  << ", " << positions[i].y << ") "
                          << "To  (" << positions[i+1].x  << ", " << positions[i+1].y << ") "
                          << "Lx, Ly  (" << std::abs(positions[i+1].x - positions[i].x)
                          << ", "        << std::abs(positions[i+1].y - positions[i].y)<< ") "
                          << "Time  (" << time_on_off[0]  << ", " << time_on_off[1] << ") "
                          << "MU  "  << positions[i+1].meterset << "\n";

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

        ///< MU to Particles
        T mid_x  = spot_to_emittance_.begin() -> second.cross_line ;
        T mid_y  = spot_to_emittance_.begin() -> second.in_line ;
        T mid_xp = spot_to_emittance_.begin() -> second.angular_cross_line ;
        T mid_yp = spot_to_emittance_.begin() -> second.angular_in_line ;
        
        auto spot_up = spot_to_emittance_.lower_bound(s0.e);
        if (spot_up != spot_to_emittance_.begin() ){
            auto spot_down = std::prev(spot_up,1);
            auto down = spot_down->second;
            auto up   = spot_up -> second;

            mid_x  = intpl(s0.e, spot_down->first, spot_up->first, down.cross_line,  up.cross_line);
            mid_y  = intpl(s0.e, spot_down->first, spot_up->first, down.in_line,  up.in_line);
            mid_xp = intpl(s0.e,
                           spot_down->first, spot_up->first,
                           down.angular_cross_line,  up.angular_cross_line);
            mid_yp = intpl(s0.e,
                           spot_down->first, spot_up->first,
                           down.angular_in_line,  up.angular_in_line);

        } 
        std::cerr<<"For " << s0.e <<", "
                 << mid_x <<", "
                 << mid_y <<", "
                 << mid_xp << ", "
                 << mid_yp <<", " <<"\n";
        std::array<T,6> spot_pos_range = {pos0.x, pos1.x, pos0.y, pos1.y, pos0.z, pos1.z};
        std::array<T,6> spot_sigma     = {mid_x,
                                          mid_y,
                                          0.0,
                                          mid_xp,
                                          mid_yp,
                                          0.0};
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
        T mid_dose_correction = mu_to_particle_.begin()->second.proton_per_dose_correction;
        T mid_mu_correction   = mu_to_particle_.begin()->second.dose_per_mu_count_correction;
        
        auto mu_up = mu_to_particle_.lower_bound(s.e);
        if ( mu_up != mu_to_particle_.begin() ){
            auto mu_down = std::prev(mu_up,1);
            auto down = mu_down->second;
            auto up   = mu_up -> second;
        
            mid_dose_correction    = intpl(s.e,
                                           mu_down->first,
                                           mu_up->first,
                                           down.proton_per_dose_correction,
                                           up.proton_per_dose_correction);
        
            mid_mu_correction    = intpl(s.e,
                                         mu_down->first,
                                         mu_up->first,
                                         down.dose_per_mu_count_correction,
                                         up.dose_per_mu_count_correction);
        }
        
        std::cerr<<"For " << s.e <<", "
                 << mid_dose_correction <<", "
                 << mid_mu_correction <<", "
                 << s.meterset << ", "
                 << (s.meterset*mid_dose_correction*mid_mu_correction/scale) <<", " <<"\n";
        return (s.meterset*mid_dose_correction*mid_mu_correction)/scale;
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

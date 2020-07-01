#ifndef RTI_TREATMENT_MACHINE_RBE_1P1_H
#define RTI_TREATMENT_MACHINE_RBE_1P1_H

/// \file
///
/// Treatment machine for Lunder proton center and TOPAS MC code

#include <rti/base/rti_treatment_machine_ion.hpp>

namespace rti{

namespace rbe{

/// \class rbe_1p1
/// represents beam model for Lunder IMPT machine in MGH
template <typename T>
class rbe_1p1 : public treatment_machine_ion<T> {
protected:
    
public: 
    /// Default constructor
    rbe_1p1()
    {
        treatment_machine_ion<T>::SAD_[0] = 1800.0;
        treatment_machine_ion<T>::SAD_[1] = 2000.0;
    }

    ~rbe_1p1(){;}

    /// Characterize MODULATED_SPEC/UNIFORM
    /// \note not yet implemented cause we treat patients with spot-scanning only
    virtual
    size_t
    characterize_history(
        const rti::beam_module_ion::spot& s0, 
        const rti::beam_module_ion::spot& s1, 
        float scale)
    {
        return 0;
    }

    /// Characterize MODULATED_SPEC/UNIFORM
    /// \note not yet implemented cause we treat patients with spot-scanning only
    virtual
    size_t
    characterize_history(
        const rti::beam_module_ion::spot& s, 
        float scale)
    {
        return s.meterset / scale;
    }
    
    
    /// Characterize UNIFORM beam delivery continous scan
    /// define a continuous position distribution in rti_distributions   
    /// this uses "phsp6d_fanbeam" distribution to sample position between x0 to x1.
    virtual
    rti::beamlet<T>
    characterize_beamlet(
        const rti::beam_module_ion::spot& s0,
        const rti::beam_module_ion::spot& s1)
    {
        return rti::beamlet<T>();
    } 
    
    /// Characterize beamlet for MODULATED beam
    virtual
    rti::beamlet<T>
    characterize_beamlet(const rti::beam_module_ion::spot& s)
    {
    
	    //Energy distribution
        // 1. constant energy
        auto energy = new rti::const_1d<T>({s.e},{0});

        // 2. X-Y position at at z
        rti::vec3<T> iso(s.x, s.y, 0);
        rti::vec3<T> beam =
			this->beam_starting_position(iso,
								   rti::treatment_machine_ion<T>::source_to_isocenter_mm_);
        rti::vec3<T> dir = iso - beam;
        dir.normalize(); 
		
        // 3. Complete fluence distribution
        // this samples x, x', y, y', z, z'
        std::array<T,6> spot_mean = {beam.x, beam.y, beam.z, dir.x, dir.y, dir.z};
        std::array<T,6> spot_sigm = {s.fwhm_x/T(2.355), s.fwhm_y/T(2.355), 0.0, 0.0, 0.0, 0.0};
        std::array<T,2> corr      = {0.0, 0.0};
        auto fluence= new rti::phsp_6d<T>(spot_mean, spot_sigm, corr);

        return rti::beamlet<T>(energy, fluence);
    }

    /// Characterize Rangeshifter
    rti::rangeshifter*
    characterize_rangeshifter(
        const rti::dataset* ds,
        rti::modality_type m)
    {
        auto seq_tags = &rti::seqtags_per_modality.at(m);

        //1. rangeshifter sequence
        // get a DICOM tag for given modality.
        auto  rs_ds = (*ds)( seq_tags->at("rs")) ; 
        assert(rs_ds.size() >=1);

        //2. Snout position from control point 0
        std::vector<float> ftmp;
        auto layer0    = (*ds)(seq_tags->at("ctrl"))[0]; //layer0 for snout position
        layer0->get_values( "SnoutPosition", ftmp);

        rti::vec3<float>    lxyz(300.0, 300.0, 0.0) ;
        rti::vec3<float>    pxyz(0.0, 0.0, ftmp[0]) ;
        rti::mat3x3<float>  rxyz(0.0, 0.0, 0.0) ;

        //3. range shifter sequence
        for(auto r : rs_ds){
            std::vector<std::string> rs_id(0);
            r->get_values("RangeShifterID" , rs_id);

            if ( !rs_id[0].compare("10mm")){
                lxyz.z += 10.0;
            }else if ( !rs_id[0].compare("20mm")){
                lxyz.z += 20.0;
            }else if ( !rs_id[0].compare("40mm")){
                lxyz.z += 40.0;
            }
            assert(lxyz.z>0);
        }

        pxyz.z -= (lxyz.z + 5.0) ; //gap is 5 mm between snout downstream to rangeshifter upstream
        std::cout<<"Range shifter thickness: " << lxyz.z 
        << " (mm) and position: " << pxyz.z <<" (mm)" << std::endl;

        return new rti::rangeshifter(lxyz, pxyz, rxyz);   
    }

    /// Characterize aperture
    rti::aperture*
    characterize_aperture(
        const rti::dataset* ds,
        rti::modality_type m
    ){
        //1. aperture opening points.
        auto xypts = this->characterize_aperture_opening(ds,m);

        //2. block sequence
        auto seq_tags = &rti::seqtags_per_modality.at(m);
        auto  blk_ds = (*ds)( seq_tags->at("blk")) ; 
        assert(blk_ds.size() >=1);

        std::vector<float> ftmp;

        blk_ds[0]->get_values( "BlockThickness", ftmp);
        rti::vec3<float> lxyz(400.0, 400.0, ftmp[0]); 

        blk_ds[0]->get_values("IsocenterToBlockTrayDistance", ftmp);
        rti::vec3<float> pxyz(0.0, 0.0, ftmp[0]+lxyz.z*0.5);

        rti::mat3x3<float> rxyz(0.0, 0.0, 0.0);

        return new rti::aperture(xypts, lxyz, pxyz, rxyz);
    }
};

}
}
#endif

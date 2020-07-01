#ifndef RTI_TREATMENT_SESSION_HPP
#define RTI_TREATMENT_SESSION_HPP

/// \file
///
/// Treatment session

#include <type_traits>
#include <algorithm>

#include "gdcmReader.h"

#include <rti/base/rti_dataset.hpp>
#include <rti/base/rti_utils.hpp>
#include <rti/base/rti_treatment_machine.hpp>
#include <rti/base/rti_ct.hpp>
//generic treatment machine
#include <rti/base/rti_treatment_machine_pbs.hpp>


//Custom treatment machines
#include <rti/treatment_machines/rbe/rbe_1p1.hpp>

namespace rti {

/// \class treatment_session
/// \tparam T type of phase-space variables
/// Reads RT-Ion file, creates treatment_machine,
/// and returns machine objects, geometry, source, and coordinate system.
/// treatment_session is an entry point to RT-Ion interface to a MC engine.
/// \note we are considering to include patient, and dosegrid here.
template<typename T>
class treatment_session {
protected:

    ///< Sequence tag dictionary for modality specific
    const std::map<
        const std::string,
        const gdcm::Tag
        >* seq_tags_ ;

    ///< RT Modality type, e.g., RTPLAN, RTRECORD, IONPLAN, IONRECORD
    rti::modality_type  mtype_ ;

    ///< machine name, e.g., institution_name:machine_name
    std::string machine_name_ ;

    rti::treatment_machine<T>* tx_machine_ = nullptr;

    ///< top level DICOM dataset, either RTIP or RTIBTR
    rti::dataset* rti_ds_ = nullptr;

public:

    /// Constructs treatment machine based on DICOM or specific file name.
    /// It reads in RT file recursively and construct a dataset tree
    /// Depending on RTIP or RTIBTR, it copies a propriate DICOM tag dictionaries
    /// (seqtags_per_modality).
    ///
    /// \param file_for_tx_machine : RTPLAN, IONPLAN, RTRECORD, IONRECORD
    /// Currently IONPLAN and IONRECORD are supported only.
    treatment_session(
        std::string file_for_tx_machine, //for beamline and source,
        std::string m_name  = "",
        std::string mc_code = "topas:3.x")
    {

        gdcm::Reader reader;

        reader.SetFileName(file_for_tx_machine.c_str());
        
        const bool is_file_valid = reader.Read();
        if(!is_file_valid) throw std::runtime_error("Invalid DICOM file is given to treatment_session.");

        gdcm::MediaStorage ms;
        ms.SetFromFile(reader.GetFile());

        switch(ms){
        case gdcm::MediaStorage::RTPlanStorage :
            mtype_ = RTPLAN;
            break;
        //case gdcm::MediaStorage::RTBeamsTreatmentRecordStorage :
        //gdcm doesn't have definition of RTBeamsTreatmentRecordStorge
        //case gdcm::MediaStorage::RTBeamsTreatmentRecordStorage : //1.2.840.10008.5.1.4.1.1.481.4
        //    mtype_ = RTRECORD;
        //    break;
        case gdcm::MediaStorage::RTIonPlanStorage :
            mtype_ = IONPLAN;
            break;
        case gdcm::MediaStorage::RTIonBeamsTreatmentRecordStorage :
            mtype_ = IONRECORD;
            break;
        default:
            throw std::runtime_error("treatment_session does not supports given RTMODALITY");
        }

        rti_ds_   = new rti::dataset(reader.GetFile().GetDataSet(), true);
        seq_tags_ = &rti::seqtags_per_modality.at(mtype_);

        ///< machine name is set
        if( !m_name.compare("") ){
            ///< we assume machine name is in 1-st beam sequence
            ///< This part doesn't work for RTIBTR
            const rti::dataset* beam_ds;
            if (mtype_ == IONPLAN){
                beam_ds = (*rti_ds_)( seq_tags_->at("beam") )[0];
            }else if( mtype_ == IONRECORD){
                beam_ds = (*rti_ds_)( seq_tags_->at("machine") )[0];
            }else{
                throw std::runtime_error("treatment_session can't find machine name and institute.");
            }

            std::vector<std::string> in ; beam_ds->get_values("InstitutionName", in);
            std::vector<std::string> bn ; beam_ds->get_values("TreatmentMachineName" , bn);

            if(in.size()>0 && bn.size()>0)
            {
                machine_name_ = rti::trim_copy(in[0]) + ":" + rti::trim_copy(bn[0]);
            }
            else if(in.size()>0)
            {
                machine_name_ = rti::trim_copy(in[0]);
            }
            else if(bn.size()>0)
            {
                machine_name_ = rti::trim_copy(bn[0]);
            }

        }else{
            machine_name_ = m_name;
			std::cout<<machine_name_<<"\n";
        }
	
        if (!this->create_machine(machine_name_, mc_code)){
            std::runtime_error("No MC machine is registered for "+machine_name_);
        }

    }


    /// Creates rti::machine and return true for sucessful creation or false.
    /// \param machine_name for machine name
    /// \param mc_code for mc engine, e.g., code:version
    /// \note it takes itype_ member variables. caution, itype_ shouldn't be changed after creation.
    /// \note we have a branch for machines based on "string" comparison.
    ///  Looking for better way to determine during 'ideally' pre-processing.
    /// type_traits allows to branch the logic flow based on the type of variables.
    bool
    create_machine
	(std::string machine_name,
     std::string mc_code)
    {
		if(tx_machine_){
				throw std::runtime_error("Preexisting machine.");
		}
        std::cout<<"machine_name: "<<machine_name<<", mc_code: "<<mc_code<<std::endl;
		const size_t deli = machine_name.find(":");
		std::string site  = machine_name.substr(0, deli);
        std::transform(site.begin(),site.end(),site.begin(),::tolower);
        std::string model = machine_name.substr(deli+1,machine_name.size());
        
        if( site.compare("pbs")) std::transform(model.begin(), model.end(), model.begin(), ::tolower);
        
		std::cout<<site<< ":" << model <<"\n";

		if( !site.compare("pbs") ){
			//Generic PBS beam model
			//expecting file
			std::cout<<"Creating a generic PBS machine from : " << model << "\n";
			tx_machine_ = new rti::pbs<T>(model);
			return true;
		}else if(!site.compare("rbe") ){
			if ( !model.compare("1.1") ){
				tx_machine_ = new rti::rbe::rbe_1p1<T>;
				return true;
			}
       	    else{
				throw std::runtime_error("Valid machine is not available.");
			}

		}else{
			throw std::runtime_error("Valid site is not available.");
		}

        return false;
    }

    /// Default destructor
    ~treatment_session(){
	delete tx_machine_;
        delete rti_ds_;
    }

   /*
    const
    std::shared_ptr<rti::dosegrid>
    get_dosegrid(std::shared_ptr<const rti::patient> p)
    {   dosegrid_->set_from_patient(p);
        return dosegrid_;
    }
    */

   /// Returns a list of beam names present in the plan file.
   std::vector<std::string>
   get_beam_names() {
       auto beam_sequence = (*rti_ds_)( seq_tags_->at("beam") );
       std::vector<std::string> beam_names;
       for (const auto &beam : beam_sequence) {
            std::vector<std::string> tmp;
            beam->get_values("BeamName", tmp);
            beam_names.push_back(tmp[0]);
       }
       return beam_names;
   }


    /// Search and return a beam (DICOM dataset) in BeamSequence for given beam name
    /// \param bnm for beam name
    /// \return rti::dataset* constant pointer.
    const
    rti::dataset*
    get_beam_dataset(std::string bnm)
    {
        auto bseq = (*rti_ds_)( seq_tags_->at("beam") );
        for(auto i : bseq){
            std::vector<std::string> bn ;
            i->get_values("BeamName" , bn); //works for RTIP & RTIBTR
            if(bn.size()==0) continue;
            if (bnm == rti::trim_copy(bn[0])) {
                return i ;
            }
        }
        throw std::runtime_error("Invalid beam name.");
    }

    /// Search and return a beam (DICOM dataset) in BeamSequence for given beam number
    /// \param bnm for beam number
    /// \return rti::dataset* constant pointer.
    const
    rti::dataset*
    get_beam_dataset(int bnb)
    {
        auto bseq = (*rti_ds_)( seq_tags_->at("beam") );
        for(auto i : bseq){
            //i->dump();
            std::vector<int> bn ;
            i->get_values( "BeamNumber" , bn); //works only for RTIP
            assert(bn.size()>0);
            if (bnb == bn[0]){
                return i ;
            }
        }
        throw std::runtime_error("Invalid beam number.");
    }


    /// Get beamline object for given beam id, e.g, beam name or beam number
    /// \param beam_id for beam number or beam name
    /// \return beamline object.
    template<typename S>
    rti::beamline<T>
    get_beamline(S beam_id)
    {
        return tx_machine_->create_beamline(
                    this->get_beam_dataset(beam_id),
                    mtype_);
    }

    /// Gets beam source object for given beam id, e.g, beam name or beam number
    /// \param beam_id for beam number or beam name
    /// \param coord   for coordinate transformation
    /// \param sid     for source to isocenter distance in mm
    /// \param scale   for calculating number of histories to be simulated from beamlet weight
    /// \return beamsource object.
    template<typename S>
    rti::beamsource<T>
    get_beamsource
	( S beam_id,
      const rti::coordinate_transform<T> coord,
      float scale,
      T     sid)
    {
        return  tx_machine_->create_beamsource(
                    this->get_beam_dataset(beam_id),
                    mtype_,
                    coord,
                    scale,
                    sid);
    }

   /// Gets time line object for given beam id, e.g., beam name or number
   /// \param beam_id
   template<typename S>
   std::map<T, int32_t>
   get_timeline(S beam_id)
   {
       return tx_machine_->create_timeline(
           this->get_beam_dataset(beam_id),
           mtype_);
   }

    /// Gets beam coordinate object for given beam id, e.g, beam name or beam number
    /// \param beam_id for beam number or beam name
    /// \return beamsource object.
    template<typename S>
    rti::coordinate_transform<T>
    get_coordinate(S beam_id){
        return tx_machine_->create_coordinate_transform(
            this->get_beam_dataset(beam_id),
            mtype_);
    }

    /// Summarize plan
    /// \param bnb for beam number
    void
    summary(void){
        //plan_ds
        rti_ds_->dump();
    }

    /// Return ion type
    /// \return rti::m_type_
    rti::modality_type
    get_modality_type(void){
        return     mtype_  ;
    }

};

}

#endif

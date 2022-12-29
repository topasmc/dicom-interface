
#ifndef TsRTIonSource_hh
#define TsRTIonSource_hh

#include <memory>
#include <queue>

#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"

#include "globals.hh"

#include "TsSource.hh"
#include "TsPrimaryParticle.hh"

#include <rti/base/rti_treatment_session.hpp>

class TsRTIonSource : public TsSource {
protected:
    
    G4ParticleDefinition* particle_definition_ ;

    std::unique_ptr<rti::treatment_session<float>> rti_session_  ; 

    rti::beamsource<float> beam_source_;

    //current can be replaced by start. however, we keep this.
    struct {
        int start;
        int stop ;
        int      current; //beamlet id, where all histories to this beamlet were passed
        long int history;    //number of histories to be passed to generaters so far.
    } counter_ ;

    void ResolveParameters();

    //Get corrdinate system first
    rti::coordinate_transform<float> rt_coordinate;

    void ExportDICOMCoordinateToParameters(
        rti::coordinate_transform<float>& p
    );

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

public:
    TsRTIonSource(TsParameterManager* pM, TsSourceManager* psM, G4String sourceName);
    ~TsRTIonSource();
    
    //function object, which will be coupled one of 
    //  : ReadHistoriesPerBeamlet for spot-scanning
    //  : ReadHistoriesPerBCM     for IBA's DS
    std::function<G4bool(std::queue<TsPrimaryParticle>*)> ReadHistories;

    //G4bool ReadSomeDataFrombeam(std::queue<TsPrimaryParticle>* );
    G4bool ReadHistoriesPerBeamlet(std::queue<TsPrimaryParticle>* );

    void UpdateForNewRun(G4bool){;}
    
};

#endif

//RTI beamlet for spot-scanning, rasterscanning, continuous scanning
#ifndef TsRTIonSourceGenerator_hh
#define TsRTIonSourceGenerator_hh

#include "globals.hh"

#include "TsVGenerator.hh"
#include "TsPrimaryParticle.hh"

#include <rti/base/rti_treatment_session.hpp>

class TsParameterManager;
class TsGeometryManager ;
class TsGeneratorManager;

class TsRTIonSourceGenerator : public TsVGenerator {
protected:
    std::queue<TsPrimaryParticle>* fParticleBuffer;
    void ResolveParameters();
    void GeneratePrimaries(G4Event* );

public:
    TsRTIonSourceGenerator(
        TsParameterManager* pM, 
        TsGeometryManager* gM, 
        TsGeneratorManager* psM, 
        G4String sourceName);

    ~TsRTIonSourceGenerator();

};

#endif

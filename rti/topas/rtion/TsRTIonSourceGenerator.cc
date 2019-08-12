// Particle Generator for TsRTIonSource
// This code implementation is the intellectual property of the MGH.
// this is an user extension for RT-Ion plan (RTIP) and beam treatment record (RTIBTR)
// by Jungwook Shin jshin13@mgh.harvard.edu

#include <iostream>

//#include <vector>
//#include <tuple>
//#include <array>

#include <map>
#include <vector>
#include <array>
#include <tuple>
#include <algorithm>

#include "G4SystemOfUnits.hh"

//Headers for TOPAS
#include "TsParameterManager.hh"

//Headers for extension
#include "TsRTIonSource.hh"
#include "TsRTIonSourceGenerator.hh"


TsRTIonSourceGenerator::TsRTIonSourceGenerator(
    TsParameterManager* pM, 
    TsGeometryManager* gM, 
    TsGeneratorManager* psM, 
    G4String n)
    :TsVGenerator(pM, gM, psM,n)
{                                   
    fParticleBuffer = new std::queue<TsPrimaryParticle>;
    this->ResolveParameters();
}

TsRTIonSourceGenerator::~TsRTIonSourceGenerator() {;}

void TsRTIonSourceGenerator::ResolveParameters(){
    //TsVGenerator::ResolveParameters();
    //TODO: beam number? or beam name?, e.g., beam number tentatively
    CacheGeometryPointers();
}

void TsRTIonSourceGenerator::GeneratePrimaries(G4Event *evt){
    
    if(CurrentSourceHasGeneratedEnough()) return;

    if (fParticleBuffer->empty()){
        ((TsRTIonSource*) fPs)->ReadHistories(fParticleBuffer);
    }

    //IMPORTANT: assign histories from the buffer to single event.
    //this reduce overhead to create G4Event per history
    while( !fParticleBuffer->empty() ){
        
        TsPrimaryParticle p = fParticleBuffer->front();

        //TransformPrimaryForComponent is a part of topas version higher 3.1
        //TransformPrimaryForComponent(&p);

        GenerateOnePrimary(evt, p); //in this loop, topas add p to vertex
        AddPrimariesToEvent(evt);
        fParticleBuffer->pop();
    }

}
#ifndef __ParticlePhiFilter_h__
#define __ParticlePhiFilter_h__

#include <vector>
#include <string>
#include "StarGenerator/FILT/StarFilterMaker.h"
#include "StarGenerator/EVENT/StarGenEvent.h"

class StarGenParticleMaster;
class StarGenParticle;
class StarGenEvent;

class ParticlePhiFilter : public StarFilterMaker{
public:
    ParticlePhiFilter();
    ~ParticlePhiFilter();

    Int_t Filter(StarGenEvent* event);
private:
    ClassDef(ParticlePhiFilter, 1)
};

#endif

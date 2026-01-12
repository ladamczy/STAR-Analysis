#ifndef __ParticleKstarFilter_h__
#define __ParticleKstarFilter_h__

#include <vector>
#include <string>
#include "StarGenerator/FILT/StarFilterMaker.h"
#include "StarGenerator/EVENT/StarGenEvent.h"

class StarGenParticleMaster;
class StarGenParticle;
class StarGenEvent;

class ParticleKstarFilter : public StarFilterMaker{
public:
    ParticleKstarFilter();
    ~ParticleKstarFilter();

    Int_t Filter(StarGenEvent* event);
private:
    ClassDef(ParticleKstarFilter, 1)
};

#endif

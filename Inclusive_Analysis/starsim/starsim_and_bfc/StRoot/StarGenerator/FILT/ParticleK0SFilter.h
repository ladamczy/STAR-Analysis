#ifndef __ParticleK0SFilter_h__
#define __ParticleK0SFilter_h__

#include <vector>
#include <string>
#include "StarGenerator/FILT/StarFilterMaker.h"
#include "StarGenerator/EVENT/StarGenEvent.h"

class StarGenParticleMaster;
class StarGenParticle;
class StarGenEvent;

class ParticleK0SFilter : public StarFilterMaker{
public:
    ParticleK0SFilter();
    ~ParticleK0SFilter();

    Int_t Filter(StarGenEvent* event);
private:
    ClassDef(ParticleK0SFilter, 1)
};

#endif

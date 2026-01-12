#ifndef __ParticleLambda0Filter_h__
#define __ParticleLambda0Filter_h__

#include <vector>
#include <string>
#include "StarGenerator/FILT/StarFilterMaker.h"
#include "StarGenerator/EVENT/StarGenEvent.h"

class StarGenParticleMaster;
class StarGenParticle;
class StarGenEvent;

class ParticleLambda0Filter : public StarFilterMaker{
public:
    ParticleLambda0Filter();
    ~ParticleLambda0Filter();

    Int_t Filter(StarGenEvent* event);
private:
    ClassDef(ParticleLambda0Filter, 1)
};

#endif

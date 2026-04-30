#ifndef __GeneralHadronFilter_h__
#define __GeneralHadronFilter_h__

#include <vector>
#include <string>
#include "StarGenerator/FILT/StarFilterMaker.h"
#include "StarGenerator/EVENT/StarGenEvent.h"

class StarGenParticleMaster;
class StarGenParticle;
class StarGenEvent;

class GeneralHadronFilter : public StarFilterMaker{
public:
    GeneralHadronFilter();
    ~GeneralHadronFilter();

    Int_t Filter(StarGenEvent* event);
private:
    ClassDef(GeneralHadronFilter, 1)
};

#endif

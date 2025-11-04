#ifndef __StrangeHadronFilter_h__
#define __StrangeHadronFilter_h__

#include <vector>
#include <string>
#include "StarGenerator/FILT/StarFilterMaker.h"
#include "StarGenerator/EVENT/StarGenEvent.h"

class StarGenParticleMaster;
class StarGenParticle;
class StarGenEvent;

class StrangeHadronFilter : public StarFilterMaker{
public:
    StrangeHadronFilter();
    ~StrangeHadronFilter();

    Int_t Filter(StarGenEvent* event);
private:
    ClassDef(StrangeHadronFilter, 1)
};

#endif

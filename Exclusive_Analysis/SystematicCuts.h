#ifndef SYSTEMATIC_CUTS_H
#define SYSTEMATIC_CUTS_H

#include <string>
#include <vector>
#include <iostream>

// Holds all variable cut values for one systematic variation.
// Cuts that are fixed by physics / detector (e.g. BEAM_ENERGY, MASS_PION)
// are NOT included here — only the ones we actually vary.
struct CutConfig {
    std::string name;

    // SC 4.1 — K0S mass window
    double massWinLow;
    double massWinHigh;

    // SC 4.2 — missing pT
    double ptMissMax;

    // SC 5 — topological (4-TOF case)
    double dcaBeamline4;
    double dcaDaughters4;
    double dcaBeamline23;
    double dcaDaughters23;
    double decayLength;
    double cosPA;       // COS_PNT_ANG  (4-TOF and N=5)
    double cosPA23;     // COS_PNT_ANG for N=2/3 (tighter)

    // SC 7 — kinematic correlations
    double ksiCorrMax;      // SC 7.1
    double etaCorrMax;      // SC 7.2
    double separateCorrMax; // SC 7.3 (corrE, corrW)

    // SC 8 — decay vertex positions
    double zmeanMax;
    double zdiffMax;

    // SC 9 — cos(theta*)
    double cosThetaStarMax;

    // Default constructor = nominal thesis values
    CutConfig(std::string n = "nominal") : name(n) {
        massWinLow      = 0.48;
        massWinHigh     = 0.52;
        ptMissMax       = 0.15;
        dcaBeamline4    = 2.5;
        dcaDaughters4   = 2.5;
        dcaBeamline23   = 1.5;
        dcaDaughters23  = 1.5;
        decayLength     = 3.0;
        cosPA           = 0.925;
        cosPA23         = 0.95;
        ksiCorrMax      = 0.004;
        etaCorrMax      = 0.9;
        separateCorrMax = 0.004;
        zmeanMax        = 80.0;
        zdiffMax        = 15.0;
        cosThetaStarMax = 0.8;
    }

    void Print() const {
        std::cout << "Config: " << name << std::endl;
        std::cout << "  SC4.1  mass window     : [" << massWinLow << ", " << massWinHigh << "] GeV" << std::endl;
        std::cout << "  SC4.2  ptMiss max      : " << ptMissMax << " GeV" << std::endl;
        std::cout << "  SC5    dcaBeamline4    : " << dcaBeamline4 << " cm" << std::endl;
        std::cout << "  SC5    dcaDaughters4   : " << dcaDaughters4 << " cm" << std::endl;
        std::cout << "  SC5    dcaBeamline23   : " << dcaBeamline23 << " cm" << std::endl;
        std::cout << "  SC5    dcaDaughters23  : " << dcaDaughters23 << " cm" << std::endl;
        std::cout << "  SC5    decayLength     : " << decayLength << " cm" << std::endl;
        std::cout << "  SC5    cosPA           : " << cosPA << std::endl;
        std::cout << "  SC5    cosPA23         : " << cosPA23 << std::endl;
        std::cout << "  SC7.1  ksiCorr max     : " << ksiCorrMax << std::endl;
        std::cout << "  SC7.2  etaCorr max     : " << etaCorrMax << std::endl;
        std::cout << "  SC7.3  separateCorr max: " << separateCorrMax << std::endl;
        std::cout << "  SC8.1  zmean max       : " << zmeanMax << " cm" << std::endl;
        std::cout << "  SC8.2  zdiff max       : " << zdiffMax << " cm" << std::endl;
        std::cout << "  SC9    cosThetaStar max: " << cosThetaStarMax << std::endl;
    }
};

inline std::vector<CutConfig> CreateSystematicConfigs() {
    std::vector<CutConfig> configs;

    // ---- 0: Nominal ----
    configs.push_back(CutConfig("nominal"));

    // ---- SC 4.1: K0S mass window ----
    CutConfig massWide("massWin_wide");
    massWide.massWinLow  = 0.47;
    massWide.massWinHigh = 0.53;
    configs.push_back(massWide);

    CutConfig massTight("massWin_tight");
    massTight.massWinLow  = 0.49;
    massTight.massWinHigh = 0.51;
    configs.push_back(massTight);

    // ---- SC 4.2: ptMiss ----
    CutConfig ptMissLow("ptMiss_low");
    ptMissLow.ptMissMax = 0.12;
    configs.push_back(ptMissLow);

    CutConfig ptMissHigh("ptMiss_high");
    ptMissHigh.ptMissMax = 0.18;
    configs.push_back(ptMissHigh);

    // ---- SC 5: DCA beamline (4-TOF case) ----
    CutConfig dcaBL4tight("dcaBeamline4_tight");
    dcaBL4tight.dcaBeamline4 = 2.0;
    configs.push_back(dcaBL4tight);

    CutConfig dcaBL4loose("dcaBeamline4_loose");
    dcaBL4loose.dcaBeamline4 = 3.0;
    configs.push_back(dcaBL4loose);

    // ---- SC 5: DCA daughters (4-TOF case) ----
    CutConfig dcaDau4tight("dcaDaughters4_tight");
    dcaDau4tight.dcaDaughters4 = 2.0;
    configs.push_back(dcaDau4tight);

    CutConfig dcaDau4loose("dcaDaughters4_loose");
    dcaDau4loose.dcaDaughters4 = 3.0;
    configs.push_back(dcaDau4loose);

    // ---- SC 5: cos pointing angle ----
    CutConfig cosPAtight("cosPA_tight");
    cosPAtight.cosPA   = 0.94;
    cosPAtight.cosPA23 = 0.96;
    configs.push_back(cosPAtight);

    CutConfig cosPAloose("cosPA_loose");
    cosPAloose.cosPA   = 0.91;
    cosPAloose.cosPA23 = 0.94;
    configs.push_back(cosPAloose);

    // ---- SC 5: decay length ----
    CutConfig decLenTight("decayLength_tight");
    decLenTight.decayLength = 2.5;
    configs.push_back(decLenTight);

    CutConfig decLenLoose("decayLength_loose");
    decLenLoose.decayLength = 3.5;
    configs.push_back(decLenLoose);

    // ---- SC 7.1: ksi correlation ----
    CutConfig ksiTight("ksiCorr_tight");
    ksiTight.ksiCorrMax = 0.003;
    configs.push_back(ksiTight);

    CutConfig ksiLoose("ksiCorr_loose");
    ksiLoose.ksiCorrMax = 0.005;
    configs.push_back(ksiLoose);

    // ---- SC 7.2: eta correlation ----
    CutConfig etaTight("etaCorr_tight");
    etaTight.etaCorrMax = 0.7;
    configs.push_back(etaTight);

    CutConfig etaLoose("etaCorr_loose");
    etaLoose.etaCorrMax = 1.1;
    configs.push_back(etaLoose);

    // ---- SC 7.3: separate corr (corrE, corrW) ----
    CutConfig sepTight("separateCorr_tight");
    sepTight.separateCorrMax = 0.003;
    configs.push_back(sepTight);

    CutConfig sepLoose("separateCorr_loose");
    sepLoose.separateCorrMax = 0.005;
    configs.push_back(sepLoose);

    // ---- SC 8.1: zmean ----
    CutConfig zmeanTight("zmean_tight");
    zmeanTight.zmeanMax = 70.0;
    configs.push_back(zmeanTight);

    CutConfig zmeanLoose("zmean_loose");
    zmeanLoose.zmeanMax = 90.0;
    configs.push_back(zmeanLoose);

    // ---- SC 8.2: zdiff ----
    CutConfig zdiffTight("zdiff_tight");
    zdiffTight.zdiffMax = 12.0;
    configs.push_back(zdiffTight);

    CutConfig zdiffLoose("zdiff_loose");
    zdiffLoose.zdiffMax = 18.0;
    configs.push_back(zdiffLoose);

    // ---- SC 9: cos(theta*) ----
    CutConfig cosThtight("cosThetaStar_tight");
    cosThtight.cosThetaStarMax = 0.7;
    configs.push_back(cosThtight);

    CutConfig cosThloose("cosThetaStar_loose");
    cosThloose.cosThetaStarMax = 0.9;
    configs.push_back(cosThloose);

    return configs;
}

#endif // SYSTEMATIC_CUTS_H
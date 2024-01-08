#ifndef Util_h
#define Util_h

   // Enumerations - very helpful and convenient !
   enum STUDYMAP { kMAINANA = 0, kEMBEDQA, kVERTEXSTUDY, kTOFQA, kTRIGEFF, kFULLZB, kELASTICANA, kRPMCANA, kALIGNMENT, nStudies };
   enum SIDE { E=0, East=0, W=1, West=1, nSides };
   enum RPPORIENTATIONS { Up=0, Down=1, nRpOrientations };
   enum RPCONFIGURATION { IT=0, ET=1, nRpConfigurations};
   enum XY_COORDINATE { X = 0, Y, Z, nCoordinates, nXYCoordinates = Z};
   enum TRIGGER_ID { UPC_v1, UPC_v2, SDT, ET_v1, ET_v2, CPT2_v1, CPT2_v2, CPT2noBBCL, Zerobias, JPsiHTTP_v1, 
   JPsiHTTP_v2, JPsiHTTP_v3, SDT_RHICf, ET_RHICf, CPT2_RHICf, CPT2noBBCL_RHICf_v1, CPT2noBBCL_RHICf_v2, nTriggers };
   enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };
   enum PLANE_ID {A, B, C, D, nPlanes };
   enum BRANCH_ID { EU, ED, WU, WD, nBranches };
   enum BRANCHES_CONFIGURATION_ID { CONF_EU_WU=0, CONF_ED_WD, CONF_EU_WD, CONF_ED_WU, nBranchesConfigurations };
   enum ARM_ID { EU_WD, ED_WU, nArms };
   enum STATION_ID { E1, E2, W1, W2, nStations };
   enum STATION_ORDER { RP1, RP2, nStationPerSide, nRpPerStation = nStationPerSide};
   enum PARTICLE { PION = 0, KAON, PROTON, nParticles };
   enum TPC_TRACK_TYPE { GLO, PRI, TOF, QUA, nTpcTrkTypes }; // GLO=global(all), PRI=primary, TOF=PRI&TofMatched, QUA=TOF&QualityCuts
   enum BUNCH_CROSSING { CB, AG, nBnchXngsTypes }; // CB=colliding bunches, AG=abort gaps
   enum QSUM_2TRKS { OPPO, SAME, nCharges2Trks };
   enum SIGN { PLUS, MINUS, nSigns };
   enum LIST_OF_EFF_CORRECTIONS { RPACC, TPCRECOEFF, TOFMATCHEFF, nEffCorrections };
   enum ANALYSIS_CUT { ALL = 1, TRIG, TWORPTRKS, INFID, ONETOFVX, ZVERTEX, TWOTOFTRKS, ETA, OPPOSITE, EXCLUSIVE, PIPI, KK, PPBAR, nAnalysisCuts };
   enum RANGE_LIMIT { MIN, MAX };
   enum DATASET { MC = 0, MCZB, DATA, nDataSets };
   enum DATATAG { TRUEMC = 0, RECO, nDataTag };

 #endif

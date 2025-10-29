/*
    please check the documentation how to use this class in you analysis -> https://github.com/ladamczy/STAR-Analysis/blob/main/share/HowToUseEfficiencyCorrection.txt
*/

        
#ifndef StEfficiencyCorrector3D_h
#define StEfficiencyCorrector3D_h

#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>
#include <map>
#include <string>

// ROOT
#include "TH3F.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TAxis.h"

class StEfficiencyCorrector3D {
public:
    enum DetectorType { 
        TPC = 0, 
        TOF = 1, 
        N_DETECTORS = 2  
    };
    
    enum ChargeType { POSITIVE = 0, NEGATIVE = 1, N_CHARGES = 2 };
    enum ParticleType { PION = 0, KAON = 1, PROTON = 2, N_PARTICLES = 3 };

    StEfficiencyCorrector3D() : mUseInterpolation(true), mVerbose(false), mMinEfficiency(1e-4) {
        // Initialize all pointers to nullptr
        for(int det = 0; det < N_DETECTORS; det++) {
            for(int charge = 0; charge < N_CHARGES; charge++) {
                for(int particle = 0; particle < N_PARTICLES; particle++) {
                    mEfficiencyHist[det][charge][particle] = nullptr;
                    mTEfficiencyObj[det][charge][particle] = nullptr;
                    mOwnsHistogram[det][charge][particle] = false;
                }
            }
        }
        
        // Initialize scaling parameters to defaults
        for(int charge = 0; charge < N_CHARGES; charge++) {
            for(int particle = 0; particle < N_PARTICLES; particle++) {
                mTpcPlateauShift[charge][particle] = 0.0;
                mTpcMultiplicationFactor[charge][particle] = 1.0;
                mTpcInflationShift[charge][particle] = 0.0;
                mTpcIsUnchanged[charge][particle] = true;
                
                mTofPlateauShift[charge][particle] = 0.0;
                mTofMultiplicationFactor[charge][particle] = 1.0;
                mTofInflationShift[charge][particle] = 0.0;
                mTofIsUnchanged[charge][particle] = true;
            }
        }
    }
    
    ~StEfficiencyCorrector3D() {
        cleanup();
    }
    
    StEfficiencyCorrector3D(const StEfficiencyCorrector3D& other) 
        : mUseInterpolation(other.mUseInterpolation)
        , mVerbose(other.mVerbose)
        , mMinEfficiency(other.mMinEfficiency) {
        
        // Initialize pointers
        for(int det = 0; det < N_DETECTORS; det++) {
            for(int charge = 0; charge < N_CHARGES; charge++) {
                for(int particle = 0; particle < N_PARTICLES; particle++) {
                    mEfficiencyHist[det][charge][particle] = nullptr;
                    mTEfficiencyObj[det][charge][particle] = nullptr;
                    mOwnsHistogram[det][charge][particle] = false;
                }
            }
        }
        
        // Copy efficiency objects
        for(int det = 0; det < N_DETECTORS; det++) {
            for(int charge = 0; charge < N_CHARGES; charge++) {
                for(int particle = 0; particle < N_PARTICLES; particle++) {
                    if(other.mEfficiencyHist[det][charge][particle]) {
                        mEfficiencyHist[det][charge][particle] = (TH3F*)other.mEfficiencyHist[det][charge][particle]->Clone();
                        mOwnsHistogram[det][charge][particle] = true;
                    }
                    if(other.mTEfficiencyObj[det][charge][particle]) {
                        mTEfficiencyObj[det][charge][particle] = other.mTEfficiencyObj[det][charge][particle];
                        mOwnsHistogram[det][charge][particle] = false;
                    }
                }
            }
        }
        
        // Copy scaling parameters
        for(int charge = 0; charge < N_CHARGES; charge++) {
            for(int particle = 0; particle < N_PARTICLES; particle++) {
                mTpcPlateauShift[charge][particle] = other.mTpcPlateauShift[charge][particle];
                mTpcMultiplicationFactor[charge][particle] = other.mTpcMultiplicationFactor[charge][particle];
                mTpcInflationShift[charge][particle] = other.mTpcInflationShift[charge][particle];
                mTpcIsUnchanged[charge][particle] = other.mTpcIsUnchanged[charge][particle];
                
                mTofPlateauShift[charge][particle] = other.mTofPlateauShift[charge][particle];
                mTofMultiplicationFactor[charge][particle] = other.mTofMultiplicationFactor[charge][particle];
                mTofInflationShift[charge][particle] = other.mTofInflationShift[charge][particle];
                mTofIsUnchanged[charge][particle] = other.mTofIsUnchanged[charge][particle];
            }
        }
    }
    
    StEfficiencyCorrector3D& operator=(const StEfficiencyCorrector3D& other) {
        if(this != &other) {
            cleanup();
            
            mUseInterpolation = other.mUseInterpolation;
            mVerbose = other.mVerbose;
            mMinEfficiency = other.mMinEfficiency;
            
            // Copy efficiency objects and scaling parameters (same as copy constructor)
            for(int det = 0; det < N_DETECTORS; det++) {
                for(int charge = 0; charge < N_CHARGES; charge++) {
                    for(int particle = 0; particle < N_PARTICLES; particle++) {
                        if(other.mEfficiencyHist[det][charge][particle]) {
                            mEfficiencyHist[det][charge][particle] = (TH3F*)other.mEfficiencyHist[det][charge][particle]->Clone();
                            mOwnsHistogram[det][charge][particle] = true;
                        }
                        if(other.mTEfficiencyObj[det][charge][particle]) {
                            mTEfficiencyObj[det][charge][particle] = other.mTEfficiencyObj[det][charge][particle];
                            mOwnsHistogram[det][charge][particle] = false;
                        }
                        
                        mTpcPlateauShift[charge][particle] = other.mTpcPlateauShift[charge][particle];
                        mTpcMultiplicationFactor[charge][particle] = other.mTpcMultiplicationFactor[charge][particle];
                        mTpcInflationShift[charge][particle] = other.mTpcInflationShift[charge][particle];
                        mTpcIsUnchanged[charge][particle] = other.mTpcIsUnchanged[charge][particle];
                        
                        mTofPlateauShift[charge][particle] = other.mTofPlateauShift[charge][particle];
                        mTofMultiplicationFactor[charge][particle] = other.mTofMultiplicationFactor[charge][particle];
                        mTofInflationShift[charge][particle] = other.mTofInflationShift[charge][particle];
                        mTofIsUnchanged[charge][particle] = other.mTofIsUnchanged[charge][particle];
                    }
                }
            }
        }
        return *this;
    }

    //=== SETUP METHODS ===
    
    void setTpcEfficiency(TH3F* hist, Int_t charge, Int_t particle, Bool_t takeOwnership = false) {
        setEfficiency(TPC, hist, charge, particle, takeOwnership);
    }
    
    void setTpcEfficiency(TEfficiency* teff, Int_t charge, Int_t particle) {
        setEfficiency(TPC, teff, charge, particle);
    }
    
    void setTofEfficiency(TH3F* hist, Int_t charge, Int_t particle, Bool_t takeOwnership = false) {
        setEfficiency(TOF, hist, charge, particle, takeOwnership);
    }
    
    void setTofEfficiency(TEfficiency* teff, Int_t charge, Int_t particle) {
        setEfficiency(TOF, teff, charge, particle);
    }

    //=== LEGACY SCALING PARAMETER SETTERS ===
    
    /**
     * \brief Set TPC plateau shift (affects efficiency scaling)
     * This is the efficiency offset added to high-pT plateau for scaling calculations
     */
    void setTpcPlateauShift(Int_t charge, Int_t particle, Double_t shift) {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        mTpcPlateauShift[chargeType][particleType] = shift;
    }
    
    /**
     * \brief Set TPC multiplication factor (final efficiency scaling)
     * Applied as: finalEfficiency *= multiplicationFactor
     */
    void setTpcMultiplicationFactor(Int_t charge, Int_t particle, Double_t factor) {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        mTpcMultiplicationFactor[chargeType][particleType] = factor;
    }
    
    /**
     * \brief Set TPC inflation shift (pT clamping adjustment)
     * CLARIFICATION: This shifts pT clamping boundaries, NOT the same as original TF1 inflation.
     * It adjusts where pT values get clamped to histogram edges, not the efficiency fit itself.
     */
    void setTpcInflationShift(Int_t charge, Int_t particle, Double_t shift) {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        mTpcInflationShift[chargeType][particleType] = shift;
    }
    
    /**
     * \brief Set TPC unchanged flag (algorithm selection)
     * true = simple interpolation, false = complex legacy scaling with plateau shift
     */
    void setTpcUnchangedFlag(Int_t charge, Int_t particle, Bool_t unchanged) {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        mTpcIsUnchanged[chargeType][particleType] = unchanged;
    }
    
    /**
     * \brief Set TOF plateau shift (affects efficiency scaling)
     */
    void setTofPlateauShift(Int_t charge, Int_t particle, Double_t shift) {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        mTofPlateauShift[chargeType][particleType] = shift;
    }
    
    /**
     * \brief Set TOF multiplication factor (final efficiency scaling)
     */
    void setTofMultiplicationFactor(Int_t charge, Int_t particle, Double_t factor) {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        mTofMultiplicationFactor[chargeType][particleType] = factor;
    }
    
    /**
     * \brief Set TOF inflation shift (pT clamping adjustment)
     * CLARIFICATION: This shifts pT clamping boundaries, NOT the original TF1 inflation.
     */
    void setTofInflationShift(Int_t charge, Int_t particle, Double_t shift) {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        mTofInflationShift[chargeType][particleType] = shift;
    }
    
    /**
     * \brief Set TOF unchanged flag (algorithm selection)
     */
    void setTofUnchangedFlag(Int_t charge, Int_t particle, Bool_t unchanged) {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        mTofIsUnchanged[chargeType][particleType] = unchanged;
    }

    //=== EFFICIENCY GETTERS ===
    
    Double_t getTpcEfficiency(Double_t eta, Double_t pT, Double_t Vz, Int_t charge, Int_t particle) const {
        if(!validateParameters(eta, pT, Vz, charge, particle)) {
            return mMinEfficiency;
        }
        
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        
        // Try TH3F histogram first
        if(mEfficiencyHist[TPC][chargeType][particleType]) {
            return getTpcEfficiencyLegacy(mEfficiencyHist[TPC][chargeType][particleType], eta, pT, Vz, charge, particle);
        }
        
        // Try TEfficiency object
        if(mTEfficiencyObj[TPC][chargeType][particleType]) {
            return getEfficiencyFromTEfficiency(mTEfficiencyObj[TPC][chargeType][particleType], eta, pT, Vz);
        }
        
        return mMinEfficiency;
    }
    
    Double_t getTofEfficiency(Double_t eta, Double_t pT, Double_t Vz, Int_t charge, Int_t particle) const {
        if(!validateParameters(eta, pT, Vz, charge, particle)) {
            return mMinEfficiency;
        }
        
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        
        // Try TH3F histogram first
        if(mEfficiencyHist[TOF][chargeType][particleType]) {
            return getTofEfficiencyLegacy(mEfficiencyHist[TOF][chargeType][particleType], eta, pT, Vz, charge, particle);
        }
        
        // Try TEfficiency object
        if(mTEfficiencyObj[TOF][chargeType][particleType]) {
            return getEfficiencyFromTEfficiency(mTEfficiencyObj[TOF][chargeType][particleType], eta, pT, Vz);
        }
        
        return mMinEfficiency;
    }
    
    Double_t getCombinedEfficiency(Double_t eta, Double_t pT, Double_t Vz, Int_t charge, Int_t particle) const {
        Double_t tpcEff = getTpcEfficiency(eta, pT, Vz, charge, particle);
        Double_t tofEff = getTofEfficiency(eta, pT, Vz, charge, particle);
        return tpcEff * tofEff;
    }
    
    Double_t getCorrectedYield(Double_t rawYield, Double_t eta, Double_t pT, Double_t Vz,
                               Int_t charge, Int_t particle, Bool_t useTpcOnly = false) const {
        if(rawYield <= 0) return rawYield;
        
        Double_t efficiency = useTpcOnly ? getTpcEfficiency(eta, pT, Vz, charge, particle) 
                                        : getCombinedEfficiency(eta, pT, Vz, charge, particle);
        
        if(efficiency < mMinEfficiency) {
            efficiency = mMinEfficiency;
        }
        
        return rawYield / efficiency;
    }

    //=== UTILITY METHODS ===
    
    Bool_t hasEfficiency(DetectorType detector, Int_t charge, Int_t particle) const {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        return (mEfficiencyHist[detector][chargeType][particleType] != nullptr) || 
               (mTEfficiencyObj[detector][chargeType][particleType] != nullptr);
    }
    
    Bool_t isValid() const {
        for(int det = 0; det < N_DETECTORS; det++) {
            for(int charge = 0; charge < N_CHARGES; charge++) {
                for(int particle = 0; particle < N_PARTICLES; particle++) {
                    if(mEfficiencyHist[det][charge][particle] || mTEfficiencyObj[det][charge][particle]) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
    
    void setVerbose(Bool_t verbose = true) { mVerbose = verbose; }
    void setUseInterpolation(Bool_t useInterpolation = true) { mUseInterpolation = useInterpolation; }
    void setMinEfficiency(Double_t minEff) { mMinEfficiency = minEff; }
    
    void printDebugInfo(Double_t eta, Double_t pT, Double_t Vz, Int_t charge, Int_t particle) const {
        std::cout << "\n=== StEfficiencyCorrector3D Debug Info ===" << std::endl;
        std::cout << "Parameters: eta=" << eta << ", pT=" << pT << ", Vz=" << Vz 
                  << ", charge=" << charge << ", particle=" << particle << std::endl;
        
        if(!validateParameters(eta, pT, Vz, charge, particle)) {
            std::cout << "Invalid parameters!" << std::endl;
            return;
        }
        
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        const char* chargeName = (chargeType == POSITIVE) ? "Positive" : "Negative";
        const char* particleNames[] = {"Pion", "Kaon", "Proton"};
        
        std::cout << "Particle: " << particleNames[particleType] << " (" << chargeName << ")" << std::endl;
        
        // TPC info
        if(hasEfficiency(TPC, charge, particle)) {
            Double_t tpcEff = getTpcEfficiency(eta, pT, Vz, charge, particle);
            std::cout << "TPC efficiency: " << tpcEff << std::endl;
        } else {
            std::cout << "TPC: No efficiency data available" << std::endl;
        }
        
        // TOF info
        if(hasEfficiency(TOF, charge, particle)) {
            Double_t tofEff = getTofEfficiency(eta, pT, Vz, charge, particle);
            std::cout << "TOF efficiency: " << tofEff << std::endl;
        } else {
            std::cout << "TOF: No efficiency data available" << std::endl;
        }
        
        // Combined
        if(hasEfficiency(TPC, charge, particle) && hasEfficiency(TOF, charge, particle)) {
            Double_t combinedEff = getCombinedEfficiency(eta, pT, Vz, charge, particle);
            std::cout << "Combined efficiency: " << combinedEff << std::endl;
            std::cout << "Correction factor: " << (combinedEff > 0 ? 1.0/combinedEff : 0) << std::endl;
        }
        
        std::cout << "========================================\n" << std::endl;
    }
    
    void printHistogramInfo() const {
        std::cout << "\n=== StEfficiencyCorrector3D Histogram Information ===" << std::endl;
        std::cout << "Axis order: X=pT, Y=eta, Z=Vz" << std::endl;
        
        const char* detNames[N_DETECTORS] = {"TPC", "TOF"};
        const char* chargeNames[N_CHARGES] = {"Positive", "Negative"};
        const char* particleNames[N_PARTICLES] = {"Pion", "Kaon", "Proton"};
        
        for(int det = 0; det < N_DETECTORS; det++) {
            for(int charge = 0; charge < N_CHARGES; charge++) {
                for(int particle = 0; particle < N_PARTICLES; particle++) {
                    TH3F* hist = mEfficiencyHist[det][charge][particle];
                    
                    std::cout << "\n" << detNames[det] << " (" << chargeNames[charge] 
                              << " " << particleNames[particle] << "):" << std::endl;
                    
                    if(hist) {
                        std::cout << "  ✓ Available" << std::endl;
                        std::cout << "  Histogram: " << hist->GetName() << std::endl;
                        std::cout << "  pT range:  [" << hist->GetXaxis()->GetXmin() << ", " 
                                  << hist->GetXaxis()->GetXmax() << "] GeV/c" << std::endl;
                        std::cout << "  Eta range: [" << hist->GetYaxis()->GetXmin() << ", " 
                                  << hist->GetYaxis()->GetXmax() << "]" << std::endl;
                        std::cout << "  Vz range:  [" << hist->GetZaxis()->GetXmin() << ", " 
                                  << hist->GetZaxis()->GetXmax() << "] cm" << std::endl;
                    } else {
                        std::cout << "  ✗ Not available" << std::endl;
                    }
                }
            }
        }
        std::cout << "========================================================" << std::endl;
    }

private:
    //=== CORE LEGACY ALGORITHM METHODS ===
    
    Double_t getTpcEfficiencyLegacy(TH3F* hist, Double_t eta, Double_t pT, Double_t Vz,
                                    Int_t charge, Int_t particle) const {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        
        // Apply legacy boundary clamping
        Double_t clampedEta = eta, clampedPt = pT, clampedVz = Vz;
        applyLegacyBoundaryClamping(hist, clampedEta, clampedPt, clampedVz, charge, particle, TPC);
        
        // Check "unchanged" flag
        if(mTpcIsUnchanged[chargeType][particleType]) {
            // Simple interpolation when unchanged
            return mUseInterpolation ? 
                hist->Interpolate(clampedPt, clampedEta, clampedVz) :
                hist->GetBinContent(hist->FindBin(clampedPt, clampedEta, clampedVz));
        } else {
            // Complex scaling when changed
            Double_t highPtEff = calculateHighPtPlateau(hist, clampedEta, clampedVz);
            
            // Apply plateau shift scaling
            Double_t efficiencyScaling = (highPtEff + mTpcPlateauShift[chargeType][particleType]) / 
                                       TMath::Max(highPtEff, 1e-6);
            
            Double_t nominalEff = mUseInterpolation ? 
                hist->Interpolate(clampedPt, clampedEta, clampedVz) :
                hist->GetBinContent(hist->FindBin(clampedPt, clampedEta, clampedVz));
                
            Double_t finalEff = efficiencyScaling * nominalEff;
            return finalEff * mTpcMultiplicationFactor[chargeType][particleType];
        }
    }
    
    Double_t getTofEfficiencyLegacy(TH3F* hist, Double_t eta, Double_t pT, Double_t Vz,
                                    Int_t charge, Int_t particle) const {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        
        // Apply legacy boundary clamping
        Double_t clampedEta = eta, clampedPt = pT, clampedVz = Vz;
        applyLegacyBoundaryClamping(hist, clampedEta, clampedPt, clampedVz, charge, particle, TOF);
        
        // Check "unchanged" flag
        if(mTofIsUnchanged[chargeType][particleType]) {
            // Simple interpolation when unchanged
            return mUseInterpolation ? 
                hist->Interpolate(clampedPt, clampedEta, clampedVz) :
                hist->GetBinContent(hist->FindBin(clampedPt, clampedEta, clampedVz));
        } else {
            // Complex scaling when changed
            Double_t highPtEff = calculateHighPtPlateau(hist, clampedEta, clampedVz);
            
            // Apply plateau shift scaling
            Double_t efficiencyScaling = (highPtEff + mTofPlateauShift[chargeType][particleType]) / 
                                       TMath::Max(highPtEff, 1e-6);
            
            Double_t nominalEff = mUseInterpolation ? 
                hist->Interpolate(clampedPt, clampedEta, clampedVz) :
                hist->GetBinContent(hist->FindBin(clampedPt, clampedEta, clampedVz));
                
            Double_t finalEff = efficiencyScaling * nominalEff;
            return finalEff * mTofMultiplicationFactor[chargeType][particleType];
        }
    }
    
    /**
     * \brief Legacy boundary clamping with 1e-3 offsets - adapted for YOUR axis order
     * YOUR histograms: X=pT, Y=eta, Z=Vz
     * Legacy logic preserved: clamp to bin centers ± 1e-3
     * 
     * CLARIFICATION: The "inflation shift" parameter here is a pT clamp adjustment,
     * NOT identical to the original TF1 inflation parameter from legacy code.
     * It shifts the pT clamping boundaries but doesn't modify the TF1 fit itself.
     */
    void applyLegacyBoundaryClamping(TH3F* hist, Double_t& eta, Double_t& pT, Double_t& Vz,
                                     Int_t charge, Int_t particle, DetectorType detector) const {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        
        // Legacy clamping with 1e-3 offsets - YOUR axis order: X=pT, Y=eta, Z=Vz
        
        // Clamp pT (X-axis in YOUR histograms)
        if(hist->GetXaxis()->GetBinCenter(1) > pT) {
            pT = hist->GetXaxis()->GetBinCenter(1) + 1e-3;
        } else if(hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) < pT) {
            pT = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) - 1e-3;
        }
        
        // Apply inflation shift to pT clamping boundaries (pT clamp adjustment)
        // NOTE: This is NOT the same as the original TF1 inflation parameter!
        // It only shifts where we clamp pT values, not the efficiency calculation itself.
        Double_t inflationShift = (detector == TPC) ? 
            mTpcInflationShift[chargeType][particleType] : 
            mTofInflationShift[chargeType][particleType];
        
        if(inflationShift != 0.0) {
            Double_t adjustedPt = pT + inflationShift;
            if(hist->GetXaxis()->GetBinCenter(1) > adjustedPt) {
                pT = hist->GetXaxis()->GetBinCenter(1) + 1e-3 - inflationShift;
            } else if(hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) < adjustedPt) {
                pT = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) - 1e-3 - inflationShift;
            }
        }
        
        // Clamp eta (Y-axis in YOUR histograms)  
        if(hist->GetYaxis()->GetBinCenter(1) > eta) {
            eta = hist->GetYaxis()->GetBinCenter(1) + 1e-3;
        } else if(hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()) < eta) {
            eta = hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()) - 1e-3;
        }
        
        // Clamp Vz (Z-axis in YOUR histograms)
        if(hist->GetZaxis()->GetBinCenter(1) > Vz) {
            Vz = hist->GetZaxis()->GetBinCenter(1) + 1e-3;
        } else if(hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ()) < Vz) {
            Vz = hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ()) - 1e-3;
        }
    }
    
    /**
     * \brief High-pT plateau averaging (legacy method) with robustness fallback
     * YOUR histograms: X=pT, so we average last 3 X-axis bins (not Y-axis)
     * 
     * Includes fallback to bin content if interpolation returns 0/NaN at edges
     */
    Double_t calculateHighPtPlateau(TH3F* hist, Double_t clampedEta, Double_t clampedVz) const {
        // Legacy high-pT averaging of last 3 bins - YOUR axis order: X=pT
        Double_t eff1 = hist->Interpolate(hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) - 1e-3, clampedEta, clampedVz);
        Double_t eff2 = hist->Interpolate(hist->GetXaxis()->GetBinCenter(hist->GetNbinsX() - 1), clampedEta, clampedVz);
        Double_t eff3 = hist->Interpolate(hist->GetXaxis()->GetBinCenter(hist->GetNbinsX() - 2), clampedEta, clampedVz);
        
        // Robustness: fallback to bin content if interpolation fails at edges
        if(!TMath::Finite(eff1) || eff1 <= 0) {
            eff1 = hist->GetBinContent(hist->GetNbinsX(), 
                                      hist->GetYaxis()->FindBin(clampedEta), 
                                      hist->GetZaxis()->FindBin(clampedVz));
        }
        if(!TMath::Finite(eff2) || eff2 <= 0) {
            eff2 = hist->GetBinContent(hist->GetNbinsX() - 1, 
                                      hist->GetYaxis()->FindBin(clampedEta), 
                                      hist->GetZaxis()->FindBin(clampedVz));
        }
        if(!TMath::Finite(eff3) || eff3 <= 0) {
            eff3 = hist->GetBinContent(hist->GetNbinsX() - 2, 
                                      hist->GetYaxis()->FindBin(clampedEta), 
                                      hist->GetZaxis()->FindBin(clampedVz));
        }
        
        Double_t highPtEff = (eff1 + eff2 + eff3) / 3.0;
        
        // Additional robustness: if all three are bad, try single best bin
        if(!TMath::Finite(highPtEff) || highPtEff <= 0) {
            // Try the second-to-last bin (often more stable than the very last)
            highPtEff = hist->GetBinContent(hist->GetNbinsX() - 1, 
                                          hist->GetYaxis()->FindBin(clampedEta), 
                                          hist->GetZaxis()->FindBin(clampedVz));
            
            if(!TMath::Finite(highPtEff) || highPtEff <= 0) {
                // Last resort: use minimum efficiency
                highPtEff = mMinEfficiency * 10; // Give it some headroom for scaling
            }
        }
        
        return TMath::Max(mMinEfficiency, highPtEff);
    }
Double_t getEfficiencyFromTEfficiency(TEfficiency* teff, Double_t eta, Double_t pT, Double_t Vz) const {
    if(!teff) return mMinEfficiency;
    
    TH3F* totalHist = (TH3F*)teff->GetTotalHistogram();
    if(!totalHist) return mMinEfficiency;
    
    // IMPROVED: Clamp to BIN CENTERS, not boundaries, with small margins
    // This prevents interpolation failures at histogram edges
    
    // Get axis information
    TAxis* xAxis = totalHist->GetXaxis();  // pT
    TAxis* yAxis = totalHist->GetYaxis();  // eta  
    TAxis* zAxis = totalHist->GetZaxis();  // Vz
    
    // Clamp to bin centers with small margins (like legacy algorithm)
    Double_t clampedPt = pT;
    if(xAxis->GetBinCenter(1) > pT) {
        clampedPt = xAxis->GetBinCenter(1) + 1e-3;
    } else if(xAxis->GetBinCenter(xAxis->GetNbins()) < pT) {
        clampedPt = xAxis->GetBinCenter(xAxis->GetNbins()) - 1e-3;
    }
    
    Double_t clampedEta = eta;
    if(yAxis->GetBinCenter(1) > eta) {
        clampedEta = yAxis->GetBinCenter(1) + 1e-3;
    } else if(yAxis->GetBinCenter(yAxis->GetNbins()) < eta) {
        clampedEta = yAxis->GetBinCenter(yAxis->GetNbins()) - 1e-3;
    }
    
    Double_t clampedVz = Vz;
    if(zAxis->GetBinCenter(1) > Vz) {
        clampedVz = zAxis->GetBinCenter(1) + 1e-3;
    } else if(zAxis->GetBinCenter(zAxis->GetNbins()) < Vz) {
        clampedVz = zAxis->GetBinCenter(zAxis->GetNbins()) - 1e-3;
    }
    
    Double_t efficiency = 0.0;
    
    // Method 1: Try interpolation with properly clamped coordinates
    TH3F* passedHist = (TH3F*)teff->GetPassedHistogram();
    if(passedHist && totalHist) {
        // Check if coordinates are within interpolation-safe bounds
        Bool_t safeToInterpolate = (clampedPt > totalHist->GetXaxis()->GetXmin() + 1e-6 &&
                                   clampedPt < totalHist->GetXaxis()->GetXmax() - 1e-6 &&
                                   clampedEta > totalHist->GetYaxis()->GetXmin() + 1e-6 &&
                                   clampedEta < totalHist->GetYaxis()->GetXmax() - 1e-6 &&
                                   clampedVz > totalHist->GetZaxis()->GetXmin() + 1e-6 &&
                                   clampedVz < totalHist->GetZaxis()->GetXmax() - 1e-6);
        
        if(safeToInterpolate) {
            Double_t passed = passedHist->Interpolate(clampedPt, clampedEta, clampedVz);
            Double_t total = totalHist->Interpolate(clampedPt, clampedEta, clampedVz);
            
            if(total > 0) {
                efficiency = passed / total;
            }
        }
    }
    
    // Method 2: Fallback to bin content if interpolation failed or unsafe
    if(efficiency <= mMinEfficiency && passedHist && totalHist) {
        Int_t binX = totalHist->GetXaxis()->FindBin(clampedPt);
        Int_t binY = totalHist->GetYaxis()->FindBin(clampedEta);
        Int_t binZ = totalHist->GetZaxis()->FindBin(clampedVz);
        
        // Ensure bins are within valid range
        binX = TMath::Max(1, TMath::Min(totalHist->GetNbinsX(), binX));
        binY = TMath::Max(1, TMath::Min(totalHist->GetNbinsY(), binY));
        binZ = TMath::Max(1, TMath::Min(totalHist->GetNbinsZ(), binZ));
        
        Double_t passed = passedHist->GetBinContent(binX, binY, binZ);
        Double_t total = totalHist->GetBinContent(binX, binY, binZ);
        
        if(total > 0) {
            efficiency = passed / total;
        }
    }
    
    // Don't spam verbose output for expected edge cases
    if(mVerbose && efficiency <= mMinEfficiency) {
        // Only print verbose if this wasn't an extreme coordinate test
        Bool_t isExtremeCoordinate = (TMath::Abs(eta) > 5.0 || pT < 0.05 || pT > 10.0 || TMath::Abs(Vz) > 100.0);
        if(!isExtremeCoordinate) {
            std::cout << "TEfficiency lookup failed for (" << clampedPt << "," << clampedEta << "," << clampedVz 
                      << ") - using minimum efficiency" << std::endl;
        }
    }
    
    return TMath::Max(mMinEfficiency, efficiency);
}
    void setEfficiency(DetectorType detector, TH3F* hist, Int_t charge, Int_t particle, Bool_t takeOwnership) {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        
        // Clean up existing data
        if(mOwnsHistogram[detector][chargeType][particleType] && mEfficiencyHist[detector][chargeType][particleType]) {
            delete mEfficiencyHist[detector][chargeType][particleType];
        }
        
        mEfficiencyHist[detector][chargeType][particleType] = hist;
        mTEfficiencyObj[detector][chargeType][particleType] = nullptr;
        mOwnsHistogram[detector][chargeType][particleType] = takeOwnership;
    }
    
    void setEfficiency(DetectorType detector, TEfficiency* teff, Int_t charge, Int_t particle) {
        ChargeType chargeType = getChargeType(charge);
        ParticleType particleType = getParticleType(particle);
        
        // Clean up existing data
        if(mOwnsHistogram[detector][chargeType][particleType] && mEfficiencyHist[detector][chargeType][particleType]) {
            delete mEfficiencyHist[detector][chargeType][particleType];
        }
        
        mEfficiencyHist[detector][chargeType][particleType] = nullptr;
        mTEfficiencyObj[detector][chargeType][particleType] = teff;
        mOwnsHistogram[detector][chargeType][particleType] = false; // Never own TEfficiency objects
    }
    
    void cleanup() {
        for(int det = 0; det < N_DETECTORS; det++) {
            for(int charge = 0; charge < N_CHARGES; charge++) {
                for(int particle = 0; particle < N_PARTICLES; particle++) {
                    if(mOwnsHistogram[det][charge][particle] && mEfficiencyHist[det][charge][particle]) {
                        delete mEfficiencyHist[det][charge][particle];
                    }
                    mEfficiencyHist[det][charge][particle] = nullptr;
                    mTEfficiencyObj[det][charge][particle] = nullptr;
                    mOwnsHistogram[det][charge][particle] = false;
                }
            }
        }
    }
    
    ChargeType getChargeType(Int_t charge) const { 
        return (charge > 0) ? POSITIVE : NEGATIVE; 
    }
    
    ParticleType getParticleType(Int_t particle) const {
        if(particle == 0) return PION;
        if(particle == 1) return KAON;
        if(particle == 2) return PROTON;
        return PROTON; // Default fallback
    }
    
    Bool_t validateParameters(Double_t eta, Double_t pT, Double_t Vz, Int_t charge, Int_t particle) const {
        if(pT <= 0) {
            if(mVerbose) std::cerr << "StEfficiencyCorrector3D::ERROR - Invalid pT: " << pT << std::endl;
            return false;
        }
        if(TMath::Abs(eta) > 10) {
            if(mVerbose) std::cerr << "StEfficiencyCorrector3D::ERROR - Invalid eta: " << eta << std::endl;
            return false;
        }
        //if(TMath::Abs(Vz) > 200) {
        //    if(mVerbose) std::cerr << "StEfficiencyCorrector3D::ERROR - Invalid Vz: " << Vz << std::endl;
        //    return false;
        //}
        if(charge != 1 && charge != -1) {
            if(mVerbose) std::cerr << "StEfficiencyCorrector3D::ERROR - Invalid charge: " << charge << std::endl;
            return false;
        }
        if(particle < 0 || particle > 2) {
            if(mVerbose) std::cerr << "StEfficiencyCorrector3D::ERROR - Invalid particle: " << particle << std::endl;
            return false;
        }
        return true;
    }

private:
    // Data storage (3D arrays for charge/particle combinations)
    TH3F*        mEfficiencyHist[N_DETECTORS][N_CHARGES][N_PARTICLES];
    TEfficiency* mTEfficiencyObj[N_DETECTORS][N_CHARGES][N_PARTICLES];
    Bool_t       mOwnsHistogram[N_DETECTORS][N_CHARGES][N_PARTICLES];
    
    // Legacy scaling parameters
    Double_t mTpcPlateauShift[N_CHARGES][N_PARTICLES];
    Double_t mTpcMultiplicationFactor[N_CHARGES][N_PARTICLES];
    Double_t mTpcInflationShift[N_CHARGES][N_PARTICLES];
    Bool_t   mTpcIsUnchanged[N_CHARGES][N_PARTICLES];
    
    Double_t mTofPlateauShift[N_CHARGES][N_PARTICLES];
    Double_t mTofMultiplicationFactor[N_CHARGES][N_PARTICLES];
    Double_t mTofInflationShift[N_CHARGES][N_PARTICLES];
    Bool_t   mTofIsUnchanged[N_CHARGES][N_PARTICLES];
    
    // Configuration
    Bool_t   mUseInterpolation;
    Bool_t   mVerbose;
    Double_t mMinEfficiency;
};

#endif // StEfficiencyCorrector3D_h

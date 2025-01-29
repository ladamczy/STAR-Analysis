#include "ExclusiveCode.h"


void CheckN1(vector <bool *> vCuts,bool & condition, bool & conditionN1 )
{
    for (int i = 0; i < vCuts.size(); i++ )
    {
     
        if (vCuts[i]!=&condition and *(vCuts[i])==false )
        {
            conditionN1 = false;
            break;
        }
    }
}


void FindProtons(bool isMC, StRPEvent *rpEvt, StUPCEvent *upcEvt, TLorentzVector & protonE, TLorentzVector & protonW)
{
   
    if (isMC == 0 && rpEvt)  // Add null check)
    {
            if (rpEvt->getNumberOfTracks() < 2) return;  // Add safety check
            
            StUPCRpsTrack *trk = rpEvt->getTrack(0);
            StUPCRpsTrack *trk2 = rpEvt->getTrack(1);
            
            trk->setEvent(rpEvt);
            trk2->setEvent(rpEvt);
            if (!trk || !trk2) return;  // Add safety check

            if (trk->branch() < 2)
            {
                protonE.SetPxPyPzE(ExclusiveK0K0::BEAM_ENERGY*trk->thetaRp(0), ExclusiveK0K0::BEAM_ENERGY*trk->thetaRp(1), -(ExclusiveK0K0::BEAM_ENERGY-ExclusiveK0K0::BEAM_ENERGY*trk->thetaRp(1)-ExclusiveK0K0::BEAM_ENERGY*trk->thetaRp(0)), ExclusiveK0K0::BEAM_ENERGY);
                protonW.SetPxPyPzE(ExclusiveK0K0::BEAM_ENERGY*trk2->thetaRp(0), ExclusiveK0K0::BEAM_ENERGY*trk2->thetaRp(1),ExclusiveK0K0::BEAM_ENERGY-ExclusiveK0K0::BEAM_ENERGY*trk2->thetaRp(1)-ExclusiveK0K0::BEAM_ENERGY*trk2->thetaRp(0), ExclusiveK0K0::BEAM_ENERGY);    
            }

            else if (trk2->branch() < 2)
            {
                protonE.SetPxPyPzE(ExclusiveK0K0::BEAM_ENERGY*trk2->thetaRp(0), ExclusiveK0K0::BEAM_ENERGY*trk2->thetaRp(1), -(ExclusiveK0K0::BEAM_ENERGY-ExclusiveK0K0::BEAM_ENERGY*trk2->thetaRp(1)-ExclusiveK0K0::BEAM_ENERGY*trk2->thetaRp(0)), ExclusiveK0K0::BEAM_ENERGY);
                protonW.SetPxPyPzE(ExclusiveK0K0::BEAM_ENERGY*trk->thetaRp(0), ExclusiveK0K0::BEAM_ENERGY*trk->thetaRp(1), ExclusiveK0K0::BEAM_ENERGY-ExclusiveK0K0::BEAM_ENERGY*trk->thetaRp(1)-ExclusiveK0K0::BEAM_ENERGY*trk->thetaRp(0), ExclusiveK0K0::BEAM_ENERGY);    
            }
    }    

    else if (isMC == 1)
    {
        vector <TParticle*> vProtons;
        TParticle* particle;
        
        for (int i = 0; i < upcEvt->getNumberOfMCParticles(); i++)
        {
            particle = upcEvt->getMCParticle(i);
            if (particle->GetPDG()->PdgCode() == 2212 and particle->GetFirstMother() == 1)
            {
                vProtons.push_back(particle);
            }
        }

        double sigma = 0.033;
        double sigmaz = 0.51;
        double px_reco1 = gRandom->Gaus(vProtons[0]->Px(), sigma);
        double py_reco1 = gRandom->Gaus(vProtons[0]->Py(), sigma);
        double pz_reco1 = gRandom->Gaus(vProtons[0]->Pz(), sigmaz);
        protonE.SetXYZM(px_reco1,py_reco1, pz_reco1, ExclusiveK0K0::MASS_PROTON);

        double px_reco2 = gRandom->Gaus(vProtons[1]->Px(), sigma);
        double py_reco2 = gRandom->Gaus(vProtons[1]->Py(), sigma);
        double pz_reco2 = gRandom->Gaus(vProtons[1]->Pz(), sigmaz);
        protonW.SetXYZM(px_reco2,py_reco2, pz_reco2, ExclusiveK0K0::MASS_PROTON);
        
        if (pz_reco1 < 0.0)
        {
            protonE.SetXYZM(px_reco1,py_reco1, pz_reco1, ExclusiveK0K0::MASS_PROTON);
            protonW.SetXYZM(px_reco2,py_reco2, pz_reco2, ExclusiveK0K0::MASS_PROTON);
        }

        else
        {
            protonW.SetXYZM(px_reco1,py_reco1, pz_reco1, ExclusiveK0K0::MASS_PROTON);
            protonE.SetXYZM(px_reco2,py_reco2, pz_reco2, ExclusiveK0K0::MASS_PROTON);
        }

        vProtons.clear();
    }

   
}



void SeparateTracks(StUPCEvent *upcEvt, vector <StUPCTrack const*> &tracksWithTofHit, vector <StUPCTrack const*> &tracksWithoutTofHit, bool isMC, TH1D* HistNumWithTofTrakcs,TH1D* HistNumWithoutTofTrakcs)
{

    for (Int_t i = 0; i<upcEvt->getNumberOfTracks(); i++)
    {

        if (isMC == 0)
        {
            if((upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof)) and (upcEvt->getTrack(i)->getFlag(StUPCTrack::kV0)))  
            {
                tracksWithTofHit.push_back(upcEvt->getTrack(i));
            } 
            else
            {
                tracksWithoutTofHit.push_back(upcEvt->getTrack(i));
            }
        }

        else if (isMC == 1)
        {
            if((upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof)))  
            {
                tracksWithTofHit.push_back(upcEvt->getTrack(i));
            } 
            else
            {
                tracksWithoutTofHit.push_back(upcEvt->getTrack(i));
            } 
        }
    }
        if (isMC == 1)
        {

            
            Int_t BBCLEast = upcEvt->getBBCLargeEast();
            Int_t BBCLWest = upcEvt->getBBCLargeWest();

            if (tracksWithoutTofHit.size() >= 2 and BBCLWest <= 52 and BBCLEast <= 52.0)
            {
                HistNumWithTofTrakcs->Fill(tracksWithTofHit.size());
                HistNumWithoutTofTrakcs->Fill(tracksWithoutTofHit.size());
            }

        }
        else
        {
                    HistNumWithTofTrakcs->Fill(tracksWithTofHit.size());
        HistNumWithoutTofTrakcs->Fill(tracksWithoutTofHit.size());

        }
}



bool ValidNumberOfTofTracks(vector <int> vNums, vector <StUPCTrack const*>  tracksWithTofHit , vector <StUPCTrack const*>  tracksWithoutTofHit,vector<TH1D *>&HistPtPionWithTof,vector<TH1D *>&HistPtPionWithoutTof, vector<TH1D *>& HistEtaPionWithTof, vector<TH1D *> &HistEtaPionWithoutTof, vector<TH1D *>& HistNfitPionWithTof,vector<TH1D *>  &HistNfitPionWithoutTof  )
{
    
    bool isValidNumberOfTofMatchedTracks = false;
    for (int i = 0; i < vNums.size(); i++)
    {

        if (vNums[i]!=6)
        {
            if (tracksWithTofHit.size() == vNums[i])
            {
                isValidNumberOfTofMatchedTracks = true;
            }
        }
        else
        {
            if (tracksWithTofHit.size() >= vNums[i])
            {
                isValidNumberOfTofMatchedTracks = true;
            }
        }     
    }


     if (isValidNumberOfTofMatchedTracks == true)
    {
      
        int histInd = tracksWithTofHit.size()-2;
        for (int i = 0; i < tracksWithTofHit.size(); i++)
        {            
            HistPtPionWithTof[histInd]->Fill( tracksWithTofHit[i]->getPt());
            HistEtaPionWithTof[histInd]->Fill( tracksWithTofHit[i]->getEta());
            HistNfitPionWithTof[histInd]->Fill( tracksWithTofHit[i]->getNhitsFit());           
        }

        for (int i = 0; i < tracksWithoutTofHit.size(); i++)
        {            
            HistPtPionWithoutTof[histInd]->Fill( tracksWithoutTofHit[i]->getPt());
            HistEtaPionWithoutTof[histInd]->Fill( tracksWithoutTofHit[i]->getEta());
            HistNfitPionWithoutTof[histInd]->Fill( tracksWithoutTofHit[i]->getNhitsFit());           
        }
    }
    return isValidNumberOfTofMatchedTracks;
}


bool AreTofTracksGood(vector <StUPCTrack const*>  &tracksWithTofHit)
{
    bool flag = true;
    int num6 = 0;
    int num5 = 0;
    int num4 = 0;




    if (tracksWithTofHit.size() == 5)
    {    
        vector <int> index;
        for (int i = 0; i < tracksWithTofHit.size(); i++)
        {   
            StUPCTrack const* tempTrack = tracksWithTofHit[i];
            if (abs(tempTrack->getEta()) <= ExclusiveK0K0::PION_ETA and tempTrack->getPt() >= ExclusiveK0K0::PION_PT and tempTrack->getNhitsFit() >= ExclusiveK0K0::N_FIT)
            {
                num5+=1;
            }

            else
            {
                index.push_back(i);
            }
        }
 
        if (index.size() != 0)
        {
            for (int i = index.size()-1; i >=0; i--)
            {
                int ind = index[i];
                tracksWithTofHit.erase(tracksWithTofHit.begin()+ind);
            }  
        }
        
        if (num5 < 3)
        {       
            flag = false;
        }
        return flag;
    }

    if (tracksWithTofHit.size() == 4)
    {    
        vector <int> index;
        for (int i = 0; i < tracksWithTofHit.size(); i++)
        {   
            StUPCTrack const* tempTrack = tracksWithTofHit[i];
            if (abs(tempTrack->getEta()) <= ExclusiveK0K0::PION_ETA and tempTrack->getPt() >= ExclusiveK0K0::PION_PT and tempTrack->getNhitsFit() >= ExclusiveK0K0::N_FIT)
            {
                num4+=1;
            }

            else
            {
                index.push_back(i);
            }
        }
 
        if (index.size() != 0)
        {
            for (int i = index.size()-1; i >=0; i--)
            {
                int ind = index[i];
                tracksWithTofHit.erase(tracksWithTofHit.begin()+ind);
            }  
        }

        if (num4 < 3)
        {       
            flag = false;
        }
        return flag;
    }

    else  if (tracksWithTofHit.size() == 3 or tracksWithTofHit.size() == 2 )
    {
        for (int i = 0; i < tracksWithTofHit.size(); i++)
        {
            StUPCTrack const* tempTrack = tracksWithTofHit[i];  
            if (abs(tempTrack->getEta()) > ExclusiveK0K0::PION_ETA or tempTrack->getPt() < ExclusiveK0K0::PION_PT or tempTrack->getNhitsFit() < ExclusiveK0K0::N_FIT)
            {
                flag = false;
            }
        }
        return flag;
    }

    else
    {    
        vector <int> index;
        for (int i = 0; i < tracksWithTofHit.size(); i++)
        {   
            StUPCTrack const* tempTrack = tracksWithTofHit[i];
            if (abs(tempTrack->getEta()) <= ExclusiveK0K0::PION_ETA and tempTrack->getPt() >= ExclusiveK0K0::PION_PT and tempTrack->getNhitsFit() >= ExclusiveK0K0::N_FIT)
            {
                num6+=1;
            }

            else
            {
                index.push_back(i);
            }
        }
 
        if (index.size() != 0)
        {
            for (int i = index.size()-1; i >=0; i--)
            {
                int ind = index[i];
                tracksWithTofHit.erase(tracksWithTofHit.begin()+ind);
            }  
        }
        
    
        if (num6 < 3)
        {       
            flag = false;
        }

        if (num6>5)
        {
            flag = false;
        }
        return flag;
    }
}



bool FindTracks( vector <StUPCTrack const*> &tracksWithTofHit, vector <StUPCTrack const*> &tracksWithoutTofHit, vector <StUPCTrack const*> &goodTracksWithoutTofHit)
{
    bool flag = false;
    int sumCharge = 0;

    if (tracksWithTofHit.size() == 5)
    {
        for (int i =0; i < tracksWithTofHit.size(); i++)
        {
            sumCharge+=(tracksWithTofHit[i]->getCharge());
        }

        if (sumCharge == 1 or sumCharge == -1)
        {
            flag = true;
        }
    }

    if (tracksWithTofHit.size() == 4)
    {	
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {
           
            sumCharge+=(tracksWithTofHit[i]->getCharge());
        }

        if (sumCharge == 0)  
        {
            flag = true;
        }
    }

    else if (tracksWithTofHit.size() == 3)
    {
        sumCharge = tracksWithTofHit[0]->getCharge() + tracksWithTofHit[1]->getCharge() + tracksWithTofHit[2]->getCharge();    
       
        for (int i = 0; i < tracksWithoutTofHit.size(); i++)
        {
            if (tracksWithoutTofHit[i]->getCharge() == sumCharge*(-1) and abs(tracksWithoutTofHit[i]->getEta()) <= ExclusiveK0K0::PION_ETA and  tracksWithoutTofHit[i]->getPt() >= ExclusiveK0K0::PION_PT and tracksWithoutTofHit[i]->getNhitsFit() >= ExclusiveK0K0::N_FIT)// and tracksWithoutTofHit[i]->getNhitsDEdx() >= 15)
            {
                goodTracksWithoutTofHit.push_back(tracksWithoutTofHit[i]);     
            }
        }
        
        if (goodTracksWithoutTofHit.size() != 0) 
        {
            flag = true;
        }
    }


    else if (tracksWithTofHit.size() == 2)
    {
        sumCharge = tracksWithTofHit[0]->getCharge() + tracksWithTofHit[1]->getCharge();   
        if (sumCharge == 2 or sumCharge == -2)
        {
            for (int i = 0; i < tracksWithoutTofHit.size(); i++)
            {
                if (tracksWithoutTofHit[i]->getCharge() == int(sumCharge/(-2)) and abs(tracksWithoutTofHit[i]->getEta()) <= ExclusiveK0K0::PION_ETA and  tracksWithoutTofHit[i]->getPt() >= ExclusiveK0K0::PION_PT  and tracksWithoutTofHit[i]->getNhitsFit() >= ExclusiveK0K0::N_FIT)// and tracksWithoutTofHit[i]->getNhitsDEdx() >= 15)
                {
                    goodTracksWithoutTofHit.push_back(tracksWithoutTofHit[i]);     
                }
            }

            if (goodTracksWithoutTofHit.size() > 1) 
            {
                flag = true;
            }
        }

        else if (sumCharge == 0)
        {
            int sumFlag1 = 0;
            int sumFlag2 = 0;
            for (int i = 0; i < tracksWithoutTofHit.size(); i++)
            {
                if (tracksWithoutTofHit[i]->getCharge() == 1  and abs(tracksWithoutTofHit[i]->getEta()) <= ExclusiveK0K0::PION_ETA and  tracksWithoutTofHit[i]->getPt() >= ExclusiveK0K0::PION_PT  and tracksWithoutTofHit[i]->getNhitsFit() >= ExclusiveK0K0::N_FIT)
                {
                    goodTracksWithoutTofHit.push_back(tracksWithoutTofHit[i]);    
                    sumFlag1+=1;
                }
                else if (tracksWithoutTofHit[i]->getCharge() == -1 and  abs(tracksWithoutTofHit[i]->getEta()) <= ExclusiveK0K0::PION_ETA and tracksWithoutTofHit[i]->getPt() >= ExclusiveK0K0::PION_PT  and tracksWithoutTofHit[i]->getNhitsFit() >= ExclusiveK0K0::N_FIT)
                {
                    goodTracksWithoutTofHit.push_back(tracksWithoutTofHit[i]);     
                    sumFlag2+=1;
                }
            }
            if (sumFlag1!=0 and sumFlag2!=0)
            {
                flag = true;
            }
        }
    }
    return flag;
}




void FillProtons(double protonSignsY, TH1D*HistProtonsSameSide, TH1D*HistProtonsOppositeSide, int iCutFlow)
{
    if (protonSignsY >= 0.0)
    {
        HistProtonsSameSide->Fill(iCutFlow);
    }
    else 
    {
       HistProtonsOppositeSide->Fill(iCutFlow);      
    }
}

void GetBeamPar(StUPCEvent *upcEvt, double * beamPar, bool isMC)
{

 
    if (isMC == 0)
    {    
        beamPar[0] = upcEvt->getBeamXPosition();
        beamPar[1] = upcEvt->getBeamXSlope();
        beamPar[2] = upcEvt->getBeamYPosition();
        beamPar[3] = upcEvt->getBeamYSlope();
    }

    else
    {
        beamPar[0] = 0.0;
        beamPar[1] = 0.0;
        beamPar[2] = 0.0;
        beamPar[3] = 0.0;
    }

}

bool FindProtonsAndPions(StUPCEvent *upcEvt, vector <StUPCTrack const*> &tracksWithTofHit,  vector <StUPCTrack const*> &goodTracksWithoutTofHit,vector <StUPCTrack const*> &  vFinalProtons, vector <StUPCTrack const*> & vFinalPions, double * beamPar)
{
    bool flag = true;

  
    if (tracksWithTofHit.size() == 4)
    {
       
        vector <StUPCTrack const *> vProtons2Sigma;
        vector <StUPCTrack const *> vPions3Sigma;
        for (int i = 0; i < tracksWithTofHit.size(); i++)
        {

            StUPCTrack const * track = tracksWithTofHit[i];
            float sigmaProton = track->getNSigmasTPCProton();
            float sigmaPion = track->getNSigmasTPCPion();
            if (abs(sigmaProton) < 2.0)
            {
                vProtons2Sigma.push_back(track);
            }

            else if (abs(sigmaPion) < 3.0)
            {
                 vPions3Sigma.push_back(track);
            }
        }
   
        if (vProtons2Sigma.size() < 2) {flag = false; return flag;}

        else if (vProtons2Sigma.size() == 2)
        {
            if ( (vProtons2Sigma[0]->getCharge()+ vProtons2Sigma[1]->getCharge()) != 0 ) {flag = false; return flag;}

            if (vPions3Sigma.size() < 2 ) {flag = false; return flag;}
            else if ( (vPions3Sigma[0]->getCharge() +vPions3Sigma[1]->getCharge())!=0 ){flag = false; return flag;}

            vFinalPions.push_back(vPions3Sigma[0]);vFinalPions.push_back(vPions3Sigma[1]);
            vFinalProtons.push_back(vProtons2Sigma[0]);vFinalProtons.push_back(vProtons2Sigma[1]);
        }

        else if (vProtons2Sigma.size() == 3)
        {   
            if (vPions3Sigma.size() == 0) {flag = false; return flag;}
    
            double pionCharge = vPions3Sigma[0]->getCharge();
            vector <StUPCTrack const *> protonsPionSign;
            vector <StUPCTrack const *> protonsNoPionSign;
            
            for (int i = 0; i < vProtons2Sigma.size(); i++)
            {
                if (vProtons2Sigma[i]->getCharge() == pionCharge )
                {
                    protonsPionSign.push_back(vProtons2Sigma[i]);
                }
                else
                {
                    protonsNoPionSign.push_back(vProtons2Sigma[i]);
                }
            }

            vector <StUPCTrack const *> protonsNoPionSignThatCanBeAlsoPion;
            for (int i = 0; i < protonsNoPionSign.size(); i++)
            {
                float sigmaPion = protonsNoPionSign[i]->getNSigmasTPCPion();
                if (abs(sigmaPion) < 3.0)
                {
                    protonsNoPionSignThatCanBeAlsoPion.push_back(protonsNoPionSign[i]);
                }
            }

            if (protonsNoPionSignThatCanBeAlsoPion.size() == 0) {flag = false; return flag;}
            else if (protonsNoPionSignThatCanBeAlsoPion.size() == 1)
            {
                if (protonsNoPionSign[0]->getPt() == protonsNoPionSignThatCanBeAlsoPion[0]->getPt())
                {
                    vFinalPions.push_back(vPions3Sigma[0]);  vFinalPions.push_back(protonsNoPionSign[0]);
             
                    vFinalProtons.push_back(protonsPionSign[0]);  vFinalProtons.push_back(protonsNoPionSign[1]);
                    
                }
                else
                {
                    vFinalPions.push_back(vPions3Sigma[0]);  vFinalPions.push_back(protonsNoPionSign[1]);
                    vFinalProtons.push_back(protonsPionSign[0]);  vFinalProtons.push_back(protonsNoPionSign[0]);
                }   
            }
        
            else if (protonsNoPionSignThatCanBeAlsoPion.size() == 2)
            {
                TVector3 const tryVec(0,0,0);
                StUPCV0 lVecKaonComb1a(vPions3Sigma[0], protonsNoPionSign[1], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                StUPCV0 lVecKaonComb1b(protonsNoPionSign[0], protonsPionSign[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);

                StUPCV0 lVecKaonComb2a(vPions3Sigma[0], protonsNoPionSign[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                StUPCV0 lVecKaonComb2b(protonsNoPionSign[1], protonsPionSign[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
            
                double  distSumSquaredComb1 = sqrt(pow((lVecKaonComb1a.m()-ExclusiveK0K0::MASS_LAMBDA), 2) + pow((lVecKaonComb1b.m()-ExclusiveK0K0::MASS_LAMBDA), 2));
                double  distSumSquaredComb2 = sqrt(pow((lVecKaonComb2a.m()-ExclusiveK0K0::MASS_LAMBDA), 2) + pow((lVecKaonComb2b.m()-ExclusiveK0K0::MASS_LAMBDA), 2));

                if (distSumSquaredComb1 < distSumSquaredComb2)
                {
                    vFinalPions.push_back(vPions3Sigma[0]);
                    vFinalPions.push_back(protonsNoPionSign[0]);
                    vFinalProtons.push_back(protonsNoPionSign[1]);
                    vFinalProtons.push_back(protonsPionSign[0]);                  
                }
                else
                {
                    vFinalPions.push_back(vPions3Sigma[0]);
                    vFinalPions.push_back(protonsNoPionSign[1]);
                    vFinalProtons.push_back(protonsNoPionSign[0]);
                    vFinalProtons.push_back(protonsPionSign[0]);                    
                }
            }

            else
            {
                cout << "HEERERERER " <<endl;
            }
        }
        
        else
        {
            vector <StUPCTrack const *> vProtons2SigmaThatCanBeAlsoPion;
            vector <StUPCTrack const *> vProtons2SigmaThatCanNotBeAlsoPion;
            for (int i = 0; i < vProtons2Sigma.size(); i++)
            {
                float sigmaPion = vProtons2Sigma[i]->getNSigmasTPCPion();
                if (abs(sigmaPion) < 3.0)
                {
                    vProtons2SigmaThatCanBeAlsoPion.push_back(vProtons2Sigma[i]);
                }
                else
                {
                    vProtons2SigmaThatCanNotBeAlsoPion.push_back(vProtons2Sigma[i]);
                }

            }

            if (vProtons2SigmaThatCanBeAlsoPion.size() < 2)  {flag = false; return flag;}
            else if (vProtons2SigmaThatCanBeAlsoPion.size() == 2)
            {
                if ((vProtons2SigmaThatCanBeAlsoPion[0]->getCharge() + vProtons2SigmaThatCanBeAlsoPion[1]->getCharge() ) !=0) {flag = false; return flag;}
                
                TVector3 const tryVec(0,0,0);
                StUPCV0 lVecKaonComb1a(vProtons2SigmaThatCanBeAlsoPion[0], vProtons2SigmaThatCanNotBeAlsoPion[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                StUPCV0 lVecKaonComb1b(vProtons2SigmaThatCanBeAlsoPion[1], vProtons2SigmaThatCanNotBeAlsoPion[1], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);

                StUPCV0 lVecKaonComb2a(vProtons2SigmaThatCanBeAlsoPion[0], vProtons2SigmaThatCanNotBeAlsoPion[1], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                StUPCV0 lVecKaonComb2b(vProtons2SigmaThatCanBeAlsoPion[1], vProtons2SigmaThatCanNotBeAlsoPion[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
            
                double  distSumSquaredComb1 = sqrt(pow((lVecKaonComb1a.m()-ExclusiveK0K0::MASS_LAMBDA), 2) + pow((lVecKaonComb1b.m()-ExclusiveK0K0::MASS_LAMBDA), 2));
                double  distSumSquaredComb2 = sqrt(pow((lVecKaonComb2a.m()-ExclusiveK0K0::MASS_LAMBDA), 2) + pow((lVecKaonComb2b.m()-ExclusiveK0K0::MASS_LAMBDA), 2));

                if (distSumSquaredComb1 < distSumSquaredComb2)
                {
                    vFinalPions.push_back(vProtons2SigmaThatCanBeAlsoPion[0]);
                    vFinalPions.push_back(vProtons2SigmaThatCanBeAlsoPion[1]);
                    vFinalProtons.push_back(vProtons2SigmaThatCanNotBeAlsoPion[0]);
                    vFinalProtons.push_back(vProtons2SigmaThatCanNotBeAlsoPion[1]);                  
                }
                else
                {
                    vFinalPions.push_back(vProtons2SigmaThatCanBeAlsoPion[0]);
                    vFinalPions.push_back(vProtons2SigmaThatCanBeAlsoPion[1]);
                    vFinalProtons.push_back(vProtons2SigmaThatCanNotBeAlsoPion[1]);
                    vFinalProtons.push_back(vProtons2SigmaThatCanNotBeAlsoPion[0]);                    
                }
            }     

            else if (vProtons2SigmaThatCanBeAlsoPion.size() == 3)
            {
                double protonCharge = vProtons2SigmaThatCanNotBeAlsoPion[0]->getCharge();
                vector <StUPCTrack const *> pionsProtonCharge;
                vector <StUPCTrack const *> pionsNoProtonCharge;
                
                for (int i = 0; i < vProtons2SigmaThatCanBeAlsoPion.size(); i++)
                {
                    if (vProtons2SigmaThatCanBeAlsoPion[i]->getCharge() == protonCharge )
                    {
                        pionsProtonCharge.push_back(vProtons2SigmaThatCanBeAlsoPion[i]);
                    }
                    else
                    {
                        pionsNoProtonCharge.push_back(vProtons2SigmaThatCanBeAlsoPion[i]);
                    }
                }

                TVector3 const tryVec(0,0,0);
                StUPCV0 lVecKaonComb1a(pionsNoProtonCharge[0], vProtons2SigmaThatCanNotBeAlsoPion[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                StUPCV0 lVecKaonComb1b(pionsProtonCharge[0], pionsNoProtonCharge[1], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);

                StUPCV0 lVecKaonComb2a(pionsNoProtonCharge[1], vProtons2SigmaThatCanNotBeAlsoPion[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                StUPCV0 lVecKaonComb2b(pionsProtonCharge[0], pionsNoProtonCharge[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
            
                double  distSumSquaredComb1 = sqrt(pow((lVecKaonComb1a.m()-ExclusiveK0K0::MASS_LAMBDA), 2) + pow((lVecKaonComb1b.m()-ExclusiveK0K0::MASS_LAMBDA), 2));
                double  distSumSquaredComb2 = sqrt(pow((lVecKaonComb2a.m()-ExclusiveK0K0::MASS_LAMBDA), 2) + pow((lVecKaonComb2b.m()-ExclusiveK0K0::MASS_LAMBDA), 2));

                if (distSumSquaredComb1 < distSumSquaredComb2)
                {
                    vFinalPions.push_back(pionsNoProtonCharge[0]);
                    vFinalPions.push_back(pionsProtonCharge[0]);
                    vFinalProtons.push_back(vProtons2SigmaThatCanNotBeAlsoPion[0]);
                    vFinalProtons.push_back(pionsNoProtonCharge[1]);                  
                }
                else
                {
                    vFinalPions.push_back(pionsNoProtonCharge[1]);
                    vFinalPions.push_back(pionsProtonCharge[0]);
                    vFinalProtons.push_back(vProtons2SigmaThatCanNotBeAlsoPion[0]);
                    vFinalProtons.push_back(pionsNoProtonCharge[0]);                    
                }            
            }

                   else
            {
                cout << "HEERERERER " <<endl;
            }

        }

            
    }

    else if (tracksWithTofHit.size() == 3)
    {
      
        vector <StUPCTrack const *> vProtons2Sigma_TOF;
        vector <StUPCTrack const *> vProtons2SigmaThatCanBePions_TOF;
        vector <StUPCTrack const *> vProtons2SigmaThatCanNOTBePions_TOF;        
        vector <StUPCTrack const *> vPions3Sigma_TOF;

        vector <StUPCTrack const *> vProtons2Sigma_NOTOF;
        vector <StUPCTrack const *> vProtons2SigmaThatCanBePions_NOTOF;
        vector <StUPCTrack const *> vPions3Sigma_NOTOF;

        for (int i = 0; i < tracksWithTofHit.size(); i++)
        {
            StUPCTrack const * track = tracksWithTofHit[i];
            double sigmaProton = track->getNSigmasTPCProton();
            double sigmaPion = track->getNSigmasTPCPion();
            if (abs(sigmaProton)<2.0)
            {
                vProtons2Sigma_TOF.push_back(track);
                if (abs(sigmaPion) < 3.0)
                {
                   vProtons2SigmaThatCanBePions_TOF.push_back(track);
                }
                else
                {
                    vProtons2SigmaThatCanNOTBePions_TOF.push_back(track);
                }
            }
             else if (abs(sigmaPion) < 3.0)
            {
                vPions3Sigma_TOF.push_back(track);
            }
        }

        for (int i = 0; i < goodTracksWithoutTofHit.size(); i++)
        {
            StUPCTrack const * track = goodTracksWithoutTofHit[i];
            double sigmaProton = track->getNSigmasTPCProton();
            double sigmaPion = track->getNSigmasTPCPion();
            if (abs(sigmaProton)<2.0)
            {
                vProtons2Sigma_NOTOF.push_back(track);
                if (abs(sigmaPion) < 3.0)
                {
                   vProtons2SigmaThatCanBePions_NOTOF.push_back(track);
                   vPions3Sigma_NOTOF.push_back(track);
                }
            }
            else if (abs(sigmaPion) < 3.0)
            {
                vPions3Sigma_NOTOF.push_back(track);
            }
        }
      

        if (vProtons2Sigma_TOF.size() == 0) {flag = false; return flag;}
        if (vProtons2Sigma_TOF.size() == 1)
        {
       
        
            if (vPions3Sigma_TOF.size() != 2) {flag = false; return flag;} // jeżeli nie ma dwóch pionów z TOF - odrzucam przypadek
            if (vPions3Sigma_TOF[0]->getCharge() + vPions3Sigma_TOF[1]->getCharge() != 0 ) {flag = false; return flag;} // odrzucam przypadek jeżeli te piony nie są przecwincyh znaków
            //jednoznaczne dopasowanie na podstawie ładunnku
            int indexOtherPion;

            if (vProtons2Sigma_TOF[0]->getCharge() != vPions3Sigma_TOF[0]->getCharge()) // jeżeli proton i pion mają przeciwne znaki to tworzą lambdę
            {
                vFinalProtons.push_back(vProtons2Sigma_TOF[0]);
                vFinalPions.push_back(vPions3Sigma_TOF[0]);
                indexOtherPion = 1;
            }
            else
            {   
                vFinalProtons.push_back(vProtons2Sigma_TOF[0]);
                vFinalPions.push_back(vPions3Sigma_TOF[1]);
                indexOtherPion = 0;              
            }

            if (vProtons2Sigma_NOTOF.size() == 0)  {flag = false; return flag;} // if there are not tracks w/o TOF that can be Protons
          
            StUPCTrack const * otherPion = vPions3Sigma_TOF[indexOtherPion];
            int chargeOtherPion = otherPion->getCharge();
            vector <StUPCTrack const*> protonsOppoisteToPion_NOTOF;
        
            for (int i = 0; i < vProtons2Sigma_NOTOF.size(); i++ )
            {
                StUPCTrack const * track = vProtons2Sigma_NOTOF[i];
                if (track->getCharge()!=chargeOtherPion ) // jeżeli mają przeciwne ładunki 
                {
                    protonsOppoisteToPion_NOTOF.push_back(track);
                }
            }

            if (protonsOppoisteToPion_NOTOF.size() == 0 )  {flag = false; return flag;} // jeżeli isteniją protony bez TOF które mają przeciny ładunek do pionu
            vector <double> masses;
            TVector3 const tryVec(0,0,0);
            for (int i = 0; i < protonsOppoisteToPion_NOTOF.size(); i++)
            {
                StUPCV0 lambda(otherPion, protonsOppoisteToPion_NOTOF[i], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                masses.push_back( abs(lambda.m() -ExclusiveK0K0::MASS_LAMBDA ));
            }

            auto minElement = std::min_element(masses.begin(), masses.end());
            int minIndex = std::distance(masses.begin(), minElement);
        
            vFinalProtons.push_back(protonsOppoisteToPion_NOTOF[minIndex]);
            vFinalPions.push_back(otherPion);
        }

        else if(vProtons2Sigma_TOF.size() == 2)
        {
            int sumProtonCharge = vProtons2Sigma_TOF[0]->getCharge() + vProtons2Sigma_TOF[1]->getCharge();

            // PODZIAŁ NA PRZYPADKI SUMY ŁADUNKU PROTONÓW
            if (sumProtonCharge == 0)
            {
                if (vPions3Sigma_TOF.size() == 0)  {flag = false; return flag;} //czy istnieje trzeci ślad który jest pionem
                int chargePion = vPions3Sigma_TOF[0]->getCharge();
                int indexOtherProton;
                if (chargePion != vProtons2Sigma_TOF[0]->getCharge())
                {
                    indexOtherProton = 1;
                    vFinalProtons.push_back(vProtons2Sigma_TOF[0]);
                    vFinalPions.push_back(vPions3Sigma_TOF[0]);
                }
                else
                {
                    indexOtherProton = 0;
                    vFinalProtons.push_back(vProtons2Sigma_TOF[1]);
                    vFinalPions.push_back(vPions3Sigma_TOF[0]);
                }

                StUPCTrack const * otherProton = vProtons2Sigma_TOF[indexOtherProton];
                vector <StUPCTrack const *> pionsOppositeToProton_NOTOF;
                for (int i = 0; i < vPions3Sigma_NOTOF.size(); i++)
                {   
                    StUPCTrack const * track = vPions3Sigma_NOTOF[i];
                    if (otherProton->getCharge() != track->getCharge() )
                    {
                        pionsOppositeToProton_NOTOF.push_back(track);
                    }
                }

                if (pionsOppositeToProton_NOTOF.size() == 0 ) {flag = false; return flag;} // jeżeli nie ma pionów o znaku przeciwnym do pozostałego protnu  - opóść

                vector <double> masses;
                TVector3 const tryVec(0,0,0);
                for (int i = 0; i < pionsOppositeToProton_NOTOF.size(); i++)
                {
                    StUPCV0 lambda(pionsOppositeToProton_NOTOF[i], otherProton, ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                    masses.push_back(abs(lambda.m()-ExclusiveK0K0::MASS_LAMBDA));
                }

                auto minElement = std::min_element(masses.begin(), masses.end());
                int minIndex = std::distance(masses.begin(), minElement);

                vFinalPions.push_back(pionsOppositeToProton_NOTOF[minIndex]);
                vFinalProtons.push_back(otherProton);           
            }
            
            else  // protony maja ten sam znak - ale można zachowa przypadek gdy przynajmniej jeden z nich jest pionem
            {   
                if (vProtons2SigmaThatCanBePions_TOF.size() == 0) {flag = false; return flag;} // odzucam przypadek gdy żaden z protonów nie może być pionem
                else if (vProtons2SigmaThatCanBePions_TOF.size() == 1) // przypadek gdy jest dokładnie jeden proton który może być pionem       
                {
                    StUPCTrack const * protonThatCanBePion = vProtons2SigmaThatCanBePions_TOF[0];
                    StUPCTrack const * protonThatCanNOTBePion = vProtons2SigmaThatCanNOTBePions_TOF[0];        
                    if (vPions3Sigma_TOF.size() == 0)  {flag = false; return flag;}  // jeżeli nie ma pionu odrzucam przypadek
                    StUPCTrack const * pion = vPions3Sigma_TOF[0];
                    vFinalPions.push_back(pion); 
                    vFinalProtons.push_back(protonThatCanNOTBePion); //ONE MUSZĄ Z ZAŁOŻENIA MIEĆ PRZECIWNY ŁADUNEK MAM -1 LUB 1 sumaryczny
                    
                    if (vProtons2Sigma_NOTOF.size() == 0)   {flag = false; return flag;}         // jeżeli nie ma kandydatów na proton bez TOF opóść  - opóść

                    vector <StUPCTrack const *> protonsOppositeSignToPPion;
                    for (int i = 0; i <vProtons2Sigma_NOTOF.size(); i++)
                    {
                        StUPCTrack const * track = vProtons2Sigma_NOTOF[i];
                        if (track->getCharge()!=protonThatCanBePion->getCharge())
                        {
                            protonsOppositeSignToPPion.push_back(track);
                        }
                    }

                    if (protonsOppositeSignToPPion.size() == 0) {flag = false; return flag;} // jeżeli nie ma kandydatów na proton o dobrym ładunku to opóść
                    vector <double> masses;
                    TVector3 const tryVec(0,0,0);
                    for (int i = 0; i < protonsOppositeSignToPPion.size(); i++)
                    {
                        StUPCV0 lambda(protonThatCanBePion, protonsOppositeSignToPPion[i], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                        masses.push_back(abs(lambda.m()-ExclusiveK0K0::MASS_LAMBDA));
                    }

                    auto minElement = std::min_element(masses.begin(), masses.end());
                    int minIndex = std::distance(masses.begin(), minElement);
                    vFinalPions.push_back(protonThatCanBePion);
                    vFinalProtons.push_back(protonsOppositeSignToPPion[minIndex]);                    
                }   

                else if (vProtons2SigmaThatCanBePions_TOF.size() == 2) // przypadek gdy jest dokładnie jeden proton który może być pionem       
                {
                    
          
                    if (vPions3Sigma_TOF.size() == 0)  {flag = false; return flag;} // jeżeli nie ma pionu prawdziewgo to bye bye.

                    StUPCTrack const * pion = vPions3Sigma_TOF[0];
                    TVector3 const tryVec(0,0,0);
                    StUPCV0 lambda1(pion, vProtons2SigmaThatCanBePions_TOF[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                    StUPCV0 lambda2(pion, vProtons2SigmaThatCanBePions_TOF[1], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                    
                    int indexOtherPP;
                    if( abs(lambda1.m()-ExclusiveK0K0::MASS_LAMBDA) < abs(lambda2.m()-ExclusiveK0K0::MASS_LAMBDA)  )
                    {
                        vFinalProtons.push_back(vProtons2SigmaThatCanBePions_TOF[0]);
                        vFinalPions.push_back(pion);
                        indexOtherPP = 1;
                    }
                    else
                    {
                        vFinalPions.push_back(pion);
                        vFinalProtons.push_back(vProtons2SigmaThatCanBePions_TOF[1]);
                        indexOtherPP = 0;
                    }
                        

                    StUPCTrack const * otherPP = vProtons2SigmaThatCanBePions_TOF[indexOtherPP];
                    vector <StUPCTrack const *> vOtherProtons;  vector <StUPCTrack const *> vOtherPions;

                   
                    for (int i = 0; i < vPions3Sigma_NOTOF.size(); i++)
                    {
                        if (vPions3Sigma_NOTOF[i]->getCharge() != otherPP->getCharge() )
                        {
                            vOtherPions.push_back(vPions3Sigma_NOTOF[i]) ;
                        }
                    }

                    for (int i = 0; i < vProtons2Sigma_NOTOF.size(); i++)
                    {
                        if (vProtons2Sigma_NOTOF[i]->getCharge() != otherPP->getCharge() )
                        {
                            vOtherProtons.push_back(vProtons2Sigma_NOTOF[i]) ;
                        }
                    }

                 
                    if (vOtherProtons.size() == 0 and vOtherPions.size() == 0)   {flag = false; return flag;}// nie ma komplemenatracyh pionow ani protonow

                    if (vOtherProtons.size() != 0 and vOtherPions.size() != 0)
                    {
                        vector <double> massesPions;
                        vector <double> massesProtons;   
                        for (int i = 0; i < vOtherPions.size(); i++)
                        {
                            StUPCTrack const * track = vOtherPions[i];
                            StUPCV0 lambda1(track, otherPP, ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                            massesPions.push_back(  abs(lambda1.m() - ExclusiveK0K0::MASS_LAMBDA)  );
                        }
                        for (int i = 0; i < vOtherProtons.size(); i++)
                        {
                            StUPCTrack const * track = vOtherProtons[i];
                            StUPCV0 lambda2(otherPP, track, ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                            massesProtons.push_back(  abs(lambda2.m() - ExclusiveK0K0::MASS_LAMBDA)  );
                        }
                      
                        auto minElementPion = std::min_element(massesPions.begin(), massesPions.end());
                        int minIndexPion = std::distance(massesPions.begin(), minElementPion);

                        auto minElementProton = std::min_element(massesProtons.begin(), massesProtons.end());
                        int minIndexProton = std::distance(massesProtons.begin(), minElementProton);

                        if (massesProtons[minIndexProton] <massesPions[minIndexPion] )
                        {
                            vFinalPions.push_back(vOtherPions[minIndexPion]);
                            vFinalProtons.push_back(otherPP);
                        }
                        else
                        {
                            vFinalPions.push_back(otherPP);
                            vFinalProtons.push_back(vOtherProtons[minIndexProton]);
                        }
                     
                    }

                    else if (vOtherPions.size() != 0)
                    {
                        vector <double> massesPions;
                        for (int i = 0; i < vOtherPions.size(); i++)
                        {
                            StUPCTrack const * track = vOtherPions[i];
                            StUPCV0 lambda1(track, otherPP, ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                            massesPions.push_back(  abs(lambda1.m() - ExclusiveK0K0::MASS_LAMBDA)  );
                        }
                        auto minElementPion = std::min_element(massesPions.begin(), massesPions.end());
                        int minIndexPion = std::distance(massesPions.begin(), minElementPion);
                        vFinalProtons.push_back(otherPP);
                        vFinalPions.push_back(vOtherPions[minIndexPion]);
                    }

                    else if (vOtherProtons.size() != 0)
                    {  vector <double> massesProtons;   
                        for (int i = 0; i < vOtherProtons.size(); i++)
                        {
                            StUPCTrack const * track =  vOtherProtons[i];
                            StUPCV0 lambda2(otherPP, track, ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                            massesProtons.push_back(  abs(lambda2.m() - ExclusiveK0K0::MASS_LAMBDA)  );
                        }
                        auto minElementProton = std::min_element(massesProtons.begin(), massesProtons.end());
                        int minIndexProton = std::distance(massesProtons.begin(), minElementProton);
                        vFinalPions.push_back(otherPP);
                        vFinalProtons.push_back(vOtherProtons[minIndexProton]);
                    }
                    
                }
                
            } // nieparzysty ładunek dla 2 protonów mogących być pioniami
        }

        else if (vProtons2Sigma_TOF.size() == 3)
        {   
           
        
            // wszystkie TOF są protnami, czyli potzebuję przynajmniej jedneo pionu
            if (vProtons2SigmaThatCanBePions_TOF.size() == 0) {flag = false; return flag;}
            int sumCharge = vProtons2Sigma_TOF[0]->getCharge() +  vProtons2Sigma_TOF[1]->getCharge() +  vProtons2Sigma_TOF[2]->getCharge();
            if (vProtons2SigmaThatCanBePions_TOF.size() == 0) {flag = false; return flag;}

            else if (vProtons2SigmaThatCanBePions_TOF.size() == 1)
            {
                StUPCTrack const * pion = vProtons2SigmaThatCanBePions_TOF[0];
                StUPCTrack const * proton1 = vProtons2SigmaThatCanNOTBePions_TOF[0];
                StUPCTrack const * proton2 = vProtons2SigmaThatCanNOTBePions_TOF[1];
                                
                int sumCharge = pion->getCharge() + proton1->getCharge() + proton2->getCharge();


                 if (pion->getCharge() == sumCharge) // 1 1 -1 np. pion ma 1, to jest parowany z protonem -1,  czyli szukamy pionu bez tof o ładunku przeciwnym do sumarycznego! do pionu z TOF o ładu sum
                {
                        int ind;
                        if (proton1->getCharge() == sumCharge)
                        {
                            ind = 0;
                        }

                        else
                        {
                            ind = 1;
                        }
                        StUPCTrack const * protonOfPionCharge = vProtons2SigmaThatCanNOTBePions_TOF[ind];
                        if (vPions3Sigma_NOTOF.size() == 0) {flag = false; return flag;} // brak śladów bez TOF kótre mogłybyb być pionami 
                        vector <StUPCTrack const *> vPions3Sigma_NOTOF_GOODcharge;
                        
                        for (int i = 0; i < vPions3Sigma_NOTOF.size(); i++)
                        {
                            StUPCTrack const * track = vPions3Sigma_NOTOF[i];
                            if (track->getCharge() != sumCharge )
                            {
                                vPions3Sigma_NOTOF_GOODcharge.push_back(track);
                            }
                        }    
                        if (vPions3Sigma_NOTOF_GOODcharge.size() == 0) {flag = false; return flag;}  // jeżeli brak pionów bez TOF mających dobry ładunek                   

                        vector <double> masses;
                        TVector3 const tryVec(0,0,0);
                        for (int i = 0; i < vPions3Sigma_NOTOF_GOODcharge.size(); i++)
                        {
                            StUPCV0 lambda(vPions3Sigma_NOTOF_GOODcharge[i], protonOfPionCharge, ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                            masses.push_back( abs(lambda.m() -ExclusiveK0K0::MASS_LAMBDA ));
                        }

                        auto minElement = std::min_element(masses.begin(), masses.end());
                        int minIndex = std::distance(masses.begin(), minElement);

                        vFinalPions.push_back(vPions3Sigma_NOTOF_GOODcharge[minIndex]);
                        vFinalPions.push_back(pion);
                        vFinalProtons.push_back(proton1);
                        vFinalProtons.push_back(proton2);
                }

                else
                {
                    flag = false; return flag;
                }
            }


            else if (vProtons2SigmaThatCanBePions_TOF.size() == 2)
            {
                StUPCTrack const * proton = vProtons2SigmaThatCanNOTBePions_TOF[0];
                StUPCTrack const * pion1 = vProtons2SigmaThatCanBePions_TOF[0];
                StUPCTrack const * pion2 = vProtons2SigmaThatCanBePions_TOF[1];
                 
                if (pion1->getCharge() != pion2->getCharge() )
                {

                        int ind;
                        if (pion1->getCharge() == sumCharge)
                        {
                            ind = 0;
                        }

                        else
                        {
                            ind = 1;
                        }

                        StUPCTrack const * otherPion = vProtons2SigmaThatCanBePions_TOF[ind];
                        if (vProtons2Sigma_NOTOF.size() == 0) {flag = false; return flag;} // brak śladów bez TOF kótre mogłybyb być pionami 
                        vector <StUPCTrack const *> vProtons2Sigma_NOTOF_GOODcharge;
                        
                        for (int i = 0; i < vPions3Sigma_NOTOF.size(); i++)
                        {
                            StUPCTrack const * track = vPions3Sigma_NOTOF[i];
                            if (track->getCharge() != sumCharge )
                            {
                                vProtons2Sigma_NOTOF_GOODcharge.push_back(track);
                            }
                        }    
                        if (vProtons2Sigma_NOTOF_GOODcharge.size() == 0) {flag = false; return flag;}  // jeżeli brak pionów bez TOF mających dobry ładunek                   

                        vector <double> masses;
                        TVector3 const tryVec(0,0,0);
                        for (int i = 0; i < vProtons2Sigma_NOTOF_GOODcharge.size(); i++)
                        {
                            StUPCV0 lambda(otherPion, vProtons2Sigma_NOTOF_GOODcharge[i], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                            masses.push_back( abs(lambda.m() -ExclusiveK0K0::MASS_LAMBDA ));
                        }

                        auto minElement = std::min_element(masses.begin(), masses.end());
                        int minIndex = std::distance(masses.begin(), minElement);

                        vFinalProtons.push_back(vProtons2Sigma_NOTOF_GOODcharge[minIndex]);
                        vFinalProtons.push_back(proton);
                        vFinalPions.push_back(pion1);
                        vFinalPions.push_back(pion2);            
                }

                else // pion 1, pion1, proton -1
                {
                    
          
                    TVector3 const tryVec(0,0,0);
                    StUPCV0 lambda1(pion1, proton, ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                    StUPCV0 lambda2(pion2, proton, ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                    

                    int indexOtherPP;
                    if( abs(lambda1.m()-ExclusiveK0K0::MASS_LAMBDA) < abs(lambda2.m()-ExclusiveK0K0::MASS_LAMBDA)  )
                    {
                        vFinalProtons.push_back(proton);
                        vFinalPions.push_back(pion1);
                        indexOtherPP = 1;
                    }
                    else
                    {
                        vFinalPions.push_back(pion2);
                        vFinalProtons.push_back(proton);
                        indexOtherPP = 0;
                    }
                        

                    StUPCTrack const * otherPP = vProtons2SigmaThatCanBePions_TOF[indexOtherPP];
                    vector <StUPCTrack const *> vOtherProtons;  vector <StUPCTrack const *> vOtherPions;

                   
                    for (int i = 0; i < vPions3Sigma_NOTOF.size(); i++)
                    {
                        if (vPions3Sigma_NOTOF[i]->getCharge() != otherPP->getCharge() )
                        {
                            vOtherPions.push_back(vPions3Sigma_NOTOF[i]) ;
                        }
                    }

                    if (vOtherPions.size() == 0) {flag = false; return false;}
                    
                    vector <double> massesPions;
                    for (int i = 0; i < vOtherPions.size(); i++)
                    {
                        StUPCTrack const * track = vOtherPions[i];
                        StUPCV0 lambda1(track, otherPP, ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                        massesPions.push_back(  abs(lambda1.m() - ExclusiveK0K0::MASS_LAMBDA)  );
                    }
                    auto minElementPion = std::min_element(massesPions.begin(), massesPions.end());
                    int minIndexPion = std::distance(massesPions.begin(), minElementPion);
                    vFinalProtons.push_back(otherPP);
                    vFinalPions.push_back(vOtherPions[minIndexPion]);
                    
                }


            }

            else 
            {  
                int charge = vProtons2SigmaThatCanBePions_TOF[0]->getCharge() + vProtons2SigmaThatCanBePions_TOF[1]->getCharge() + vProtons2SigmaThatCanBePions_TOF[2]->getCharge() ;
                vector <StUPCTrack const *> vProtonsThatCanBePionsSameCharge; 
                vector <StUPCTrack const *> vProtonsThatCanBePionsUniqeCharge; 

                for (int i = 0; i < 3; i++)
                {
                    if (vProtons2SigmaThatCanBePions_TOF[i]->getCharge() == charge)
                    {
                        vProtonsThatCanBePionsSameCharge.push_back(vProtons2SigmaThatCanBePions_TOF[i]);
                    }

                    else
                    {
                        vProtonsThatCanBePionsUniqeCharge.push_back(vProtons2SigmaThatCanBePions_TOF[i]);
                    }
                }
                
                 TVector3 const tryVec(0,0,0);
                StUPCV0 lVecKaonComb1a(vProtonsThatCanBePionsSameCharge[0], vProtonsThatCanBePionsUniqeCharge[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                StUPCV0 lVecKaonComb1b(vProtonsThatCanBePionsSameCharge[1], vProtonsThatCanBePionsUniqeCharge[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                StUPCV0 lVecKaonComb2a(vProtonsThatCanBePionsSameCharge[0], vProtonsThatCanBePionsUniqeCharge[0], ExclusiveK0K0::MASS_PROTON, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                StUPCV0 lVecKaonComb2b(vProtonsThatCanBePionsSameCharge[1], vProtonsThatCanBePionsUniqeCharge[0], ExclusiveK0K0::MASS_PROTON, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);

                double m1 = abs(lVecKaonComb1a.m() - ExclusiveK0K0::MASS_LAMBDA);
                double m2 = abs(lVecKaonComb1b.m() - ExclusiveK0K0::MASS_LAMBDA);
                double m3 = abs(lVecKaonComb2a.m() - ExclusiveK0K0::MASS_LAMBDA);
                double m4 = abs(lVecKaonComb2b.m() - ExclusiveK0K0::MASS_LAMBDA);

                vector <double> vmas; vmas.push_back(m1); vmas.push_back(m2); vmas.push_back(m3); vmas.push_back(m4);
      
                
                auto minElement = std::min_element(vmas.begin(), vmas.end());
                int minIndex = std::distance(vmas.begin(), minElement);
                int isProt = 0;
                int indexOtherProton;
                if (minIndex == 0)
                {
                    vFinalPions.push_back(vProtonsThatCanBePionsSameCharge[0]);
                    vFinalProtons.push_back(vProtonsThatCanBePionsUniqeCharge[0]);
                    indexOtherProton=1;
                }
                else if (minIndex == 1)
                {
                    vFinalPions.push_back(vProtonsThatCanBePionsSameCharge[1]);
                    vFinalProtons.push_back(vProtonsThatCanBePionsUniqeCharge[0]);
                    indexOtherProton=0;
                }
                else if (minIndex == 2)
                {
                    vFinalProtons.push_back(vProtonsThatCanBePionsSameCharge[0]);
                    vFinalPions.push_back(vProtonsThatCanBePionsUniqeCharge[0]);
                    isProt=1;
                    indexOtherProton=1;
                }

                else
                {
                    vFinalProtons.push_back(vProtonsThatCanBePionsSameCharge[1]);
                    vFinalPions.push_back(vProtonsThatCanBePionsUniqeCharge[0]);
                    isProt=1;
                    indexOtherProton=0;
                }

                vector <double> masses;             
                vector <StUPCTrack const *> vPions3Sigma_NOTOF_GOODcharge;
                 vector <StUPCTrack const *> vProtons3Sigma_NOTOF_GOODcharge;
                for (int i = 0; i < vPions3Sigma_NOTOF.size(); i++)
                {
                    if (vPions3Sigma_NOTOF[i]->getCharge()!= charge)
                    {
                        vPions3Sigma_NOTOF_GOODcharge.push_back(vPions3Sigma_NOTOF[i]);
                    }
                }

                if (vPions3Sigma_NOTOF_GOODcharge.size()== 0) {flag = false; return flag;}

                if (isProt == 0)
                {

                    for (int i = 0; i < vPions3Sigma_NOTOF_GOODcharge.size(); i++)
                    {
                        StUPCV0 lambda(vPions3Sigma_NOTOF_GOODcharge[i], vProtonsThatCanBePionsSameCharge[indexOtherProton], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PROTON, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                        masses.push_back( abs(lambda.m() -ExclusiveK0K0::MASS_LAMBDA ));
                    }

                    auto minElement2 = std::min_element(masses.begin(), masses.end());
                    int minIndex2 = std::distance(masses.begin(), minElement2);
                    vFinalProtons.push_back(vProtonsThatCanBePionsSameCharge[indexOtherProton]);
                    vFinalPions.push_back(vPions3Sigma_NOTOF_GOODcharge[minIndex2]);
                }

                else
                {
                    for (int i = 0; i < vProtons2Sigma_NOTOF.size(); i++)
                    {
                        if (vProtons2Sigma_NOTOF[i]->getCharge()!= charge)
                        {
                            vProtons3Sigma_NOTOF_GOODcharge.push_back(vProtons2Sigma_NOTOF[i]);
                        }
                    }

                    for (int i = 0; i < vProtons3Sigma_NOTOF_GOODcharge.size(); i++)
                    {
                        StUPCV0 lambda(vProtons3Sigma_NOTOF_GOODcharge[i], vProtonsThatCanBePionsSameCharge[indexOtherProton], ExclusiveK0K0::MASS_PROTON, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                        masses.push_back( abs(lambda.m() -ExclusiveK0K0::MASS_LAMBDA ));
                    }

                    auto minElement2 = std::min_element(masses.begin(), masses.end());
                    int minIndex2 = std::distance(masses.begin(), minElement2);
                    vFinalPions.push_back(vProtonsThatCanBePionsSameCharge[indexOtherProton]);
                    vFinalProtons.push_back(vProtons3Sigma_NOTOF_GOODcharge[minIndex2]);

                }
                cout << vFinalPions[0]->getCharge() << vFinalPions[1]->getCharge() << vFinalProtons[0]->getCharge() << vFinalProtons[1]->getCharge() << endl;

                flag = true;                
            }         
        }   
        else
        {
            flag = false;
        }

    }
    else
    {
        flag = false;return flag;
    }

    return flag;
}



bool  FindPions(  StUPCEvent *upcEvt, vector <StUPCTrack const*> &tracksWithTofHit, vector <StUPCTrack const*> &tracksWithoutTofHit, vector <StUPCTrack const*> &goodTracksWithoutTofHit,vector <StUPCTrack const*> &  vPosNegPionLeadingKaon, vector <StUPCTrack const*> & vPosNegPionSubLeadingKaon, double * beamPar)
{
    bool flag = true;
    vector <StUPCTrack const *> posPion;
    vector <StUPCTrack const *> negPion;

    if (tracksWithTofHit.size() == 5)
    {
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {
            if( tracksWithTofHit[i]->getCharge() == 1)
            {
                posPion.push_back(tracksWithTofHit[i]);
            }
            else if (tracksWithTofHit[i]->getCharge() == -1)
            {
                negPion.push_back(tracksWithTofHit[i]);
            }
        }  

        vector  <StUPCTrack const *> vec3;
        vector  <StUPCTrack const *> vec2; 

        if (negPion.size() == 3) { vec3 = negPion; vec2 = posPion;}    
        else {vec2 = negPion; vec3 = posPion;}    

        vector <vector <int>> indices;
        vector <double> masses;
        vector <vector<StUPCTrack const *>> tracks5;
        TVector3 const tryVec(0,0,0);

        for (int i = 0; i < vec3.size(); i++)
        {
            for (int j = 0; j < vec3.size(); j++)
            {
                if (i!=j)
                {
                    vector <StUPCTrack const *> t;
                    StUPCV0 *lVecKaonComb1 = new StUPCV0(vec3[i], vec2[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(),  true);
                    StUPCV0 *lVecKaonComb2 = new StUPCV0(vec3[j], vec2[1], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
                    double dcaDau1 = lVecKaonComb1->dcaDaughters();
                    double dcaDau2 = lVecKaonComb2->dcaDaughters();
                    double dcaBeam1 = lVecKaonComb1->DCABeamLine();
                    double dcaBeam2 = lVecKaonComb2->DCABeamLine();

                    double cos1 = lVecKaonComb1->pointingAngleHypo();
                    double cos2 = lVecKaonComb2->pointingAngleHypo();
                    double len1 = lVecKaonComb1->decayLengthHypo();
                    double len2 = lVecKaonComb2->decayLengthHypo();

                    if (dcaDau1 < 2.5 and dcaDau2 < 2.5 and dcaBeam1 < 2.5 and dcaBeam2 < 2.5  and (cos1>0.925 or len1 < 3.0) and  (cos2>0.925 or len2 < 3.0) )
                    {
                        double distMass = sqrt(pow((lVecKaonComb1->m()-ExclusiveK0K0::MASS_KAON), 2) + pow((lVecKaonComb2->m()-ExclusiveK0K0::MASS_KAON), 2));
                        masses.push_back(distMass);
                        t.push_back(vec3[i]);
                        t.push_back(vec3[j]);
                        tracks5.push_back(t);
                    }
                }
            }
        }

        if (masses.size() == 0)
        {
            flag = false;
            return flag;
        }
           
        auto minElement = std::min_element(masses.begin(), masses.end());
        int minIndex = std::distance(masses.begin(), minElement);

        vector <StUPCTrack const *> d = tracks5[minIndex];
        StUPCV0 kaon1(d[0], vec2[0],ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(),  true);                
        StUPCV0 kaon2(d[1], vec2[1],ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(),  true);   
        
        if (kaon1.pt() > kaon2.pt())
        {
            vPosNegPionLeadingKaon.push_back(d[0]); vPosNegPionLeadingKaon.push_back(vec2[0]);
            vPosNegPionSubLeadingKaon.push_back(d[1]); vPosNegPionSubLeadingKaon.push_back(vec2[1]);
        }

        else
        {
            vPosNegPionSubLeadingKaon.push_back(d[0]); vPosNegPionSubLeadingKaon.push_back(vec2[0]);
            vPosNegPionLeadingKaon.push_back(d[1]); vPosNegPionLeadingKaon.push_back(vec2[1]);
        }
    }

    else if (tracksWithTofHit.size() == 4)
    {
        double distSumSquaredComb1, distSumSquaredComb2;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {
            if( tracksWithTofHit[i]->getCharge() == 1)
            {
                posPion.push_back(tracksWithTofHit[i]);
            }
            else if (tracksWithTofHit[i]->getCharge() == -1)
            {
                negPion.push_back(tracksWithTofHit[i]);
            }
        }

        TVector3 const tryVec(0,0,0);
        StUPCV0 lVecKaonComb1a(posPion[0], negPion[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(),true);
        StUPCV0 lVecKaonComb1b(posPion[1], negPion[1], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
        StUPCV0 lVecKaonComb2a(posPion[0], negPion[1], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
        StUPCV0 lVecKaonComb2b(posPion[1], negPion[0], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
    
        distSumSquaredComb1 = sqrt(pow((lVecKaonComb1a.m()-ExclusiveK0K0::MASS_KAON), 2) + pow((lVecKaonComb1b.m()-ExclusiveK0K0::MASS_KAON), 2));
        distSumSquaredComb2 = sqrt(pow((lVecKaonComb2a.m()-ExclusiveK0K0::MASS_KAON), 2) + pow((lVecKaonComb2b.m()-ExclusiveK0K0::MASS_KAON), 2));

        if (distSumSquaredComb1 < distSumSquaredComb2)
        {
            if (lVecKaonComb1a.pt() > lVecKaonComb1b.pt())
            {
                vPosNegPionLeadingKaon.push_back(posPion[0]);
                vPosNegPionLeadingKaon.push_back(negPion[0]);
                vPosNegPionSubLeadingKaon.push_back(posPion[1]);
                vPosNegPionSubLeadingKaon.push_back(negPion[1]);
            }
            else
            {
                vPosNegPionLeadingKaon.push_back(posPion[1]);
                vPosNegPionLeadingKaon.push_back(negPion[1]);
                vPosNegPionSubLeadingKaon.push_back(posPion[0]);
                vPosNegPionSubLeadingKaon.push_back(negPion[0]);
            }
        }

        else
        {
            if (lVecKaonComb2a.pt() > lVecKaonComb2b.pt())
            {
                vPosNegPionLeadingKaon.push_back(posPion[0]);
                vPosNegPionLeadingKaon.push_back(negPion[1]);
                vPosNegPionSubLeadingKaon.push_back(posPion[1]);
                vPosNegPionSubLeadingKaon.push_back(negPion[0]);
            }
            else
            {
                vPosNegPionLeadingKaon.push_back(posPion[1]);
                vPosNegPionLeadingKaon.push_back(negPion[0]);
                vPosNegPionSubLeadingKaon.push_back(posPion[0]);
                vPosNegPionSubLeadingKaon.push_back(negPion[1]);
            }
        }
    }

    else if (tracksWithTofHit.size() == 3)
    {
        int complementaryCharge = goodTracksWithoutTofHit[0]->getCharge();
        int sumTofCharge = tracksWithTofHit[0]->getCharge() + tracksWithTofHit[1]->getCharge() + tracksWithTofHit[2]->getCharge();

        vector <StUPCTrack const *> twoFoundTracks;
        vector <StUPCTrack const *> oneFoundTrack;

        for (int i = 0; i < tracksWithTofHit.size(); i++)
        {
            if (tracksWithTofHit[i]->getCharge() == complementaryCharge)
            {
                oneFoundTrack.push_back(tracksWithTofHit[i]);
            }
            
            else
            {
                twoFoundTracks.push_back(tracksWithTofHit[i]);
            }
        }

        vector <double> dist1;
        vector <double> dist2;
        TVector3 const tryVec(0,0,0);

        for (int i = 0; i < goodTracksWithoutTofHit.size(); i++)
        {
            StUPCTrack const* tempTrack =  goodTracksWithoutTofHit[i];
            StUPCV0 lVecKaonComb1a(twoFoundTracks[0], oneFoundTrack[0],ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(),  true);
            StUPCV0 lVecKaonComb1b(twoFoundTracks[1], tempTrack,ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(),  true); 
           
            StUPCV0 lVecKaonComb2a(twoFoundTracks[1], oneFoundTrack[0],ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(),  true);
            StUPCV0 lVecKaonComb2b(twoFoundTracks[0], tempTrack,ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(),  true); 
                  
            double massDist1 = sqrt(pow((lVecKaonComb1a.m()-ExclusiveK0K0::MASS_KAON),2) + pow((lVecKaonComb1b.m()-ExclusiveK0K0::MASS_KAON),2));   
            double massDist2 = sqrt(pow((lVecKaonComb2a.m()-ExclusiveK0K0::MASS_KAON),2) + pow((lVecKaonComb2b.m()-ExclusiveK0K0::MASS_KAON),2));   
            dist1.push_back(massDist1);
            dist2.push_back(massDist2);   
        } 

        auto minElement1 = std::min_element(dist1.begin(), dist1.end());
        int minIndex1 = std::distance(dist1.begin(), minElement1);

        auto minElement2 = std::min_element(dist2.begin(), dist2.end());
        int minIndex2 = std::distance(dist2.begin(), minElement2);

        if (dist1[minIndex1] < dist2[minIndex2])
        {
            vPosNegPionLeadingKaon.push_back(twoFoundTracks[0]); vPosNegPionLeadingKaon.push_back(oneFoundTrack[0]);
            vPosNegPionSubLeadingKaon.push_back(twoFoundTracks[1]); vPosNegPionSubLeadingKaon.push_back(goodTracksWithoutTofHit[minIndex1]);
        }

        else
        {
            vPosNegPionLeadingKaon.push_back(twoFoundTracks[1]); vPosNegPionLeadingKaon.push_back(oneFoundTrack[0]);
            vPosNegPionSubLeadingKaon.push_back(twoFoundTracks[0]); vPosNegPionSubLeadingKaon.push_back(goodTracksWithoutTofHit[minIndex2]);
        }
    }

    else if (tracksWithTofHit.size() == 2)
    {
        int sumCharge = tracksWithTofHit[0]->getCharge() + tracksWithTofHit[1]->getCharge();
        if (sumCharge == 0)
        {    
            vector <StUPCTrack const *> firstGoodTracksWithoutTofHit; vector <StUPCTrack const *> secondGoodTracksWithoutTofHit;
            for (int i = 0; i < goodTracksWithoutTofHit.size(); i++)
            {
                if (goodTracksWithoutTofHit[i]->getCharge() == -1*tracksWithTofHit[0]->getCharge())
                {
                    firstGoodTracksWithoutTofHit.push_back(goodTracksWithoutTofHit[i]);
                }
                else
                {
                    secondGoodTracksWithoutTofHit.push_back(goodTracksWithoutTofHit[i]);
                }
            }


            vector <double> dist1; vector <double> dist2; TVector3 const tryVec(0,0,0);
            vector <double> pt1; vector <double> pt2;
            for (int i = 0; i < firstGoodTracksWithoutTofHit.size(); i++) 
            {
                StUPCV0 lVecKaonComb1(tracksWithTofHit[0], firstGoodTracksWithoutTofHit[i],ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(),  true); 
                double massDist = pow((lVecKaonComb1.m()-ExclusiveK0K0::MASS_KAON),2);
                dist1.push_back(massDist);
                pt1.push_back(lVecKaonComb1.pt());
            }

            for (int i = 0; i < secondGoodTracksWithoutTofHit.size(); i++)
            {
                StUPCV0 lVecKaonComb2(tracksWithTofHit[1], secondGoodTracksWithoutTofHit[i],ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(),  true); 
                double massDist = pow((lVecKaonComb2.m()-ExclusiveK0K0::MASS_KAON),2);
                dist2.push_back(massDist);
                pt2.push_back(lVecKaonComb2.pt());
            }

            auto minElement1 = std::min_element(dist1.begin(), dist1.end());
            int minIndex1 = std::distance(dist1.begin(), minElement1);

            auto minElement2 = std::min_element(dist2.begin(), dist2.end());
            int minIndex2 = std::distance(dist2.begin(), minElement2);     

            double pt1Max = pt1[minIndex1]; double pt2Max = pt2[minIndex2];
            
            if (pt1Max >= pt2Max)
            {
                vPosNegPionLeadingKaon.push_back(tracksWithTofHit[0]);  vPosNegPionLeadingKaon.push_back(firstGoodTracksWithoutTofHit[minIndex1]);
                vPosNegPionSubLeadingKaon.push_back(tracksWithTofHit[1]); vPosNegPionSubLeadingKaon.push_back(secondGoodTracksWithoutTofHit[minIndex2]);
            }

            else
            {
                vPosNegPionSubLeadingKaon.push_back(tracksWithTofHit[0]);  vPosNegPionSubLeadingKaon.push_back(firstGoodTracksWithoutTofHit[minIndex1]);
                vPosNegPionLeadingKaon.push_back(tracksWithTofHit[1]); vPosNegPionLeadingKaon.push_back(secondGoodTracksWithoutTofHit[minIndex2]);
            }
       }

        if (sumCharge == 2 or sumCharge == -2)
        {
            vector <double> dist;
            vector <vector <int>> ind; 
            vector <vector <double>> pt; 
            TVector3 const tryVec(0,0,0);
            for (int i = 0 ; i < goodTracksWithoutTofHit.size(); i++)
            {
                for (int j = 0; j < goodTracksWithoutTofHit.size(); j++)
                {
                    if (i!=j)
                    {
                        StUPCV0 lVecKaonComb1(tracksWithTofHit[0], goodTracksWithoutTofHit[i],ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(),  true); 
                        StUPCV0 lVecKaonComb2(tracksWithTofHit[1], goodTracksWithoutTofHit[j],ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(),  true); 
                        double massDist = sqrt(pow((lVecKaonComb1.m()-ExclusiveK0K0::MASS_KAON),2) + pow((lVecKaonComb2.m()-ExclusiveK0K0::MASS_KAON),2));   
                        dist.push_back(massDist);
                        vector <int> in; in.push_back(i); in.push_back(j);
                        vector <double> p; p.push_back(lVecKaonComb1.pt()); p.push_back(lVecKaonComb2.pt());
                        ind.push_back(in);
                        pt.push_back(p);
                    }
                }
            }
            auto minElement = std::min_element(dist.begin(), dist.end());
            int minIndex = std::distance(dist.begin(), minElement);
            int iMin = ind[minIndex][0]; int jMin = ind[minIndex][1];
            double pt1Max = pt[minIndex][0]; double pt2Max = pt[minIndex][1];

            if (pt1Max >= pt2Max)
            {
                vPosNegPionLeadingKaon.push_back(tracksWithTofHit[0]);  vPosNegPionLeadingKaon.push_back(goodTracksWithoutTofHit[iMin]);
                vPosNegPionSubLeadingKaon.push_back(tracksWithTofHit[1]); vPosNegPionSubLeadingKaon.push_back(goodTracksWithoutTofHit[jMin]);         
            }

            else
            {
                vPosNegPionSubLeadingKaon.push_back(tracksWithTofHit[0]);  vPosNegPionSubLeadingKaon.push_back(goodTracksWithoutTofHit[iMin]);
                vPosNegPionLeadingKaon.push_back(tracksWithTofHit[1]); vPosNegPionLeadingKaon.push_back(goodTracksWithoutTofHit[jMin]);
            }
        }
    }
    return flag;
}


bool CheckPtMiss(StUPCV0 &leadingKaon, StUPCV0 &subLeadingKaon, TLorentzVector protonE, TLorentzVector protonW, double & pTmiss)
{
    double pXmiss = leadingKaon.px() + subLeadingKaon.px() + protonE.X() + protonW.X();
    double pYmiss = leadingKaon.py() + subLeadingKaon.py()  + protonE.Y() + protonW.Y();
    pTmiss = sqrt(pow(pXmiss,2) + pow(pYmiss,2));   
    bool flag = 0;
    if (pTmiss <= ExclusiveK0K0::PT_MISS) 
    {
        flag = 1;
    }
    return flag;
}

bool CheckNumberOfClusters(  StUPCEvent *upcEvt, vector <StUPCTrack const*> &tracksWithTofHit, int &  totalCluster)
{
 
    Int_t nTofHits = upcEvt->getNumberOfHits();           
    vector <Int_t> vTray;
    vector <Int_t> vTrayUniqueVector;
    vector <Int_t> vTrayUniqueSet;
    vector <Int_t> vMmodule;
    
    for (Int_t i = 0; i < nTofHits; i++)
    {
        vTray.push_back(Int_t(upcEvt->getHit(i)->getTray()));
        vTrayUniqueVector.push_back(Int_t(upcEvt->getHit(i)->getTray()));
        vMmodule.push_back(Int_t(upcEvt->getHit(i)->getModule()));
    }

    sort(vTrayUniqueVector.begin(), vTrayUniqueVector.end());   
    auto last = unique(vTrayUniqueVector.begin(), vTrayUniqueVector.end());   
    vTrayUniqueVector.erase(last, vTrayUniqueVector.end()); 
    

    vector <vector <Int_t>> vModuleUniqueTray;
    for (int i = 0; i < vTrayUniqueVector.size(); i++)
    {
        vector <Int_t> vModuleUnique;
        for (int j = 0; j < nTofHits; j++)
        {
            if (vTrayUniqueVector[i] == vTray[j] )
            {
                vModuleUnique.push_back(vMmodule[j]);
            }
        }
        vModuleUniqueTray.push_back(vModuleUnique);
        vModuleUnique.clear();
    }

    
    for (int i  = 0; i < vModuleUniqueTray.size(); i++)
    {
        vector <Int_t> vec =  vModuleUniqueTray[i] ;  
        sort(vec.begin(), vec.end());   
        auto last = unique(vec.begin(), vec.end());   
        vec.erase(last, vec.end()); 
    
        if (vec.size() == 1)
        {
            totalCluster+=1;
        
        } 

        
        for (int j = 0; j < vec.size()-1; j++)
        {
            Int_t modNum = vec[j];
            if (vec.size() == 1)
            cout << vec.size() << endl;
            int diff = 1;
            int num = 0;
            for (int z = j+1; z < vec.size(); z++)
            {
                if (modNum+diff == vec[z])
                {
                    num+=0;
                    diff+=1;

                    if (z == j+1)
                    {
                        num+=1;
                    }
                }

                else if ( j == (vec.size()-2) and  vec[vec.size()-2]+1 != vec[vec.size()-1])
                {
                    num+=2; 
                    continue;
                }

                else
                {
                    num+=1;
                }

                j+=(diff);
            }
        
            totalCluster+=num;
        }
    
    }

    bool flag = 0;

    if (totalCluster  <= 2*tracksWithTofHit.size() + 1)
    {
        flag = 1;
    }
    return flag;
}
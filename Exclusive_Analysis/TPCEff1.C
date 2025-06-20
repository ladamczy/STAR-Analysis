// Script to check loading and interpolation of 3D histograms from a ROOT file
#include <TFile.h>
#include <TH3F.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include <iomanip>

namespace Util {
  enum SIGN { PLUS, MINUS, nSigns };
  enum PARTICLE_NAME { PION, KAON, PROTON, nDefinedParticles };
}

class EfficiencyChecker {
private:
  TH3F* mhTpcRecoEff3D[Util::nSigns][Util::nDefinedParticles];
  TString histNames[Util::nSigns][Util::nDefinedParticles];
  
  void InitializeHistNames() {
    histNames[Util::MINUS][Util::PION] = TString("0");
    histNames[Util::MINUS][Util::KAON] = TString("1");
    histNames[Util::MINUS][Util::PROTON] = TString("2");
    histNames[Util::PLUS][Util::PION] = TString("3");
    histNames[Util::PLUS][Util::KAON] = TString("4");
    histNames[Util::PLUS][Util::PROTON] = TString("5");
  }
  
public:
  EfficiencyChecker() {
    InitializeHistNames();
    // Initialize histogram pointers to NULL
    for(int i=0; i<Util::nSigns; ++i) {
      for(int j=0; j<Util::nDefinedParticles; ++j) {
        mhTpcRecoEff3D[i][j] = NULL;
      }
    }
  }
  
  ~EfficiencyChecker() {
    // Cleaning up histograms
    for(int i=0; i<Util::nSigns; ++i) {
      for(int j=0; j<Util::nDefinedParticles; ++j) {
        if(mhTpcRecoEff3D[i][j]) delete mhTpcRecoEff3D[i][j];
      }
    }
  }
  
  bool LoadHistograms(const char* fileName) {
    TFile* file = TFile::Open(fileName);
    if(!file || file->IsZombie()) {
      std::cerr << "Error: Could not open file " << fileName << std::endl;
      return false;
    }
    
    bool success = true;
    for(int i=0; i<Util::nSigns; ++i) {
      for(int j=0; j<Util::nDefinedParticles; ++j) {
        TString histName = "hTPCEffiCD" + histNames[i][j] + "121";
        mhTpcRecoEff3D[i][j] = dynamic_cast<TH3F*>(file->Get(histName));
        
        if(!mhTpcRecoEff3D[i][j]) {
          std::cerr << "Error: Could not find histogram " << histName << std::endl;
          success = false;
          continue;
        }
        
        mhTpcRecoEff3D[i][j]->SetDirectory(0);
        std::cout << "Successfully loaded: " << histName << std::endl;
      }
    }
    
    file->Close();
    return success;
  }
  
  void PrintHistogramInfo() {
    std::cout << "\n=== Histogram Information ===" << std::endl;
    for(int i=0; i<Util::nSigns; ++i) {
      for(int j=0; j<Util::nDefinedParticles; ++j) {
        TH3F* hist = mhTpcRecoEff3D[i][j];
        if(!hist) continue;
        
        std::cout << "Histogram [" << i << "][" << j << "] - " << hist->GetName() << std::endl;
        std::cout << "  Dimensions: " 
                  << hist->GetNbinsX() << " x " 
                  << hist->GetNbinsY() << " x " 
                  << hist->GetNbinsZ() << std::endl;
        
        std::cout << "  X-axis (zVx): [" << hist->GetXaxis()->GetXmin() 
                  << ", " << hist->GetXaxis()->GetXmax() << "]" << std::endl;
        std::cout << "  Y-axis (pt): [" << hist->GetYaxis()->GetXmin() 
                  << ", " << hist->GetYaxis()->GetXmax() << "]" << std::endl;
        std::cout << "  Z-axis (eta): [" << hist->GetZaxis()->GetXmin() 
                  << ", " << hist->GetZaxis()->GetXmax() << "]" << std::endl;
        
        double integral = hist->Integral();
        std::cout << "  Integral: " << integral << std::endl;
        std::cout << "  Mean: "
                  << "X=" << hist->GetMean(1) << ", "
                  << "Y=" << hist->GetMean(2) << ", "
                  << "Z=" << hist->GetMean(3) << std::endl;
        std::cout << std::endl;
      }
    }
}
  
  Double_t tpcRecoEff3D(UInt_t sign, UInt_t pid, Double_t zVx, Double_t eta, Double_t pt) const {
    if(sign >= Util::nSigns || pid >= Util::nDefinedParticles) {
      std::cerr << "Error: Invalid sign or pid value" << std::endl;
      return 0.0;
    }
    
    TH3F *hist = mhTpcRecoEff3D[sign][pid];
    if(!hist) {
      std::cerr << "Error: Histogram not loaded" << std::endl;
      return 0.0;
    }
    
    // Check if values are in range, otherwise move to nearest bin
    double newZVx = (hist->GetXaxis()->GetBinCenter(1) > zVx ? 
                    hist->GetXaxis()->GetBinCenter(1)+1e-3 : 
                    (hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) < zVx ? 
                     hist->GetXaxis()->GetBinCenter(hist->GetNbinsX())-1e-3 : zVx));
                     
    double newEta = (hist->GetZaxis()->GetBinCenter(1) > eta ? 
                    hist->GetZaxis()->GetBinCenter(1)+1e-3 : 
                    (hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ()) < eta ? 
                     hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ())-1e-3 : eta));
                     
    double newPt = (hist->GetYaxis()->GetBinCenter(1) > pt ? 
                   hist->GetYaxis()->GetBinCenter(1)+1e-3 : 
                   (hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()) < pt ? 
                    hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3 : pt));
    
    return hist->Interpolate(newZVx, newPt, newEta);
  }
  
  void TestEfficiency() {
    std::cout << "\n=== Testing Efficiency Interpolation ===" << std::endl;
    
    // Test Point Definition
    struct TestPoint {
      UInt_t sign;
      UInt_t pid;
      Double_t zVx;
      Double_t eta;
      Double_t pt;
      const char* desc;
    };
    
    // Test Point Table
    TestPoint testPoints[] = {
      {Util::MINUS, Util::PION, 0.0, 0.0, 1.0, "Central pion negative"},
      {Util::PLUS, Util::PION, 0.0, 0.0, 1.0, "Central pion positive"},
      {Util::MINUS, Util::KAON, 0.0, 0.0, 1.0, "Central kaon negative"},
      {Util::PLUS, Util::KAON, 0.0, 0.0, 1.0, "Central kaon positive"},
      {Util::MINUS, Util::PROTON, 0.0, 0.0, 1.0, "Central proton negative"},
      {Util::PLUS, Util::PROTON, 0.0, 0.0, 1.0, "Central proton positive"},
      // Additional test points
      {Util::MINUS, Util::PION, 10.0, 0.5, 2.0, "Shifted pion negative"},
      {Util::PLUS, Util::PROTON, -10.0, -0.5, 3.0, "Shifted proton positive"}
    };
    
    const int numTestPoints = sizeof(testPoints) / sizeof(TestPoint);
    
    std::cout << std::left << std::setw(30) << "Description" 
              << std::setw(10) << "Sign" 
              << std::setw(10) << "PID" 
              << std::setw(10) << "zVx" 
              << std::setw(10) << "eta" 
              << std::setw(10) << "pt" 
              << std::setw(15) << "Efficiency" << std::endl;
    std::cout << std::string(85, '-') << std::endl;
    
    for(int i = 0; i < numTestPoints; ++i) {
      const TestPoint& tp = testPoints[i];
      double eff = tpcRecoEff3D(tp.sign, tp.pid, tp.zVx, tp.eta, tp.pt);
      
      std::cout << std::left << std::setw(30) << tp.desc 
                << std::setw(10) << tp.sign 
                << std::setw(10) << tp.pid 
                << std::setw(10) << tp.zVx 
                << std::setw(10) << tp.eta 
                << std::setw(10) << tp.pt 
                << std::setw(15) << eff << std::endl;
    }
  }
  
  void CreateEfficiencyProjections() {
    std::cout << "\n=== Creating Efficiency Projections ===" << std::endl;
    
    // Create canvases for different projections
    TCanvas* canvasPt = new TCanvas("cProjectionsPt", "Efficiency vs pT", 1048, 1500);//800, 1200
    canvasPt->Divide(2, 3);
    
    TCanvas* canvasEta = new TCanvas("cProjectionsEta", "Efficiency vs Eta", 1048, 1500);
    canvasEta->Divide(2, 3);
    
    TCanvas* canvasZ = new TCanvas("cProjectionsZ", "Efficiency vs Z-Vertex", 1048, 1500);//1500, 1048
    canvasZ->Divide(2, 3);
    
    // Particles to display
    const char* particleNames[] = {"#pi^{-}", "#pi^{+}", "K^{-}", "K^{+}", "p^{-}", "p^{+}"};
    int particleConfigs[][2] = {
      {Util::MINUS, Util::PION},
      {Util::PLUS, Util::PION},
      {Util::MINUS, Util::KAON},
      {Util::PLUS, Util::KAON},
      {Util::MINUS, Util::PROTON},
      {Util::PLUS, Util::PROTON}
    };


    // Colors for different particle types
    int colors[] = {kBlue, kRed, kGreen+2, kMagenta, kCyan+2, kOrange+7};
    
    for(int i = 0; i < 6; ++i) {
      int sign = particleConfigs[i][0];
      int pid = particleConfigs[i][1];
      TH3F* hist = mhTpcRecoEff3D[sign][pid];
      
      if(!hist) continue;
      
      //------------------------
      // 1. pT Projections
      //------------------------
      canvasPt->cd(i+1);
      
      // Efficiency projection vs pt
      TH1D* projPt = hist->ProjectionY(Form("projPt_%d", i), 
                                      hist->GetXaxis()->FindBin(-80.0),  // zVx bin
                                      hist->GetXaxis()->FindBin(80.0),   // zVx bin
                                      hist->GetZaxis()->FindBin(-0.9),  // eta bin
                                      hist->GetZaxis()->FindBin(0.9));  // eta bin
      
      projPt->SetTitle(Form("%s Efficiency vs p_{T} (|#eta| < 0.9, |Vtx_{z}| < 80.0)", particleNames[i]));
      projPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      projPt->GetYaxis()->SetTitle("Efficiency");
      projPt->SetLineColor(colors[i]);
      projPt->SetLineWidth(2);
      projPt->SetStats(0);
      
      // Check if we need to normalize
      if(projPt->GetMaximum() > 1.05) {
        std::cout << "Note: pT Histogram for " << particleNames[i] << " appears to need normalization." << std::endl;
        std::cout << "  Maximum value: " << projPt->GetMaximum() << std::endl;
        
        // Calculate number of integrated bins
        int nBinsX = hist->GetXaxis()->FindBin(80.0) - hist->GetXaxis()->FindBin(-80.0) + 1;
        int nBinsZ = hist->GetZaxis()->FindBin(0.9) - hist->GetZaxis()->FindBin(-0.9) + 1;
        int totalIntegratedBins = nBinsX * nBinsZ;

        // Scale to get average instead of sum
        projPt->Scale(1.0 / totalIntegratedBins);
        std::cout << "  Applied scaling to average over " << totalIntegratedBins << " bins." << std::endl;
      }
      
      projPt->GetYaxis()->SetRangeUser(0.0, 1.1);  // Adjusted for proper efficiency range
      projPt->GetXaxis()->SetRangeUser(0.0, 3.0);  // Focus on reasonable pT range
      projPt->SetMarkerColor(colors[i]);
      projPt->Draw("HIST P");
      
      //------------------------
      // 2. Eta Projections
      //------------------------
      canvasEta->cd(i+1);
      
      // Efficiency projection vs eta
      TH1D* projEta = hist->ProjectionZ(Form("projEta_%d", i), 
                                       hist->GetXaxis()->FindBin(-80.0),  // zVx bin
                                       hist->GetXaxis()->FindBin(80.0),   // zVx bin
                                       hist->GetYaxis()->FindBin(0.0),    // pt bin
                                       hist->GetYaxis()->FindBin(2.0));   // pt bin
      
      projEta->SetTitle(Form("%s Efficiency vs #eta (0.0 < p_{T} < 2.0 GeV/c, |Vtx_{z}| < 80.0)", particleNames[i]));
      projEta->GetXaxis()->SetTitle("#eta");
      projEta->GetYaxis()->SetTitle("Efficiency");
      projEta->SetLineColor(colors[i]);
      projEta->SetLineWidth(2);
      projEta->SetStats(0);
      
      // Normalize if needed
      if(projEta->GetMaximum() > 1.05) {
        std::cout << "Note: Eta Histogram for " << particleNames[i] << " appears to need normalization." << std::endl;
        std::cout << "  Maximum value: " << projEta->GetMaximum() << std::endl;  
      
        // Calculate number of integrated bins
        int nBinsX = hist->GetXaxis()->FindBin(80.0) - hist->GetXaxis()->FindBin(-80.0) + 1;
        int nBinsY = hist->GetYaxis()->FindBin(2.0) - hist->GetYaxis()->FindBin(0.0) + 1;
        int totalIntegratedBins = nBinsX * nBinsY;

        // Scale to get average instead of sum
        projEta->Scale(1.0 / totalIntegratedBins);
        std::cout << "  Applied scaling to average over " << totalIntegratedBins << " bins." << std::endl;
      }
      
      projEta->GetYaxis()->SetRangeUser(0.0, 1.1);
      projEta->SetMarkerColor(colors[i]);
      projEta->Draw("HIST P");
      
      //------------------------
      // 3. Z-Vertex Projections
      //------------------------
      canvasZ->cd(i+1);
      
      // Efficiency projection vs z-vertex
      TH1D* projVerZ = hist->ProjectionX(Form("projVerZ_%d", i), 
                                    hist->GetYaxis()->FindBin(0.0),    // pt bin
                                    hist->GetYaxis()->FindBin(2.0),    // pt bin
                                    hist->GetZaxis()->FindBin(-0.9),   // eta bin
                                    hist->GetZaxis()->FindBin(0.9));   // eta bin
      
      projVerZ->SetTitle(Form("%s Efficiency vs Vtx_{z} (0.0 < p_{T} < 2.0 GeV/c, |#eta| < 0.9)", particleNames[i]));
      projVerZ->GetXaxis()->SetTitle("Vtx_{z} (cm)");
      projVerZ->GetYaxis()->SetTitle("Efficiency");
      projVerZ->SetLineColor(colors[i]);
      projVerZ->SetLineWidth(2);
      projVerZ->SetStats(0);
      
      // Normalize if needed
      if(projVerZ->GetMaximum() > 1.05) {
        std::cout << "Note: Z-Vertex Histogram for " << particleNames[i] << " appears to need normalization." << std::endl;
        std::cout << "  Maximum value: " << projVerZ->GetMaximum() << std::endl;
        
        // Calculate number of integrated bins
        int nBinsY = hist->GetYaxis()->FindBin(2.0) - hist->GetYaxis()->FindBin(0.0) + 1;
        int nBinsZ = hist->GetZaxis()->FindBin(0.9) - hist->GetZaxis()->FindBin(-0.9) + 1;
        int totalIntegratedBins = nBinsY * nBinsZ;

        // Scale to get average instead of sum
        projVerZ->Scale(1.0 / totalIntegratedBins);
        std::cout << "  Applied scaling to average over " << totalIntegratedBins << " bins." << std::endl;
      }
      
      projVerZ->GetYaxis()->SetRangeUser(0.0, 1.1);
      projVerZ->GetXaxis()->SetRangeUser(-90.5, 90.5);  // Focus on reasonable z-vertex range
      projVerZ->SetMarkerColor(colors[i]);
      projVerZ->Draw("HIST P");
    }
    
    // Save all canvases
    canvasPt->Update();
    canvasPt->SaveAs("plots/EfficiencyVsPt.png");
    std::cout << "Saved pT projections to EfficiencyVsPt.png" << std::endl;
    
    canvasEta->Update();
    canvasEta->SaveAs("plots/EfficiencyVsEta.png");
    std::cout << "Saved eta projections to EfficiencyVsEta.png" << std::endl;
    
    canvasZ->Update();
    canvasZ->SaveAs("plots/EfficiencyVsZ.png");
    std::cout << "Saved z-vertex projections to EfficiencyVsZ.png" << std::endl;
  }
  
  void DebugEfficiencyValues() {
    std::cout << "\n=== Debugging Efficiency Values ===" << std::endl;
    
    // Select a specific histogram for debugging
    TH3F* hist = mhTpcRecoEff3D[Util::MINUS][Util::PION];
    if(!hist) {
      std::cout << "Histogram not available for debugging." << std::endl;
      return;
    }
    
    double testPt = 1.0;
    double testZVx = 0.0;
    double testEta = 0.0;
    
    // Find bin indices for our test point
    int binX = hist->GetXaxis()->FindBin(testZVx);
    int binY = hist->GetYaxis()->FindBin(testPt);
    int binZ = hist->GetZaxis()->FindBin(testEta);
    
    // Get bin content directly
    double binContent = hist->GetBinContent(binX, binY, binZ);
    
    // Get interpolated value
    double interpValue = tpcRecoEff3D(Util::MINUS, Util::PION, testZVx, testEta, testPt);
    
    std::cout << "Test point: (zVx=" << testZVx << ", pt=" << testPt << ", eta=" << testEta << ")" << std::endl;
    std::cout << "Bin coordinates: (x=" << binX << ", y=" << binY << ", z=" << binZ << ")" << std::endl;
    std::cout << "Direct bin content: " << binContent << std::endl;
    std::cout << "Interpolated value: " << interpValue << std::endl;
    
    if(fabs(binContent) > 1.1 && binContent != 0) {
      std::cout << "\nIMPORTANT: Bin content is > 1.1 which suggests these histograms " << std::endl;
      std::cout << "contain counts rather than normalized efficiency values." << std::endl;
      std::cout << "You may need to divide by the corresponding reference histogram." << std::endl;
      
      // Check several bins to find max value
      double maxValue = 0.0;
      double sumValues = 0.0;
      int nonZeroBins = 0;
      
      for(int x = 1; x <= hist->GetNbinsX(); x++) {
        for(int y = 1; y <= hist->GetNbinsY(); y++) {
          for(int z = 1; z <= hist->GetNbinsZ(); z++) {
            double value = hist->GetBinContent(x, y, z);
            if(value > maxValue) maxValue = value;
            if(value > 0) {
              sumValues += value;
              nonZeroBins++;
            }
          }
        }
      }
      
      double avgValue = nonZeroBins > 0 ? sumValues / nonZeroBins : 0;
      
      std::cout << "Maximum bin value found: " << maxValue << std::endl;
      std::cout << "Average non-zero bin value: " << avgValue << std::endl;
      std::cout << "Suggested normalization factor: " << (maxValue > 0 ? 1.0/maxValue : 1.0) << std::endl;
    }
  }
};

void TPCEff1() {
  std::cout << "Starting TPC efficiency for etaPhiEfficiency_16_01_19_delta015_twoRuns" << std::endl;
  
  EfficiencyChecker checker;
  
  // Try to find the histogram file
  const char* filePaths[] = {
    "../share/etaPhiEfficiency_16_01_19_delta015_twoRuns.root"
  };
  
  bool loaded = false;
  for(const char* filePath : filePaths) {
    std::cout << "Attempting to load from: " << filePath << std::endl;
    if(checker.LoadHistograms(filePath)) {
      loaded = true;
      std::cout << "Successfully loaded from: " << filePath << std::endl;
      break;
    }
  }
  
  if(!loaded) {
    std::cerr << "Failed to load histograms from any path. Exiting." << std::endl;
    return;
  }
  
  // Print information about histograms
  checker.PrintHistogramInfo();
  
  // Test efficiency function
  checker.TestEfficiency();
  
  // Create efficiency projections
  checker.CreateEfficiencyProjections();
  
  // Debug: Compare bin content with interpolation
  checker.DebugEfficiencyValues();
  
  std::cout << "TPC efficiency for etaPhiEfficiency_16_01_19_delta015_twoRuns completed successfully" << std::endl;
  
}
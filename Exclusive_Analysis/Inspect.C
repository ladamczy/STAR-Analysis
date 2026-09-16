void ExtractMassCuts(TH2D* h2D, TString className) {
    if (!h2D || h2D->GetEntries() == 0) {
        std::cerr << "Error: Histogram for " << className << " is null or empty!" << std::endl;
        return;
    }

    std::cout << "\n=================================================" << std::endl;
    std::cout << "   ROBUST MASS RESOLUTION: " << className << std::endl;
    std::cout << "=================================================" << std::endl;

    // 1. Project the 2D histogram into 1D
    TH1D* hLead = h2D->ProjectionX(Form("hLead_%s", className.Data()));
    TH1D* hSub = h2D->ProjectionY(Form("hSub_%s", className.Data()));

    // 2. Define a Signal (gaus) + Background (pol1) function
    // Expand the fit range slightly to capture the background sidebands
    TF1* fitLead = new TF1(Form("fitLead_%s", className.Data()), "gaus(0) + pol1(3)", 0.45, 0.55);
    TF1* fitSub  = new TF1(Form("fitSub_%s", className.Data()), "gaus(0) + pol1(3)", 0.45, 0.55);

    // 3. Set Smart Initial Guesses and STRICT Parameter Limits
    double pdgMass = 0.4976;
    
    // Setup Leading Fit
    fitLead->SetParameters(hLead->GetMaximum(), pdgMass, 0.010, 0, 0); // Amplitude, Mean, Sigma, p0, p1
    fitLead->SetParLimits(1, 0.485, 0.510); // Force Mean to stay near K0s mass
    fitLead->SetParLimits(2, 0.003, 0.025); // Force Sigma to be physically reasonable (3 - 25 MeV)
    
    // Setup Subleading Fit
    fitSub->SetParameters(hSub->GetMaximum(), pdgMass, 0.010, 0, 0);
    fitSub->SetParLimits(1, 0.485, 0.510);
    fitSub->SetParLimits(2, 0.003, 0.025);

    // 4. Perform the fits (Remove 'Q' to see minimizer warnings if it fails again)
    std::cout << "-> Fitting Leading K0s..." << std::endl;
    hLead->Fit(fitLead, "R0"); 
    
    std::cout << "-> Fitting Subleading K0s..." << std::endl;
    hSub->Fit(fitSub, "R0");

    // 5. Extract purely the Gaussian parameters (Params 0, 1, 2)
    double muLead  = fitLead->GetParameter(1);
    double sigLead = fitLead->GetParameter(2);
    
    double muSub  = fitSub->GetParameter(1);
    double sigSub = fitSub->GetParameter(2);

    // 6. Output the results
    std::cout << "\n--- LEADING K0s ---" << std::endl;
    std::cout << "  Mean   (mu): " << muLead << " GeV/c^2" << std::endl;
    std::cout << "  Width (sig): " << sigLead << " GeV/c^2" << std::endl;
    std::cout << "  3-Sigma Cut: [ " << muLead - 3*sigLead << " , " << muLead + 3*sigLead << " ]" << std::endl;

    std::cout << "\n--- SUBLEADING K0s ---" << std::endl;
    std::cout << "  Mean   (mu): " << muSub << " GeV/c^2" << std::endl;
    std::cout << "  Width (sig): " << sigSub << " GeV/c^2" << std::endl;
    std::cout << "  3-Sigma Cut: [ " << muSub - 3*sigSub << " , " << muSub + 3*sigSub << " ]" << std::endl;
    std::cout << "=================================================\n" << std::endl;
    
    // Optional: Save the projected histograms with their fits so you can physically look at them
    // hLead->Write();
    // hSub->Write();
}
void ProjectionNFit(TH2D* h2, double massLow, double massHigh, double fitmin, double fitmax, const char* suffix, double globalMean = -1.0, double globalSigma = -1.0) {
    // ROOT needs to know which BINS correspond to these axes values
    int binXlow  = h2->GetXaxis()->FindBin(massLow);
    int binXhigh = h2->GetXaxis()->FindBin(massHigh);
    int binYlow  = h2->GetYaxis()->FindBin(massLow);
    int binYhigh = h2->GetYaxis()->FindBin(massHigh);

    // 4. Print statements to understand the data
    //cout << "\n========================================================" << endl;
    //cout << " 🕵️ INSPECTING: " << h2->GetName() << endl;
    //cout << "========================================================" << endl;
    //cout << "Total entries (all cuts applied EXCEPT mass window) : " << h2->GetEntries() << endl;
    //cout << "X-axis (Leading) bin range for [" << massLow << ", " << massHigh << "] GeV : " << binXlow << " to " << binXhigh << endl;
    //cout << "Y-axis (Sublead) bin range for [" << massLow << ", " << massHigh << "] GeV : " << binYlow << " to " << binYhigh << endl;
    
    // 5. Create the projections
    
    // Scenario A: Cut on Leading mass (X-axis), plot Subleading mass (Y-axis)
    TH1D* hSubleading = h2->ProjectionY("hSubleading", binXlow, binXhigh);
    hSubleading->SetTitle("Subleading Mass (Leading in Mass Window);m_{sublead} [GeV];Events");
    hSubleading->SetLineColor(kBlue);
    hSubleading->SetLineWidth(2);
    
    // Scenario B: Cut on Subleading mass (Y-axis), plot Leading mass (X-axis)
    TH1D* hLeading = h2->ProjectionX("hLeading", binYlow, binYhigh);
    hLeading->SetTitle("Leading Mass (Subleading in Mass Window);m_{lead} [GeV];Events");
    hLeading->SetLineColor(kRed);
    hLeading->SetLineWidth(2);

    //cout << "\n--- Projection Results ---" << endl;
    //cout << "Scenario A: Entries where Leading is in mass window    = " << hSubleading->GetEntries() << endl;
    //cout << "Scenario B: Entries where Subleading is in mass window = " << hLeading->GetEntries() << endl;
    ///cout << "========================================================\n" << endl;

    // 6. Draw everything
    TCanvas* c1 = new TCanvas("c1", "N-1 Mass Inspection", 1500, 500);
    c1->Divide(3, 1);

    // Panel 1: The 2D Histogram
    c1->cd(1);
    h2->SetTitle("2D N-1 Mass Distribution;m_{lead} [GeV];m_{sublead} [GeV]");
    h2->GetZaxis()->SetRangeUser(0, h2->GetMaximum());
    h2->Draw("COLZ");
    
    // Draw a box to visually indicate the double-mass signal region
    TBox *box = new TBox(massLow, massLow, massHigh, massHigh);
    box->SetFillStyle(0); // Transparent fill
    box->SetLineColor(kRed);
    box->SetLineWidth(2);
    box->Draw("SAME");

    // Panel 2: Scenario A (Subleading mass)
    c1->cd(2);
    if (hSubleading->GetMaximum() < 9) {
        hSubleading->GetYaxis()->SetRangeUser(0, 10);
    } else {
        hSubleading->GetYaxis()->SetRangeUser(0, 18);
    } 
    hSubleading->Draw("HIST E");

    // Panel 3: Scenario B (Leading mass)
    c1->cd(3);
    if (hLeading->GetMaximum() < 9) {
        hLeading->GetYaxis()->SetRangeUser(0, 10);
    } else {
        hLeading->GetYaxis()->SetRangeUser(0, 18);
    } 
    hLeading->Draw("HIST E");

    //c1->SaveAs(Form("plots/N1_Mass_Inspection_%.2f-%.2f_%s.png", massLow, massHigh, suffix));


    hSubleading->Rebin(5);
    // 3. Set color and style for the histogram again
    hSubleading->SetTitle(Form("Subleading Mass (Lead cut: %.2f - %.2f GeV);m_{sublead} [GeV];Events", massLow, massHigh));
    hSubleading->SetLineColor(kBlack);
    hSubleading->SetMarkerStyle(20);

    // 4. Fit Signal (Gaussian) + Background (Linear)
    //TF1 *fitFunc = new TF1("fitFunc", "gaus(0) + pol1(3)", fitmin, fitmax);
    TF1 *fitFunc = new TF1("fitFunc", "gaus(0) + pol0(3)", fitmin, fitmax); //has 3 par
    
    // Set robust initial guesses for the fit
    fitFunc->SetParameters(hSubleading->GetMaximum(), 0.5, 0.005, 2, 0);
    fitFunc->SetLineColor(kRed);

    fitFunc->SetParLimits(1, 0.49, 0.51);    // 
    // ====== GLOBAL FIT LOGIC APPLIED HERE ======
    if (globalMean > 0 && globalSigma > 0) {
        cout << " -> Locking Mean to Global Value: " << globalMean << " GeV" << endl;
        cout << " -> Locking Sigma to Global Value: " << globalSigma << " GeV" << endl;
        //fitFunc->FixParameter(1, globalMean);
        fitFunc->FixParameter(2, globalSigma);
    } else {
        fitFunc->SetParLimits(1, 0.49, 0.51);    // Mean must be near K0s mass
        fitFunc->SetParLimits(2, 0.002, 0.015);  // Sigma must be realistic
    }
    
    //fitFunc->SetParLimits(2, 0.002, 0.015);  // 
    fitFunc->SetParLimits(3, 0.0, hSubleading->GetMaximum());

    // ====== ADDED THE 'I' OPTION FOR BIN INTEGRAL ======
    // R = use function range, Q = quiet, I = Use integral of function in bin
    hSubleading->Fit("fitFunc", "R Q I L"); 

    // 5. Extract the individual Signal and Background functions
    TF1 *sigFunc = new TF1("sigFunc", "gaus(0)", fitmin, fitmax);
    sigFunc->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1), fitFunc->GetParameter(2));
    sigFunc->SetLineColor(kBlue);
    sigFunc->SetLineStyle(1);

    //TF1 *bkgFunc = new TF1("bkgFunc", "pol1(0)", fitmin, fitmax);
    //bkgFunc->SetParameters(fitFunc->GetParameter(3), fitFunc->GetParameter(4));
    TF1 *bkgFunc = new TF1("bkgFunc", "pol0(0)", fitmin, fitmax);
    bkgFunc->SetParameter(0, fitFunc->GetParameter(3));
    bkgFunc->SetLineColor(kGreen+2);
    bkgFunc->SetLineStyle(2);

    // 6. Integrate functions inside our chosen mass window to get exact yields
    double binWidth = hSubleading->GetBinWidth(1);
    double S = sigFunc->Integral(massLow, massHigh) / binWidth;
    double B = bkgFunc->Integral(massLow, massHigh) / binWidth;

    double fitMean = fitFunc->GetParameter(1);
    double fitSigma = fitFunc->GetParameter(2);

    cout << Form(" Extracted Global Mean  = %.4f GeV", fitMean) << endl;
    cout << Form(" Extracted Global Sigma = %.4f GeV", fitSigma) << endl;
        
    double SoverB = S / B;
    double significance = S / sqrt(S + B);

    // 7. Print results to terminal
    cout << "\n========================================================" << endl;
    cout << " 📊 FIT RESULTS FOR WINDOW: [" << massLow << " - " << massHigh << "] GeV" << endl;
    cout << "========================================================" << endl;
    cout << Form(" Signal (S)     = %.1f events", S) << endl;
    cout << Form(" Background (B) = %.1f events", B) << endl;
    cout << Form(" S / B Ratio    = %.2f", SoverB) << endl;
    cout << Form(" Significance   = %.2f sigma", significance) << endl;
    cout << "========================================================\n" << endl;

    // 8. Draw it beautifully for presentation
    TCanvas* c2 = new TCanvas("c2", "Signal vs Background", 800, 600);
    hSubleading->GetYaxis()->SetRangeUser(0, hSubleading->GetMaximum() * 3.8);
    hSubleading->Draw("E1");
    fitFunc->Draw("SAME");
    sigFunc->Draw("SAME");
    bkgFunc->Draw("SAME");

    // Add a text box with the stats
    TPaveText *pt = new TPaveText(0.65, 0.65, 0.88, 0.88, "NDC");
    pt->SetFillColor(0);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);
    pt->AddText(Form("Mass Cut: [%.2f, %.2f]", massLow, massHigh));
    pt->AddText(Form("S = %.1f", S));
    pt->AddText(Form("B = %.1f", B));
    pt->AddText(Form("S/B = %.2f", SoverB));
    pt->AddText(Form("S/#sqrt{S+B} = %.2f", significance));
    pt->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.15, 0.78, Form("A = %.4f", fitFunc->GetParameter(0)));
    if (globalMean > 0) {
        latex.DrawLatex(0.15, 0.72, Form("fixed #mu = %.4f GeV/c^{2}", fitFunc->GetParameter(1)));
    }
    if (globalSigma > 0) {    
        latex.DrawLatex(0.15, 0.66, Form("fixed #sigma = %.4f GeV/c^{2}", fitFunc->GetParameter(2)));
    }  
    if (globalMean < 0 && globalSigma < 0) {
        latex.DrawLatex(0.15, 0.72, Form("#mu = %.4f GeV/c^{2}", fitFunc->GetParameter(1)));
        latex.DrawLatex(0.15, 0.66, Form("#sigma = %.4f GeV/c^{2}", fitFunc->GetParameter(2)));
    }
    c2->SaveAs(Form("plots/N1_Mass_Fit_Results_fixedSigma_%.2f-%.2f_%s.png", massLow, massHigh, suffix));
}

void PlotFinalMass(TH1D* hMass, const char* suffix) {
    TCanvas* c = new TCanvas("c_mass", "Mass", 800, 600);
    hMass->SetMarkerStyle(20);
    hMass->SetMarkerColor(kBlack);
    hMass->SetLineColor(kBlack);
    hMass->GetXaxis()->SetRangeUser(0.9, 3.0);
    hMass->GetXaxis()->SetTitle("m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]");
    hMass->GetYaxis()->SetTitle("events");
    //hMass->GetYaxis()->SetTitleOffset(1.1);
    //hMass->GetXaxis()->SetTitleOffset(1.1);
    hMass->GetXaxis()->SetTitleSize(0.06); //axis title size
    hMass->GetYaxis()->SetTitleSize(0.06);

    hMass->GetXaxis()->SetLabelSize(0.05);
    hMass->GetYaxis()->SetLabelSize(0.05);
    
    hMass->Draw("PE");

    TF1* fitFunc = new TF1("fitFunc", "gaus(0) + pol1(3)", 1.0, 2.5);
    fitFunc->SetParameters(15, 1.68, 0.08, 10, -2); 
    fitFunc->SetLineColor(kRed);
    hMass->Fit(fitFunc, "R"); 

    TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(hMass, "Data", "pe");
    leg->AddEntry(fitFunc, "Fit: gaus + pol1", "l");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.55, 0.68, Form("A = %.2f", fitFunc->GetParameter(0)));
    latex.DrawLatex(0.55, 0.62, Form("#mu = %.2f GeV/c^{2}", fitFunc->GetParameter(1)));
    latex.DrawLatex(0.55, 0.56, Form("#sigma = %.2f GeV/c^{2}", fitFunc->GetParameter(2)));

    c->SaveAs(Form("plots/final_K0K0_InvariantMass_%s.pdf", suffix));
    delete c;
}


void CheckScaleFactor(TH1D* h1, TH1D* h2, TH1D* h3, double sidebandLow, double sidebandHigh, const char* suffix) {
    
    // 2. Define the pure combinatorial sideband range
    int binLow = h1->FindBin(sidebandLow);
    int binHigh = h1->FindBin(sidebandHigh);

    // =========================================================
    // METHOD A: Sideband Integral Matching
    // =========================================================
    double rawSidebandInt = h1->Integral(binLow, binHigh);
    double mixSidebandInt = h2->Integral(binLow, binHigh);
    double scaleFactorA = 0.0;
    if (mixSidebandInt > 0) scaleFactorA = rawSidebandInt / mixSidebandInt;

    // =========================================================
    // METHOD B: Like-Sign Yield Matching
    // =========================================================
    // We integrate the ENTIRE mass range (0.3 to 0.7) for both LS and Mixed
    double lsTotalInt = h3->Integral();
    double mixTotalInt = h2->Integral();
    double scaleFactorB = 0.0;
    if (mixTotalInt > 0) scaleFactorB = lsTotalInt / mixTotalInt;


    char* suffixStr = const_cast<char*>(suffix);

    // 3. Print the comparison
    std::cout << "\n" << std::endl;
    cout << "BACKGROUND SCALING FOR "<< suffixStr << endl;
    std::cout << "Method A (Sideband " << sidebandLow << "-" << sidebandHigh << "): " << scaleFactorA << std::endl;
    std::cout << "Method B (Like-Sign Yield) F_scale:    " << scaleFactorB << std::endl;
    
    // To apply Method B to your plot:
    // hMix2->Scale(scaleFactorB);
}

void DrawArmenteros(TH2D* h1, TH2D* h2, TH2D* h3, const char* suffix) {

    TH2D* HistArmenteros_combined = (TH2D*)h1->Clone("HistArmenterosBefore_combined");
    HistArmenteros_combined->Add(h2);
    HistArmenteros_combined->Add(h3);

    // 2. Draw the 2D histogram
    TCanvas* c1 = new TCanvas("c1", "Armenteros-Podolansky", 800, 600);
    HistArmenteros_combined->Draw("colz"); // "colz" is best for seeing the population density

    // 3. Define the theoretical ellipse kinematics (assuming beta ~ 1)
    double K0_mass = 0.4976; 
    double pi_mass = 0.1396;
    double q = sqrt((K0_mass * K0_mass / 4.0) - (pi_mass * pi_mass)); // ~0.206
    double alpha_max = 2.0 * q / K0_mass;                             // ~0.828

    // 4. Create and format the TEllipse
    // TEllipse(x_center, y_center, x_radius, y_radius, phi_min, phi_max)
    TEllipse* perfectEllipse = new TEllipse(0.0, 0.0, alpha_max, q, 0, 180);
    perfectEllipse->SetLineColor(kRed);
    perfectEllipse->SetLineWidth(2);
    perfectEllipse->SetLineStyle(2); // Dashed line looks cleaner over data
    perfectEllipse->SetFillStyle(0); // Make it hollow    
    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);
    //latex.DrawLatex(0.15, 0.85, "extracted only from leading K0s");

    // 5. Draw the ellipse on top of the histogram
    perfectEllipse->Draw("same");

    c1->SaveAs(Form("plots/Armenteros_Podolansky_2D_%s.png", suffix));
}

void Draw3(TH1D* h1, TH1D* h2, TH1D* h3, double ymin, double ymax, const char* suffix) {
    TCanvas* c = new TCanvas("c", "Mass Distributions", 800, 600);
    h1->SetLineColor(kBlue);
    h2->SetLineColor(kRed);
    h3->SetLineColor(kGreen+2);

    h1->SetLineWidth(2);
    h2->SetLineWidth(2);
    h3->SetLineWidth(2);
    h3->GetYaxis()->SetRangeUser(ymin, ymax);
    h3->Draw("HIST");
    h2->Draw("HIST SAME");
    h1->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(h1, "2-TOF Class", "l");
    leg->AddEntry(h2, "3-TOF Class", "l");
    leg->AddEntry(h3, "4-TOF Class", "l");
    leg->Draw();

    c->SaveAs(Form("plots/Mass_Distributions_%s.png", suffix));

}


void DrawSignalExtraction(TH1D* hRaw, TH1D* hMix, double scaleFactor, TString className) {
    // 1. Clone histograms so we don't permanently modify the originals
    TH1D* hRawDraw = (TH1D*)hRaw->Clone(Form("hRawDraw_%s", className.Data()));
    TH1D* hMixScaled = (TH1D*)hMix->Clone(Form("hMixScaled_%s", className.Data()));
    
    // 2. Apply the scale factor to the mixed background
    hMixScaled->Scale(scaleFactor);

    // 3. Create the Signal histogram by subtracting background from raw
    TH1D* hSignal = (TH1D*)hRaw->Clone(Form("hSignal_%s", className.Data()));
    hSignal->Add(hMixScaled, -1.0); 

    // 4. Formatting
    hRawDraw->SetMarkerStyle(20);
    hRawDraw->SetMarkerColor(kBlack);
    hRawDraw->SetLineColor(kBlack);
    hRawDraw->SetTitle(Form("%s: Raw Data vs Scaled Mixed Background", className.Data()));
    hRawDraw->GetYaxis()->SetTitle("Events");
    
    hMixScaled->SetLineColor(kRed);
    hMixScaled->SetLineWidth(2);
    hMixScaled->SetFillColor(kRed);
    hMixScaled->SetFillStyle(3004); // Nice hatched fill
    
    hSignal->SetMarkerStyle(21);
    hSignal->SetMarkerColor(kBlue);
    hSignal->SetLineColor(kBlue);
    hSignal->SetTitle(""); // Let the top pad handle the title
    hSignal->GetXaxis()->SetTitle("m_{#pi#pi} [GeV/c^{2}]");
    hSignal->GetYaxis()->SetTitle("Pure Signal");
    hSignal->GetYaxis()->SetTitleSize(0.05);

    // 5. Canvas Setup (2 Panels)
    TCanvas* c1 = new TCanvas(Form("c_Signal_%s", className.Data()), Form("Signal Extraction %s", className.Data()), 800, 800);
    c1->Divide(1, 2);
    
    // --- TOP PAD: OVERLAY ---
    c1->cd(1);
    gPad->SetBottomMargin(0.02); // Minimize gap between pads
    // Adjust Y-axis maximum to fit both distributions and the legend comfortably
    double maxY = hRawDraw->GetMaximum();
    hRawDraw->SetMaximum(maxY * 1.4);
    hRawDraw->Draw("E1");             // Raw data with error bars
    hMixScaled->Draw("HIST SAME");    // Background as filled histogram
    
    TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hRawDraw, "Raw Data", "lep");
    leg->AddEntry(hMixScaled, Form("Mixed Bkg (F_{scale} = %.3f)", scaleFactor), "f");
    leg->Draw();

    // --- BOTTOM PAD: SUBTRACTED SIGNAL ---
    c1->cd(2);
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.25);
    
    // Add a dotted line at y=0 to guide the eye
    hSignal->Draw("E1"); 
    TLine* lineZero = new TLine(hSignal->GetXaxis()->GetXmin(), 0, hSignal->GetXaxis()->GetXmax(), 0);
    lineZero->SetLineStyle(2);
    lineZero->Draw("SAME");

    // 6. Save the canvas
    c1->SaveAs(Form("plots/SignalExtraction_%s.png", className.Data()));
}

void DrawHistsFromFile(const char* filename, const char* suffix)
{
    // This is majorly done for without mass window cut except  for HistArmenterosCore

    TFile* file = TFile::Open(filename, "READ");

    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }

    TH2D* HistArmenterosCore_2Tof = (TH2D*)file->Get("histArmenterosCore2_nominal");
    TH2D* HistArmenterosCore_3Tof = (TH2D*)file->Get("histArmenterosCore3_nominal");
    TH2D* HistArmenterosCore_4Tof = (TH2D*)file->Get("histArmenterosCore4_nominal");

    //DrawArmenteros(HistArmenterosCore_2Tof, HistArmenterosCore_3Tof, HistArmenterosCore_4Tof, "Core");

    TH2D* HistArmenterosFringe_2Tof = (TH2D*)file->Get("histArmenterosFringe2_nominal");
    TH2D* HistArmenterosFringe_3Tof = (TH2D*)file->Get("histArmenterosFringe3_nominal");
    TH2D* HistArmenterosFringe_4Tof = (TH2D*)file->Get("histArmenterosFringe4_nominal");

    //DrawArmenteros(HistArmenterosFringe_2Tof, HistArmenterosFringe_3Tof, HistArmenterosFringe_4Tof, "Fringe");

    TH2D* HistArmenterosBefore_2Tof = (TH2D*)file->Get("histArmenterosBefore2_nominal");
    TH2D* HistArmenterosBefore_3Tof = (TH2D*)file->Get("histArmenterosBefore3_nominal");
    TH2D* HistArmenterosBefore_4Tof = (TH2D*)file->Get("histArmenterosBefore4_nominal");

    //DrawArmenteros(HistArmenterosBefore_2Tof, HistArmenterosBefore_3Tof, HistArmenterosBefore_4Tof, "Before");

    TH2D* HistArmenterosMixed_2Tof = (TH2D*)file->Get("histArmenterosMixed2_nominal");
    TH2D* HistArmenterosMixed_3Tof = (TH2D*)file->Get("histArmenterosMixed3_nominal");
    TH2D* HistArmenterosMixed_4Tof = (TH2D*)file->Get("histArmenterosMixed4_nominal");

    DrawArmenteros(HistArmenterosMixed_2Tof, HistArmenterosMixed_3Tof, HistArmenterosMixed_4Tof, "Mixed_combined");

    TH1D* HistMassMixed_2Tof = (TH1D*)file->Get("histMassMixed2_nominal");
    TH1D* HistMassMixed_3Tof = (TH1D*)file->Get("histMassMixed3_nominal");
    TH1D* HistMassMixed_4Tof = (TH1D*)file->Get("histMassMixed4_nominal");

    TH1D* HistMassRaw_2Tof = (TH1D*)file->Get("histMassRaw2_nominal");
    TH1D* HistMassRaw_3Tof = (TH1D*)file->Get("histMassRaw3_nominal");
    TH1D* HistMassRaw_4Tof = (TH1D*)file->Get("histMassRaw4_nominal");

    TH1D* HistMassLS_2Tof = (TH1D*)file->Get("histMassLS2_nominal");
    TH1D* HistMassLS_3Tof = (TH1D*)file->Get("histMassLS3_nominal");
    TH1D* HistMassLS_4Tof = (TH1D*)file->Get("histMassLS4_nominal");

    //dont draw statistics box
    gStyle->SetOptStat(0);

    Draw3(HistMassRaw_2Tof, HistMassRaw_3Tof, HistMassRaw_4Tof, 0, 28, "Mass_Raw_Comparison");
    Draw3(HistMassMixed_2Tof, HistMassMixed_3Tof, HistMassMixed_4Tof, 0, 25, "Mass_Mixed_Comparison"); 
    Draw3(HistMassLS_2Tof, HistMassLS_3Tof, HistMassLS_4Tof, 0, 15, "Mass_LS_Comparison");

    CheckScaleFactor(HistMassRaw_4Tof, HistMassMixed_4Tof, HistMassLS_4Tof, 0.55, 0.70, "4-TOF");
    CheckScaleFactor(HistMassRaw_3Tof, HistMassMixed_3Tof, HistMassLS_3Tof, 0.55, 0.70, "3-TOF");
    CheckScaleFactor(HistMassRaw_2Tof, HistMassMixed_2Tof, HistMassLS_2Tof, 0.52, 0.70, "2-TOF");

    // =========================================================
    // 4-TOF: METHOD A (Sideband Matching 0.55 - 0.70)
    // =========================================================
    int binLow = HistMassRaw_4Tof->FindBin(0.55);
    int binHigh = HistMassRaw_4Tof->FindBin(0.70);
    
    double rawInt4 = HistMassRaw_4Tof->Integral(binLow, binHigh);
    double mixInt4 = HistMassMixed_4Tof->Integral(binLow, binHigh);
    double scaleFactor4_MethodA = (mixInt4 > 0) ? (rawInt4 / mixInt4) : 0;
    
    DrawSignalExtraction(HistMassRaw_4Tof, HistMassMixed_4Tof, scaleFactor4_MethodA, "4-TOF");

    // =========================================================
    // 3-TOF: METHOD B (Like-Sign Yield Matching)
    // =========================================================
    double lsInt3 = HistMassLS_3Tof->Integral();
    double mixIntTotal3 = HistMassMixed_3Tof->Integral();
    double scaleFactor3_MethodB = (mixIntTotal3 > 0) ? (lsInt3 / mixIntTotal3) : 0;
    
    DrawSignalExtraction(HistMassRaw_3Tof, HistMassMixed_3Tof, scaleFactor3_MethodB, "3-TOF");

    // =========================================================
    // 2-TOF: METHOD B (Like-Sign Yield Matching)
    // =========================================================
    double lsInt2 = HistMassLS_2Tof->Integral();
    double mixIntTotal2 = HistMassMixed_2Tof->Integral();
    double scaleFactor2_MethodB = (mixIntTotal2 > 0) ? (lsInt2 / mixIntTotal2) : 0;
    
    DrawSignalExtraction(HistMassRaw_2Tof, HistMassMixed_2Tof, scaleFactor2_MethodB, "2-TOF");

}
void Inspect() {

    // 1. Open the file
    string DataFile1 = "MixedEvMW4753.root"; //wide
    string DataFile2  = "MixedEvMW4852.root"; //narrow
   
    TFile* fwide    = TFile::Open(DataFile1.c_str(), "READ");
    TFile* fnarrow = TFile::Open(DataFile2.c_str(), "READ");

    // 2. Fetch the newly created 2D N-1 histogram.
    TH2D* histInvMass_4ToF = (TH2D*)fwide->Get("histInvMassPiPi2DN14_nominal");
    TH2D* histInvMass_3ToF = (TH2D*)fwide->Get("histInvMassPiPi2DN13_nominal");
    TH2D* histInvMass_2ToF = (TH2D*)fwide->Get("histInvMassPiPi2DN12_nominal");

    //Add 3 hists
    TH2D* histInvMass_combined = (TH2D*)histInvMass_4ToF->Clone("histInvMass_combined");
    histInvMass_combined->Add(histInvMass_3ToF);
    histInvMass_combined->Add(histInvMass_2ToF);

    /*ProjectionNFit(histInvMass_4ToF, 0.48, 0.52, 0.46, 0.54, "4ToF");
    ProjectionNFit(histInvMass_4ToF, 0.47, 0.53, 0.46, 0.54, "4ToF");
    ProjectionNFit(histInvMass_3ToF, 0.48, 0.52, 0.46, 0.54, "3ToF");
    ProjectionNFit(histInvMass_3ToF, 0.47, 0.53, 0.46, 0.54, "3ToF");
    ProjectionNFit(histInvMass_2ToF, 0.48, 0.52, 0.44, 0.56, "2ToF");
    ProjectionNFit(histInvMass_2ToF, 0.47, 0.53, 0.44, 0.56, "2ToF");*/

    //ProjectionNFit(histInvMass_combined, 0.48, 0.52, 0.44, 0.56, "combined_narrow");
   
    double global_mean_narrow = 0.4946;  // extracted value
    double global_sigma_narrow = 0.0114; // extracted value

    // Run the individual classes with locked kinematics
    //ProjectionNFit(histInvMass_4ToF, 0.48, 0.52, 0.46, 0.54, "4ToF_Locked", global_mean_narrow, global_sigma_narrow);
    //ProjectionNFit(histInvMass_3ToF, 0.48, 0.52, 0.46, 0.54, "3ToF_Locked", global_mean_narrow, global_sigma_narrow);
    //ProjectionNFit(histInvMass_2ToF, 0.48, 0.52, 0.44, 0.56, "2ToF_Locked", global_mean_narrow, global_sigma_narrow);
    

    //ProjectionNFit(histInvMass_combined, 0.47, 0.53, 0.44, 0.56, "combined_wide");
 
    double global_mean_wide = 0.4963;  // extracted value
    double global_sigma_wide = 0.0121; // extracted value

    // Run the individual classes with locked kinematics
    //ProjectionNFit(histInvMass_4ToF, 0.47, 0.53, 0.46, 0.54, "4ToF_Locked", global_mean_wide, global_sigma_wide);
    //ProjectionNFit(histInvMass_3ToF, 0.47, 0.53, 0.46, 0.54, "3ToF_Locked", global_mean_wide, global_sigma_wide);
    //ProjectionNFit(histInvMass_2ToF, 0.47, 0.53, 0.44, 0.56, "2ToF_Locked", global_mean_wide, global_sigma_wide);
    

    TH1D* HistMassK0K0_nominal_wide = (TH1D*)fwide->Get("HistMassK0K0_nominal");
    //PlotFinalMass(HistMassK0K0_nominal_wide, "wide");
    TH1D* HistMassK0K0_nominal_narrow = (TH1D*)fnarrow->Get("HistMassK0K0_nominal");
    //PlotFinalMass(HistMassK0K0_nominal_narrow, "narrow");

    TH1D* HistMassK0K0_nominal_difference = (TH1D*)HistMassK0K0_nominal_wide->Clone("HistMassK0K0_nominal_difference");
    HistMassK0K0_nominal_difference->Add(HistMassK0K0_nominal_narrow, -1.0);
    //PlotFinalMass(HistMassK0K0_nominal_difference, "difference");

    ExtractMassCuts(histInvMass_4ToF, "4-TOF Class");
    ExtractMassCuts(histInvMass_3ToF, "3-TOF Class");
    ExtractMassCuts(histInvMass_2ToF, "2-TOF Class");
    ExtractMassCuts(histInvMass_combined, "Combined (Inclusive)");

   
    DrawHistsFromFile("MixedEvMW4852.root", "narrow");
    
}
void Draw dEdxWithBichselCurves() {
    // 1. Otwórz plik i pobierz histogram dE/dx vs p (np. z PicoDst)
    TH2D* h2_dedx = (TH2D*)gDirectory->Get("h2_dedx_vs_p");
    
    TCanvas *c1 = new TCanvas("c1", "STAR TPC dE/dx PID", 800, 600);
    c1->SetLogx(); // Oś pędu zazwyczaj rysuje się w logu
    
    // Rysujemy tło z danymi eksperymentalnymi
    h2_dedx->Draw("colz");
    
    // 2. Generujemy linie Bichsela dla różnych cząstek
    double m_pion   = 0.13957;
    double m_kaon   = 0.49368;
    double m_proton = 0.93827;
    
    // Tworzymy markery (Czerwone kropki dla pionów, Zielone dla kaonów, Niebieskie dla protonów)
    TPolyMarker* pm_pion   = CreateBichselCurve(m_pion, kRed, 20);
    TPolyMarker* pm_kaon   = CreateBichselCurve(m_kaon, kGreen+2, 20);
    TPolyMarker* pm_proton = CreateBichselCurve(m_proton, kBlue, 20);
    
    // 3. Nakładamy linie Bichsela na wykres
    pm_pion->Draw("same");
    pm_kaon->Draw("same");
    pm_proton->Draw("same");
}

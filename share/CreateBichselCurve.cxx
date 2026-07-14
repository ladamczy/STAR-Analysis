#include "TPolyMarker.h"
#include "TMath.h"
#include "StBichsel.h"
#include <vector>

// Funkcja generująca obiekt TPolyMarker dla konkretnej cząstki (np. pion, proton)
TPolyMarker* CreateBichselCurve(double massHypothesis, int markerColor, int markerStyle) {
    
    // Instancja modelu Bichsela dla STAR
    StBichsel* mBichsel = StBichsel::Instance();
    
    // Definiujemy zakres pędu (Momentum p w GeV/c) dla osi X wykresu
    double p_min = 0.1;  // 100 MeV/c
    double p_max = 10.0; // 10 GeV/c
    int nPoints = 200;   // Gęstość markerów na wykresie
    
    // Tablice dynamiczne dla współrzędnych TPolyMarker
    std::vector<double> x_momentum;
    std::vector<double> y_dedx;
    
    // Założona stała długość toru próbkowania (dx) w TPC, np. 25 cm
    double log10_dx = TMath::Log10(25.0); 

    double step = (TMath::Log10(p_max) - TMath::Log10(p_min)) / nPoints;

    for (int i = 0; i < nPoints; ++i) {
        // Generujemy punkty w skali logarytmicznej pędu dla ładniejszego rozkładu
        double log10_p = TMath::Log10(p_min) + i * step;
        double p = TMath::Power(10, log10_p);
        
        // Obliczamy beta*gamma (bg) dla danej masy
        double bg = p / massHypothesis;
        double log10_bg = TMath::Log10(bg);
        
        // Pobieramy teoretyczną stratę energii z modelu Bichsela
        // GetMeandEdx zwraca log10(dE/dx w keV/cm)
        double log10_dedx_theo = mBichsel->GetMeandEdx(log10_bg, log10_dx);
        double dedx_theo = TMath::Power(10, log10_dedx_theo); // Konwersja na keV/cm
        
        x_momentum.push_back(p);
        y_dedx.push_back(dedx_theo);
    }
    
    // Inicjalizacja obiektu TPolyMarker z biblioteki ROOT
    TPolyMarker* pm = new TPolyMarker(nPoints, &x_momentum[0], &y_dedx[0]);
    
    // Ustawianie atrybutów graficznych za pomocą metod TAttMarker
    pm->SetMarkerColor(markerColor);
    pm->SetMarkerStyle(markerStyle);
    pm->SetMarkerSize(0.6);
    
    return pm;
}

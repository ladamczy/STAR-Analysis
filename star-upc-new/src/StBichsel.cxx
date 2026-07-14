#include "StBichsel.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TMath.h"
#include <iostream>
#include <fstream>

StBichsel* StBichsel::fgInstance = nullptr;

// Implementacja Singletonu
StBichsel* StBichsel::Instance() {
  if (!fgInstance) {
    fgInstance = new StBichsel();
  }
  return fgInstance;
}

// Konstruktor klasy
StBichsel::StBichsel(const char *tabName) : TNamed(tabName, "Bichsel Model for STAR TPC") {
  m_Nbg = 0; m_bg = nullptr;
  m_Ndx = 0; m_dx = nullptr;
  m_table = nullptr;

  // W STAR parametryzacja Bichsela znajduje się w ścieżce $STAR/StarDb/Amd/Bichsel
  // Poniższy kod symuluje proces mapowania bazy danych struktur dE/dx
  const char *path = gSystem->ExpandPathName("$STAR/StarDb/Amd/Bichsel/BichselWin.dat");
  std::ifstream in(path);
  
  if (!in.good()) {
    std::cerr << "!!! StBichsel::Error: Nie znaleziono pliku tablicowania Bichsela !!!" << std::endl;
    return;
  }

  // Wczytywanie wymiarów macierzy log10(bg) oraz log10(dx)
  in >> m_Nbg >> m_Ndx;
  
  m_bg = new double[m_Nbg];
  m_dx = new double[m_Ndx];
  m_table = new double[m_Nbg * m_Ndx];

  for (int i = 0; i < m_Nbg; ++i) in >> m_bg[i];
  for (int j = 0; j < m_Ndx; ++j) in >> m_dx[j];

  // Wypełnienie macierzy dE/dx
  for (int i = 0; i < m_Nbg; ++i) {
    for (int j = 0; j < m_Ndx; ++j) {
      in >> m_table[i * m_Ndx + j];
    }
  }
  in.close();
}

// Destruktor czyszczący pamięć
StBichsel::~StBichsel() {
  delete[] m_bg;
  delete[] m_dx;
  delete[] m_table;
}

// Publiczny interfejs - uśrednione dE/dx
double StBichsel::GetMeandEdx(double log10bg, double log10dx) const {
  // Bezpieczniki zakresu (Bichsel ma limity tabeli, np. log10bg od -1 do 5)
  if (log10bg < m_bg[0]) log10bg = m_bg[0];
  if (log10bg > m_bg[m_Nbg-1]) log10bg = m_bg[m_Nbg-1];
  if (log10dx < m_dx[0]) log10dx = m_dx[0];
  if (log10dx > m_dx[m_Ndx-1]) log10dx = m_dx[m_Ndx-1];

  // Wywołanie algorytmu interpolacyjnego
  return Interpolate(log10bg, log10dx);
}

// Prywatna metoda interpolacji dwuliniowej / splajnowej
double StBichsel::Interpolate(double log10bg, double log10dx) const {
  // Lokalne wyszukiwanie binarne indeksów (TMath::BinarySearch)
  int i = TMath::BinarySearch(m_Nbg, m_bg, log10bg);
  int j = TMath::BinarySearch(m_Ndx, m_dx, log10dx);

  if (i < 0) i = 0;
  if (i >= m_Nbg - 1) i = m_Nbg - 2;
  if (j < 0) j = 0;
  if (j >= m_Ndx - 1) j = m_Ndx - 2;

  // Wagi dla interpolacji liniowej pomiędzy sąsiednimi komórkami tabeli
  double t = (log10bg - m_bg[i]) / (m_bg[i+1] - m_bg[i]);
  double u = (log10dx - m_dx[j]) / (m_dx[j+1] - m_dx[j]);

  // Pobranie 4 sąsiadujących punktów z siatki Bichsela
  double z1 = m_table[i * m_Ndx + j];
  double z2 = m_table[(i+1) * m_Ndx + j];
  double z3 = m_table[(i+1) * m_Ndx + (j+1)];
  double z4 = m_table[i * m_Ndx + (j+1)];

  // Klasyczny wzór na interpolację dwuliniową (Bilinear Interpolation)
  double value = (1-t)*(1-u)*z1 + t*(1-u)*z2 + t*u*z3 + (1-t)*u*z4;

  return value; // Zwraca log10(dE/dx keV/cm)
}

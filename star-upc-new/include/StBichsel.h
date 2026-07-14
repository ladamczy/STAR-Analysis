#ifndef StBichsel_h
#define StBichsel_h

#include "TNamed.h"

class TPolyMarker;

class StBichsel : public TNamed {
 private:
  static StBichsel* fgInstance; // Wskaźnik na jedyną instancję klasy (Singleton)
  
  // Tabele struktur danych Bichsela (wczytywane z plików .dat / .root)
  int    m_Nbg;       // Liczba punktów osi beta*gamma
  double *m_bg;       // Tablica wartości log10(beta*gamma)
  int    m_Ndx;       // Liczba punktów osi grubości (dx)
  double *m_dx;       // Tablica wartości log10(charge^2 * dx)
  double *m_table;    // Dwuwymiarowa macierz przechowująca wartości dE/dx

  // Prywatne konstruktory – blokada tworzenia obiektów operatorem 'new'
  StBichsel(const char *tabName = "Bichsel");
  virtual ~StBichsel();

  // Metoda wewnętrzna do dwuwymiarowej interpolacji danych z tabeli
  double Interpolate(double log10bg, double log10dx) const;

 public:
  // Metoda dostępu do Singletonu
  static StBichsel* Instance();

  // Główna funkcja wywoływana przez użytkownika (zwraca log10(dE/dx_reprezentatywne))
  double GetMeandEdx(double log10bg, double log10dx) const;
  
  // Alternatywna funkcja zwracająca najbardziej prawdopodobną stratę energii (Most Probable Value)
  double GetMpwdEdx(double log10bg, double log10dx) const;

  ClassDef(StBichsel, 1) // Makro integracji z systemem ROOT
};

#endif

# STAR-Analysis
Repository for AGH-STAR analysis code

Repozytorium początkowo zawiera trzy projekty

Exclusive analysis (Patrycja)

Inclusive_analysis (Adam)

Data preparation (Przemysław)

# Użyteczne linki:

Diploma work of Patrycja Malinowska [PracaDyplomowa_PatrycjaMalinowska.pdf](https://github.com/user-attachments/files/18375738/PracaDyplomowa_PatrycjaMalinowska.pdf)

Praca Pana Rafała Sikory (https://cds.cern.ch/record/2747846/files/CERN-THESIS-2020-235.pdf)

Praca Pana Łukasza Fulka (https://cds.cern.ch/record/2724061/files/CERN-THESIS-2020-066.pdf)

Praca Pana Adama Wątroby (https://apd.usos.agh.edu.pl/diplomas/attachments/file/download/15182/)

Praca Pani Sverakovej (https://dspace.cvut.cz/bitstream/handle/10467/104062/F4-BP-2022-Sverakova-Michaela-bp_ejcf_22_sverakova.pdf)

Strona z poczatkiem i końcem napełnień RHIC w 2017 (https://wiki.bnl.gov/rhicspin/Run_17_polarization)

Kod do produkcji UPCDst z definicjami wszystkich klas w analizowanych danych.

https://github.com/adamjaro/star-upcDst

Potrzebne klasy z tego kodu są również w kartotekach src i include naszego kodu

Praca Pana Truhlara (https://drupal.star.bnl.gov/STAR/system/files/dp_ejcf_20_truhlar.pdf)

Nota techniczna na temat analizy starych danych STAR (https://drupal.star.bnl.gov/STAR/system/files/DiffractiveAnalyses_AnalysisNote_ver2p0.pdf)

O rekonstrukcji V0 [FengZhao_thesis.pdf](https://github.com/ladamczy/STAR-Analysis/files/13540392/FengZhao_thesis.pdf)

Praca magisterska Tomasa Truhlara https://drupal.star.bnl.gov/STAR/system/files/dp_ejcf_20_truhlar.pdf

Praca inzynierska Pana Rysia (https://github.com/ladamczy/STAR-Analysis/blob/main/share/Extention-of-the-ROOT-data-structure-describing-proton-proton-collisions-at-STAR-experiment-with-secondary-vertices-4.pdf )


# BeamPosition

Plik BeamPosition.cxx zawiera dwie funkcje ReadFillPositionData i FindPosition.

Funkcja vector <vector <double>> ReadFillPositionData(string file) służy do wczytywania danych z pliku csv https://github.com/ladamczy/STAR-Analysis/blob/main/share/Run7PolarizationWithPosition.csv (runNumberWithPosition.csv) w którym znajduje się informacja o numerze fillu (runu), pozycji wiązki oraz nachyleniu w kierunku x i y. Istnieją przypadki gdzie do jednego numeru fillu (runu) przyporządkowano dwie pozycje wiązki.
Funkcja zwraca dwuwymiarow wektor vector <vector <double>> Data z powyższymi informacjiami. W przypadku gdy na jeden numer fillu (runu) przypada jedna pozycja wiązki, wartościom drugiej pozycji i nachyenia przypisano nan.

Data[0] - fill (run) number
Data[1] - beam position x 
Data[2] - beam position y
Data[3] - beam position x2
Data[4] - beam position y2

Data[5] - beam slope x 
Data[6] - beam slope y
Data[7] - beam slope x2
Data[8] - beam slope y2

Funckja vector <double> FindPosition(int nFillNumber, double zPos, vector <double> &vFillNumber,  vector <double> &vXPosition, vector <double> &vYPosition,   vector <double> &vX2Position, vector <double> &vY2Position, vector <double> &vXSlope, vector <double> &vYSlope, vector <double> &vX2Slope, vector <double> &vY2Slope)

Funkcja FindPosition korzysta z danych wczytaych z wykorzystaniem ReadFillPositionData.
Do funkcji jest przekazywany numer fillu (runu) oraz na współrzędna z położenia werteksu zPos. W przypadku gdy na jeden numer fillu (runu) przypadają dwie pozycje wiązki obliczono ich średnią arytmetyczną. Funckja zwraca wektor położenia wiązki: vector <double> beam.  beam[0] = x, beam[1] = y


# Funkcje związane z detektorem Time-of-Flight

Klasa StUPCTrack posiada dwie nowe funkcje:
- `getT0(Double_t mass)`, która na podstawie podanej masy cząstki oblicza moment jej powstania (w ns, według tego samego zegara co ToF liczy getTofTime())
- `getMass(Double_t T0)`, która na podstawie czasu powstania cząstki oblicza jej masę (w GeV/c^2)

Klasa StUPCEvent posiada nową funkcję:
- `getMassSquared(StUPCTrack* track1, StUPCTrack* track2)`/`getMassSquared(Int_t iTrack1, Int_t iTrack2)`, która (przy założeniu, że dwie cząstki których pointery/indeksy w evencie mają taką samą masę i pochodzą z rozpadu tej samej cząstki) liczy kwadrat ich masy (w (GeV/c^2)^2)

Wszystkie powyższe funkcje w przypadku niepoprawnego czasu/długości drogi (niedodatnich) bądź braku flagi kToF przy śladzie, zwracają -999.

Dla cząstek pochodzących z tego samego zdarzenia (wierzchołek pierwotny/rozpad, przykład 1 i 2 na rysunku poniżej) o poprawnie dobranych masach różnica czasów ich powstania powinna być bliska 0
(np. dla pary proton-pion z rozpadu $\Lambda^0$ wyrażenie `protonTrack.getT0(0.938)-pionTrack.getT0(0.1395)` powinno być bliskie 0.).

Para cząstek pochodząca z różnych źródeł (przykład 3 i 4 na rysunku poniżej) bądź o źle dobranych masach tworzy wtedy gładkie continuum.

<img width="5296" height="1572" alt="TOFtiming" src="https://github.com/user-attachments/assets/7a197dc0-63dc-4947-a684-d3a951ef74f5" />

Przykład dla pary proton-pion z rozpadu $\Lambda^0$ (widać pik w okolicy 0 dla poprawnych par i tło dla reszty, z dopasowaniem piku gaussowskiego z tłem wielomianowym) poniżej. Ze względu na różne możliwości detekcji różnych typów cząstek, szerokości piku wokół 0 będą delikatnie różne.

<img width="1103" height="964" alt="image" src="https://github.com/user-attachments/assets/1baff835-a9a9-4758-a56d-891d59fc3cab" />



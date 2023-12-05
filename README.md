# STAR-Analysis
Repository for AGH-STAR analysis code

Repozytorium początkowo zawiera trzy projekty

Exclusive analysis (Patrycja)

Inclusive_analysis (Adam)

Data preparation (Przemysław)

# Użyteczne linki:

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


# BeamPosition

Plik BeamPosition.cxx zawiera dwie funkcje ReadFillPositionData i FindPosition.

Funkcja vector <vector <double>> ReadFillPositionData(string file) służy do wczytywania danych z pliku csv https://github.com/ladamczy/STAR-Analysis/blob/main/share/Run7PolarizationWithPosition.csv w którym znajduje się informacja o numerze fillu, pozycji wiązki oraz nachyleniu w kierunku x i y. Istnieją przypadki gdzie do jednego numeru fillu przyporządkowano dwie pozycje wiązki.
Funkcja zwraca dwuwymiarow wektor vector <vector <double>> Data z powyższymi informacjiami. W przypadku gdy na jeden numer fillu przypada jedna pozycja wiązki, wartościom drugiej pozycji i nachyenia przypisano nan.

Data[0] - fill number
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
Do funkcji jest przekazywany numer fillu oraz na współrzędna z położenia werteksu zPos. W przypadku gdy na jeden numer fillu przypadają dwie pozycje wiązki obliczono ich średnią arytmetyczną. Funckja zwraca wektor położenia wiązki: vector <double> beam.  beam[0] = x, beam[1] = y







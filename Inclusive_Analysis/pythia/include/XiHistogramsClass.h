//normal
#include <string>

//ROOT
#include "TH2D.h"
#include "TLorentzVector.h"

using namespace std;

class XiHistogramsClass{
private:
    TH2D* xi_e_w;
    TH1D* xi_multiplied;
    TH1D* xi_log;
public:
    XiHistogramsClass(string);
    ~XiHistogramsClass();
    void AddProton(TLorentzVector*, TLorentzVector*, double);
    void RetrieveHistograms(TH2D*, TH1D*, TH1D*);
};

XiHistogramsClass::XiHistogramsClass(string namepart){
    xi_e_w = new TH2D("xi_e_w", string("#Xi_{E} vs #Xi_{W}" + namepart + ";Xi_{E};Xi_{W}").c_str(), 100, -50, 50, 100, 50, 50);
    xi_multiplied = new TH1D("xi_multiplied", string("#Xi_{E}*#Xi_{W}" + namepart + ";Xi_{E}*Xi_{W};entries").c_str(), 100, -50, 50);
    xi_log = new TH1D("xi_log", string("ln(#frac{#Xi_{E}}{#Xi_{W}})" + namepart + ";ln(#frac{#Xi_{E}}{#Xi_{W}});entries").c_str(), 100, -50, 50);
}

XiHistogramsClass::~XiHistogramsClass(){
    delete xi_e_w;
    delete xi_multiplied;
    delete xi_log;
}

void XiHistogramsClass::AddProton(TLorentzVector* east_proton, TLorentzVector* west_proton, double s=510){
    double xi_e = (s/2-east_proton->E())/s*2;
    double xi_w = (s/2-west_proton->E())/s*2;
    xi_e_w->Fill(xi_e, xi_w);
    xi_multiplied->Fill(xi_e*xi_w);
    xi_log->Fill(TMath::Log(xi_e/xi_w));
}

void XiHistogramsClass::RetrieveHistograms(TH2D* xi_e_w_out, TH1D* xi_multiplied_out, TH1D* xi_log_out){
    xi_e_w_out = xi_e_w;
    xi_multiplied_out = xi_multiplied;
    xi_log_out = xi_log;
}

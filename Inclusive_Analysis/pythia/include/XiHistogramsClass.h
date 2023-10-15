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
    XiHistogramsClass(string, string, TH2D*, TH1D*, TH1D*);
    ~XiHistogramsClass();
    void AddProton(TLorentzVector*, TLorentzVector*, double);
    void RetrieveHistograms(TH2D*&, TH1D*&, TH1D*&);
};

XiHistogramsClass::XiHistogramsClass(string namepart, string identifier, TH2D* xi_e_w_in, TH1D* xi_multiplied_in, TH1D* xi_log_in){
    xi_e_w = xi_e_w_in;
    xi_multiplied = xi_multiplied_in;
    xi_log = xi_log_in;
    new(xi_e_w) TH2D(string("xi_e_w_"+identifier).c_str(), string("#Xi_{E} vs #Xi_{W}" + namepart + ";Xi_{E};Xi_{W}").c_str(), 50, 0, 1, 50, 0, 1);
    new(xi_multiplied) TH1D(string("xi_multiplied_"+identifier).c_str(), string("#Xi_{E}*#Xi_{W}" + namepart + ";Xi_{E}*Xi_{W};entries").c_str(), 60, 0, 0.3);
    new(xi_log) TH1D(string("xi_log_"+identifier).c_str(), string("ln(#frac{#Xi_{E}}{#Xi_{W}})" + namepart + ";ln(#frac{#Xi_{E}}{#Xi_{W}});entries").c_str(), 100, -10, 10);
}

XiHistogramsClass::~XiHistogramsClass(){
}

void XiHistogramsClass::AddProton(TLorentzVector* west_proton, TLorentzVector* east_proton, double s=510){
    double xi_e = (s/2-east_proton->E())/s*2;
    double xi_w = (s/2-west_proton->E())/s*2;
    xi_e_w->Fill(xi_e, xi_w);
    xi_multiplied->Fill(xi_e*xi_w);
    xi_log->Fill(TMath::Log(xi_e/xi_w));
}

void XiHistogramsClass::RetrieveHistograms(TH2D*& xi_e_w_out, TH1D*& xi_multiplied_out, TH1D*& xi_log_out){
    xi_e_w_out = xi_e_w;
    xi_multiplied_out = xi_multiplied;
    xi_log_out = xi_log;
}

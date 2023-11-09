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
    static int objectCount;
public:
    XiHistogramsClass(string, string);
    ~XiHistogramsClass();
    void AddProton(TLorentzVector*, TLorentzVector*, double);
    void RetrieveHistograms(TH2D&, TH1D&, TH1D&);
};

int XiHistogramsClass::objectCount;

XiHistogramsClass::XiHistogramsClass(string namepart, string identifier=""){
    if(identifier.length()==0){
        identifier = to_string(objectCount);
    }
    objectCount++;
    xi_e_w = new TH2D(string("xi_e_w_"+identifier).c_str(), string("#xi_{E} vs #xi_{W}"+namepart+";#xi_{E};#xi_{W}").c_str(), 50, 0, 1, 50, 0, 1);
    xi_multiplied = new TH1D(string("xi_multiplied_"+identifier).c_str(), string("#xi_{E}*#xi_{W}"+namepart+";#xi_{E}*#xi_{W};entries").c_str(), 60, 0, 0.3);
    xi_log = new TH1D(string("xi_log_"+identifier).c_str(), string("ln(#xi_{E}/#xi_{W})"+namepart+";ln(#xi_{E}/#xi_{W});entries").c_str(), 100, -10, 10);
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

void XiHistogramsClass::RetrieveHistograms(TH2D& xi_e_w_out, TH1D& xi_multiplied_out, TH1D& xi_log_out){
    xi_e_w_out = *xi_e_w;
    xi_multiplied_out = *xi_multiplied;
    xi_log_out = *xi_log;
}

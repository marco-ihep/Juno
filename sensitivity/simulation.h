#include <string>
#include <map>

#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooArgList.h"
#include "RooFFTConvPdf.h"
#include "RooNumIntConfig.h"
#include "RooGenericPdf.h"
#include "RooAbsReal.h"
#include "RooMsgService.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooGlobalFunc.h"
#include "RooPolyVar.h"
#include "RooWorkspace.h"
#include "RooChi2Var.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooMinuit.h"
#include "RooFitResult.h"


using namespace std;
using namespace RooFit;

enum Hierarchy {NH, IH};


class objectsToPrint {
    public:
        objectsToPrint() {};
        objectsToPrint(string arg1, map<string, string> arg2, Hierarchy arg3) : 
            multiple_spectra_name_(arg1), map_baseline_spectrumName_(arg2), thisHierarchy_(arg3) {};

        objectsToPrint(map<string, string> arg2, string arg1, Hierarchy arg3) : 
            multiple_spectra_name_(arg1), map_baseline_spectrumName_(arg2), thisHierarchy_(arg3) {};

        objectsToPrint & operator= (const objectsToPrint & other) { 
            multiple_spectra_name_ = other.multiple_spectra_name_;
            map_baseline_spectrumName_ = other.map_baseline_spectrumName_;  
            thisHierarchy_ = other.thisHierarchy_;
            return *this; 
        };

        string multiple_spectra_name_;
        map<string, string> map_baseline_spectrumName_;
        Hierarchy thisHierarchy_; 
};


class Simulation {
    private:
        RooWorkspace *w_nh;
        RooWorkspace *w_ih;
        TFile *output;
        int spectrum_index;
        //string smeared_spectrum_basename;
        string smeared_spectrum_basename_;
        string multiple_spectra_basename_;
        string path_;
        map<string, double> map_baseline_reactorPower_;

    public:
        Simulation(char *outputname);
        ~Simulation() {
            if(output && output->IsOpen())
                output->Close();
        };
        void fill_common_ws(RooWorkspace *ws);
        void init_common_vars(RooWorkspace *ws);
        RooFormulaVar *getPenalizedChi2(RooWorkspace *ws, RooChi2Var *chi2_data);
        void generate_nh_ws();
        objectsToPrint generate_multiple_spectra( Hierarchy thisHierarchy );
        void generate_ih_ws();
        string addSpectrum(RooWorkspace *ws, Hierarchy hierarchy, const char *baseline );
        //void fit(const char *nh_filename, const char *ih_filename);
        void fit(objectsToPrint & names_nh, objectsToPrint & names_ih, Hierarchy dataHierarchy, Hierarchy fcnHierarchy);
        void fit_nh(const char *nh_filename, const char *ih_filename);
        void fit_ih(const char *nh_filename, const char *ih_filename);
        void drawSpectra(objectsToPrint & names);
        void drawBaseline(objectsToPrint & names);
        void drawJointHierachies(objectsToPrint & names_nh, objectsToPrint & names_ih);

};


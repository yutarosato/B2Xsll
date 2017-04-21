#ifndef MKLIKELIHOOD
#define MKLIKELIHOOD

#include <iostream>
#include <math.h>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Set/Nominal_cut_selection.h"
#include "../Set/Branch.h"
#include "../Set/makeCut.h"

#include "draws_.h"

#include <vector>
#include <stdlib.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TF1.h>
#include <TArrow.h>

static const Char_t* dir_sig_before = "pdf_sig_before_bcs/O-U/";
static const Char_t* dir_bkg_before = "pdf_bkg_before_bcs";
static const Char_t* dir_sig_after  = "pdf_sig_after_bcs/O-U/";
static const Char_t* dir_bkg_after  = "pdf_bkg_after_bcs";
static const Char_t* event_type     = "all";
static const Char_t* stream         = "0";

Int_t main( Int_t argc, Char_t** argv );
Double_t Eval( Int_t fl_pdf, Double_t fvar, TH1D* hist, TF1* func=NULL,
	       Bool_t fl_message=!true, Bool_t fl_range=!true, Double_t fxmin=0, Double_t fxmax=0 );
Double_t Eval_hist( Double_t fvar, TH1D* hist, Bool_t fl_message=!true, Bool_t fl_range=!true, Double_t fxmin=0, Double_t fxmax=0 );
Double_t Eval_func( Double_t fvar, TF1*  func, Bool_t fl_message=!true, Bool_t fl_range=!true, Double_t fxmin=0, Double_t fxmax=0 );

#endif

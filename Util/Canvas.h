#ifndef Canvas_H
#define Canvas_H

#include <TROOT.h>
#include <TString.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TMath.h>

#include <sstream>

TCanvas* Canvas( const Char_t* name, const Char_t* title, Int_t n1=1, Int_t n2=0, Int_t width=600, Int_t height=600 );
void Deco( TH1D* hist, Int_t mode=0, Int_t color_l=1, Int_t color_m=1 );
// mode 0 : line
// mode 1 : marker w/  error bar
// mode 2 : marker w/o error bar
//void YRange( TH1D* hist, Double_t* range );
void Deco( TF1*          func,  Int_t mode=0, Int_t color_l=1, Int_t color_m=1 );
void Deco( TGraph*       graph, Int_t mode=0, Int_t color_l=1, Int_t color_m=1 );
void Deco( TGraphErrors* graph, Int_t mode=0, Int_t color_l=1, Int_t color_m=1 );

TH2D* Waku( Int_t n, TH1D**         hist,  const char* xlabel, Int_t ymin_flag=1,  Double_t space=0.15 ); // if ymin_flag=1, ymin is fixed   to zero.
TH2D* Waku( Int_t n, TGraph**       graph, const char* xlabel, Int_t ymin_flag=1,  Double_t space=0.15 ); // if ymin_flag=2, ymin is limited to positive value (for log-plot)
TH2D* Waku( Int_t n, TGraphErrors** graph, const char* xlabel, Int_t ymin_flag=1,  Double_t space=0.15 );

TH2D* Waku( Int_t n, TH1D**         hist,  const char* xlabel, const char* name, const char* title, Int_t ymin_flag=1,  Double_t space=0.15 ); // if ymin_flag=1, ymin is fixed to zero.
TH2D* Waku( Int_t n, TGraph**       graph, const char* xlabel, const char* name, const char* title, Int_t ymin_flag=1,  Double_t space=0.15 ); // if ymin_flag=1, ymin is fixed to zero.
TH2D* Waku( Int_t n, TGraphErrors** graph, const char* xlabel, const char* name, const char* title, Int_t ymin_flag=1,  Double_t space=0.15 ); // if ymin_flag=1, ymin is fixed to zero.

Double_t GetYMin( Int_t n, TH1D**         hist  );
Double_t GetYMin( Int_t n, TGraph**       graph );
Double_t GetYMin( Int_t n, TGraphErrors** graph );

Double_t GetYMax( Int_t n, TH1D**         hist  );
Double_t GetYMax( Int_t n, TGraph**       graph );
Double_t GetYMax( Int_t n, TGraphErrors** graph );

#endif

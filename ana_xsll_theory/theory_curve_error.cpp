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
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TArrow.h>
#include <TMath.h>
#include <Math/DistFunc.h>
#include <TRandom.h>
#include <TComplex.h>

const Bool_t flag_save = true; // outfile.eps and outfile.root
Int_t main( Int_t argc, Char_t** argv ){
  
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t nfile = 7;
  const Char_t* infile[nfile] = {
    "theory_error/syst_def/theory_curve_def.root",
    "theory_error/syst_mb/theory_curve_mb465.root",
    "theory_error/syst_mb/theory_curve_mb495.root",
    "theory_error/syst_ms/theory_curve_ms010.root",
    "theory_error/syst_ms/theory_curve_ms030.root",
    "theory_error/syst_mu/theory_curve_mu25.root",
    "theory_error/syst_mu/theory_curve_mu50.root",
  };
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TFile** file = new TFile*[nfile];
  TGraph** g = new TGraph*[nfile];
  Int_t npoint;
  for( Int_t i=0; i<nfile; i++ ){
    file[i] = new TFile( infile[i] );
    if( file[i]->IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file[i]->GetName() << std::endl, abort();
    g[i] = (TGraph*)file[i]->Get( "afb_true" );
    if( g[i] == NULL ) std::cerr << "[ABORT] can not find graph : " << std::endl, abort();
    if( i==0 || npoint > g[i]->GetN() ) npoint = g[i]->GetN();
  }

  npoint = g[0]->GetN(); // for uncomputable region
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TGraph* newg1 = new TGraph(); // error-sum
  TGraph* newg2 = new TGraph(); // min-max
  for( Int_t p=0; p<npoint; p++ ){
    Double_t x = g[0]->GetX()[p];
    Double_t y = g[0]->GetY()[p];
    Double_t ymax   = y;
    Double_t ymin   = y;
    Double_t yerr_L = 0.;
    Double_t yerr_H = 0.;
    Int_t    fl_cal = 0; // if can not calculate, fl_cal!=0 // for uncomputable region
    for( Int_t i=1; i<nfile; i++ ){
      if( p >= g[i]->GetN()    ){ fl_cal++; continue; } // for uncomputable region
      if( x != g[i]->GetX()[p] ) std::cout << "[WARNING] Wrong x-value" << std::endl;
      if( ymin > g[i]->GetY()[p] ) ymin = g[i]->GetY()[p];
      if( ymax < g[i]->GetY()[p] ) ymax = g[i]->GetY()[p];
      newg2->SetPoint( p,            x, ymax );
      newg2->SetPoint( 2*npoint-p-1, x, ymin );
      if     ( i==6 ) continue;
      else if( i==5 ){ // mu-scale
	if( fabs(g[5]->GetY()[p] - g[0]->GetY()[p]) > fabs(g[6]->GetY()[p] - g[0]->GetY()[p]) ){
	  yerr_L += (g[5]->GetY()[p] - g[0]->GetY()[p])*(g[5]->GetY()[p] - g[0]->GetY()[p]);
	  yerr_H += (g[5]->GetY()[p] - g[0]->GetY()[p])*(g[5]->GetY()[p] - g[0]->GetY()[p]);
	}else{
	  yerr_L += (g[6]->GetY()[p] - g[0]->GetY()[p])*(g[6]->GetY()[p] - g[0]->GetY()[p]);
	  yerr_H += (g[6]->GetY()[p] - g[0]->GetY()[p])*(g[6]->GetY()[p] - g[0]->GetY()[p]);
	}
	
      }else if( g[0]->GetY()[p] > g[i]->GetY()[p]  ) yerr_L += ( g[0]->GetY()[p]-g[i]->GetY()[p] )*( g[0]->GetY()[p]-g[i]->GetY()[p] );
      else                                           yerr_H += ( g[0]->GetY()[p]-g[i]->GetY()[p] )*( g[0]->GetY()[p]-g[i]->GetY()[p] );
    }
    yerr_L = sqrt( yerr_L );
    yerr_H = sqrt( yerr_H );

    if( fl_cal ){ // for uncomputable region
      newg2->SetPoint( p,            x, ymax );
      newg2->SetPoint( 2*npoint-p-1, x, 0    );
      yerr_L = y;
    }
    
    newg1->SetPoint( p,            x, y+yerr_H );
    newg1->SetPoint( 2*npoint-p-1, x, y-yerr_L );
  }
  newg1->SetName ( "theory_curve_error_quadsum" );
  newg2->SetName ( "theory_curve_error_minmax"  );
  newg1->SetTitle( "theory_curve_error_quadsum" );
  newg2->SetTitle( "theory_curve_error_minmax"  );
  newg1->SetFillColor(2);
  newg2->SetFillColor(3);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TCanvas* c = Canvas( "c1","c1" );
  g[0]->Draw("AL");
  newg1->Draw("F");
  newg2->Draw("F");
  g[0]->SetLineWidth(2);
  g[0]->Draw("Lsame");
  if( flag_save ){
    c->Print("pic/theory_curve_error.eps");
    TFile outfile( "pic/theory_curve_error.root", "RECREATE" );
    c    ->Write();
    g[0] ->Write(); // theory_curve
    newg1->Write(); // theory_curve with error (quad-sum)
    newg2->Write(); // theory_curve with error (min-max)
    outfile.Close();
  }

  std::cout << "finish!" << std::endl;
  app.Run();
  
  return 0;
}

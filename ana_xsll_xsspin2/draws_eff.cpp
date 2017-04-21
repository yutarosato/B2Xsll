#include <iostream>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Util/Manip.h"
#include "../Set/Nominal_cut_selection.h"
#include "../Set/Branch.h"
#include "../Set/makeCut.h"

#include "draws_.h"

#include <vector>
#include <stdlib.h>
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TArrow.h>

Int_t main( Int_t argc, Char_t** argv ){
  gROOT->SetBatch(true); // tmpppppp
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==2) ) std::cerr << "wrong input" << std::endl
			     << " Usage : ./draws_eff (int)fl_mode_ll" << std::endl
			     << "[fl_mode_ll] : 1(e), 0(mu)"           << std::endl, abort();
  Int_t fl_mode_ll = atoi( argv[1] ); // 1(e), 0(mu)
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::cout << "fl_mode_ll = " << fl_mode_ll << std::endl;
  const Int_t nc = 2;
  TChain* c_gen = new TChain("h12" );
  TChain* c_rec = new TChain("h511");
  c_gen->Add("/home/belle/syutaro/ewp/modules/xsll_sigmc_xsspin2/hbk/*.root");
  c_rec->Add("/home/belle/syutaro/ewp/modules/xsll_sigmc_xsspin2/hbk/*.root");
  std::cout << "[GEN] " << c_gen->GetEntries() << std::endl;
  std::cout << "[REC] " << c_rec->GetEntries() << std::endl;
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t add_cut_rec[nc][4096]; // spin0, spin1
  sTmp << "self==1";
  sTmp << "&& cc_m>0.2";
  sTmp << "&& gm_m_llg > 0.2 && gm_m_xs>1.10 && gm_m_xs<2.0";
  //sTmp << "&& gm_m_xs>1.1 && gm_m_xs<1.3"; // tmpppppppppppppppppppppppppp
  sTmp << Form( "&& gm_l==%d", fl_mode_ll );
  sTmp << Form( "&& rm_l==%d", fl_mode_ll );
  if     ( fl_mode_ll==1 ) sTmp << " && epp > 0.4 && emp > 0.4";
  else if( fl_mode_ll==0 ) sTmp << " && mpp > 0.8 && mmp > 0.8";
  sTmp << " && " << makeCut_5body_veto().c_str();
  sTmp << " && " << makeCut_unflavor_veto().c_str();  
  sTmp << " && " << cut_ccpi0;
  
  strcpy( add_cut_rec[0], Form("(%s)&&(abs(gm_xsid)==30343 || abs(gm_xsid)==30353)", (Char_t*)sTmp.str().c_str()) ); // Xs(spin0)
  strcpy( add_cut_rec[1], Form("(%s)&&(abs(gm_xsid)==30344 || abs(gm_xsid)==30354)", (Char_t*)sTmp.str().c_str()) ); // Xs(spin1)
  sTmp.str("");
  sTmp.clear();

  //+++++++++++
  Char_t add_cut_gen[nc][4096];
  sTmp << "llg_m > 0.2 && Xs_m>1.10 && Xs_m<2.0 && gm_fl_xs>0";
  //sTmp << "&& Xs_m>1.1 && Xs_m<1.3"; // tmpppppppppppppppppppppppppp
  sTmp << Form( "&& gm_l==%d", fl_mode_ll );

  strcpy( add_cut_gen[0], Form("(%s)&&(abs(Xs_id)==30343 || abs(Xs_id)==30353)", (Char_t*)sTmp.str().c_str()) ); // Xs(spin0)
  strcpy( add_cut_gen[1], Form("(%s)&&(abs(Xs_id)==30344 || abs(Xs_id)==30354)", (Char_t*)sTmp.str().c_str()) ); // Xs(spin1)
  sTmp.str("");
  sTmp.clear();


  std::cout << "[GEN(cut0)] " << c_gen->GetEntries(add_cut_gen[0]) << std::endl;
  std::cout << "[GEN(cut1)] " << c_gen->GetEntries(add_cut_gen[1]) << std::endl;
  std::cout << "[REC(cut0)] " << c_rec->GetEntries(add_cut_rec[0]) << std::endl;
  std::cout << "[REC(cut1)] " << c_rec->GetEntries(add_cut_rec[1]) << std::endl;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TH2D** hist_gen = new TH2D*[nc];
  TH2D** hist_rec = new TH2D*[nc];
  TH2D** hist_eff = new TH2D*[nc];


  /*
  const Int_t nslice  =   10;
  const Int_t    xbin =   20;
  const Double_t xmin =  0.0;
  const Double_t xmax = 20.0;
  */
  ///*
  const Int_t nslice  =     1;
  const Int_t    xbin =     5;
  const Double_t xmin =   1.0;
  const Double_t xmax =   6.0;
  //*/
  
  const Int_t    ybin =    5;
  const Double_t ymin =  0.0;
  const Double_t ymax =  1.0;
  /*
  const Int_t    ybin =    8;
  const Double_t ymin = -1.0;
  const Double_t ymax =  1.0;
  */
  for( Int_t ic=0; ic<nc; ic++ ){
    hist_gen[ic] = new TH2D( Form("gen_%d",ic), Form("gen_%d",ic), xbin, xmin, xmax, ybin, ymin, ymax );
    hist_rec[ic] = new TH2D( Form("rec_%d",ic), Form("rec_%d",ic), xbin, xmin, xmax, ybin, ymin, ymax );
    hist_eff[ic] = new TH2D( Form("eff_%d",ic), Form("eff_%d",ic), xbin, xmin, xmax, ybin, ymin, ymax );
    hist_gen[ic]->SetXTitle("q^{2} [GeV^{2}]");
    hist_gen[ic]->SetYTitle("cos#theta");
    hist_rec[ic]->SetXTitle("q^{2} [GeV^{2}]");
    hist_rec[ic]->SetYTitle("cos#theta");
    hist_eff[ic]->SetXTitle("q^{2} [GeV^{2}]");
    hist_eff[ic]->SetYTitle("cos#theta");

    c_gen->Project( Form("gen_%d",ic), "abs(coslp):llg_m*llg_m", add_cut_gen[ic] );
    c_rec->Project( Form("rec_%d",ic), "abs(coslp):cc_m*cc_m",   add_cut_rec[ic] );
    /*
    c_gen->Project( Form("gen_%d",ic), "coslp:llg_m*llg_m", add_cut_gen[ic] );
    c_rec->Project( Form("rec_%d",ic), "coslp:cc_m*cc_m",   add_cut_rec[ic] );
    */
    //c_rec->Project( Form("rec_%d",ic), "coslp:gm_m_llg*gm_m_llg",   add_cut_rec[ic] );
    for( Int_t ix=0; ix<xbin; ix++ ){
      for( Int_t iy=0; iy<ybin; iy++ ){
	Double_t ngen = hist_gen[ic]->GetBinContent( ix+1, iy+1 );
	Double_t nrec = hist_rec[ic]->GetBinContent( ix+1, iy+1 );
	Double_t eff  = ( ngen==0 ? 0.0 : nrec/ngen );
	hist_eff[ic]->SetBinContent( ix+1, iy+1, eff );
      }
    }
    hist_eff[ic]->SetMaximum(0.40);
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  
  TCanvas* c1 = Canvas( "c1","c1",3, 2 );
  {
    Int_t tmp_cnt = 0;
    for( Int_t ic=0; ic<nc; ic++ ){
      c1->cd(1+tmp_cnt++); hist_rec[ic]->Draw("COLZ");
      c1->cd(1+tmp_cnt++); hist_gen[ic]->Draw("COLZ");
      c1->cd(1+tmp_cnt++); hist_eff[ic]->Draw("COLZ");
    }
  }
  
  c1->Update();
  c1->Print( Form("pic/eff_xsspin2_lep%d.eps",fl_mode_ll) );

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  TH1D*** hist_gen_proj = new TH1D**[nc];
  TH1D*** hist_rec_proj = new TH1D**[nc];
  TH1D*** hist_eff_proj = new TH1D**[nc];

  for( Int_t ic=0; ic<nc; ic++ ){
    hist_gen_proj[ic] = new TH1D*[nslice];
    hist_rec_proj[ic] = new TH1D*[nslice];
    hist_eff_proj[ic] = new TH1D*[nslice];
    for( Int_t is=0; is<nslice; is++ ){
      //std::cout << is*(xbin/nslice)+1 << " : " << (is+1)*(xbin/nslice) << std::endl;
      hist_gen_proj[ic][is] = new TH1D( *hist_gen[ic]->ProjectionY( Form("_px%d",is), is*(xbin/nslice)+1, (is+1)*(xbin/nslice), "e") );
      hist_rec_proj[ic][is] = new TH1D( *hist_rec[ic]->ProjectionY( Form("_px%d",is), is*(xbin/nslice)+1, (is+1)*(xbin/nslice), "e") );
      hist_eff_proj[ic][is] = new TH1D( Form("eff_proj_%d_%d",ic,is), Form("eff_proj_%d_%d",ic,is), ybin, ymin, ymax );
      hist_rec_proj[ic][is]->SetLineColor  (ic+1);
      hist_gen_proj[ic][is]->SetLineColor  (ic+1);
      hist_eff_proj[ic][is]->SetLineColor  (ic+1);
      hist_rec_proj[ic][is]->SetMarkerColor(ic+1);
      hist_gen_proj[ic][is]->SetMarkerColor(ic+1);
      hist_eff_proj[ic][is]->SetMarkerColor(ic+1);
      if( ic==0 ) hist_eff_proj[ic][is]->SetLineWidth(2);
      hist_eff_proj[ic][is]->SetXTitle("cos#theta");
      
      for( Int_t iy=0; iy<ybin; iy++ ){
	Double_t ngen = hist_gen_proj[ic][is]->GetBinContent( iy+1 );
	Double_t nrec = hist_rec_proj[ic][is]->GetBinContent( iy+1 );
	Double_t eff  = ( ngen==0            ? 0.0 : nrec/ngen              );
	Double_t effE = ( ngen==0 || nrec==0 ? 0.0 : eff*sqrt((1-eff)/nrec) );
	hist_eff_proj[ic][is]->SetBinContent( iy+1, eff  );
	hist_eff_proj[ic][is]->SetBinError  ( iy+1, effE );
      }
      hist_gen_proj[ic][is]->SetMinimum(0.0);
      hist_rec_proj[ic][is]->SetMinimum(0.0);
      hist_eff_proj[ic][is]->SetMinimum(0.0);
    }
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  
  TCanvas* c2 = Canvas( "c2","c2", nslice );
  for( Int_t is=0; is<nslice; is++ ){
    c2->cd(is+1); 
    for( Int_t ic=0; ic<nc; ic++ ) hist_eff_proj[ic][is]->SetMaximum(0.25), hist_eff_proj[ic][is]->Draw( ic==0 ? "PE" : "PEsame" );
    //for( Int_t ic=0; ic<nc; ic++ ) hist_gen_proj[ic][is]->Draw( ic==0 ? "PE" : "PEsame" );
  }

  TLegend* leg = new TLegend( 0.75,0.75,0.99,0.99 );
  leg->AddEntry( hist_eff_proj[0][0],"Xs spin-0", "PL" );
  leg->AddEntry( hist_eff_proj[1][0],"Xs spin-1", "PL" );
  leg->Draw();
  
  c2->Update();
  c2->Print( Form("pic/eff_xsspin2_lep%d_proj.eps",fl_mode_ll) );

  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for( Int_t is=0; is<nslice; is++ ){
    for( Int_t ic=0; ic<nc; ic++ ){
      Double_t nf  = 0.0;
      Double_t nb  = 0.0;
      Double_t nfE = 0.0;
      Double_t nbE = 0.0;
      for( Int_t ibin=0; ibin<hist_eff_proj[ic][is]->GetNbinsX(); ibin++ ){
	Double_t w  = hist_eff_proj[ic][is]->GetBinWidth  ( ibin+1 );
	Double_t xl = hist_eff_proj[ic][is]->GetBinLowEdge( ibin+1 );
	Double_t xh = xl + w;
	Double_t s_f = (( xl + 1.0) + ( xh + 1.0))*w/2.0;
	Double_t s_b = ((-xl + 1.0) + (-xh + 1.0))*w/2.0;
	nf  +=     s_f * hist_eff_proj[ic][is]->GetBinContent(ibin+1);
	nb  +=     s_b * hist_eff_proj[ic][is]->GetBinContent(ibin+1);
	nfE += pow(s_f * hist_eff_proj[ic][is]->GetBinError  (ibin+1), 2);
	nbE += pow(s_b * hist_eff_proj[ic][is]->GetBinError  (ibin+1), 2);
      }
      nfE = sqrt(nfE);
      nbE = sqrt(nbE);
      Double_t x  = nb/nf;
      Double_t xE = x * sqrt( pow(nfE/nf,2) + pow(nbE/nb,2) );
      Double_t afb  = (nf - nb)/(nf + nb);
      Double_t afbE = 2/(1+x)*xE;
      std::cout << "AFB = "   << afb << " +- " << afbE << ", "
		<< "alpha = " << 0.50/afb << " +- " << 0.50/afb/afb*afbE << ", " 
		<< "NF = "  << nf  << " +- " << nfE << ", "
		<< "NB = "  << nb  << " +- " << nbE << ", "
		<< "x = "   << x   << " +- " << xE
		<< std::endl;
      
    }
    std::cout << std::endl;
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::cout << "finish" << std::flush;  
  app.Run();

  return 0;
}

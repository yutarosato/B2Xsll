#include <iostream>

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
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TArrow.h>


Int_t main( Int_t argc, Char_t** argv ){
  using namespace sigmc;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_2d (char*)setname (double)used_set (int)fl_mode_ll [(int)fl_appRun]" << std::endl
					<< "[  setname  ] A,B,AB,A-C,..." << std::endl, abort();
  
  const Char_t*  setname    = argv[1];
  const Double_t used_nset  = atof(argv[2]); 
  const Int_t    fl_mode_ll = atoi(argv[3]); // 1(e), 0(mu)
  Int_t          fl_appRun  = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t    Nchain          = 1; // fixed
  const Int_t    Nhist           = 1; // fixed
  const Int_t    nfile[Nchain]   = {0};
  const Double_t scale_event_sig = sigmc_amount*(used_nset/(Double_t)nset); // sigmc : N -> N/alpha
  const Bool_t   flag_save       = true;
  const Bool_t   flag_k4pi       = true; // 1(veto  K4pi    modes)
  const Bool_t   flag_unflavor   = true; // 1(veto unflavor modes)
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    sTmp << indir << "/sigMC_*_set[" << setname << "]";
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] = {
    { 1},
  };
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t fl_q2 = 0;
  Char_t* LRNB_cut = new Char_t[2048];
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.91, 0.36, 0.94, 0.92, 0.87, 0.57, 0.92, 0.89).c_str() ); // 2d NB_lep(bcs=bb)
  LRNB_cut = "1";
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;
  
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t j=0; j<Nchain; j++ ) add_cut[j] = new Char_t[4096];
  sTmp << "self==1";
  if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
  if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
  strcpy( add_cut[0], (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  using namespace q2_theta_nonuniform;
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH2D**    tmphist = new TH2D*  [Nchain];
  TH1D**    hist1   = new TH1D*  [Nhist]; // for A_{FB}
  TH2D**    hist2   = new TH2D*  [Nhist]; // for q2 v.s. cos-theta
  TH1D**    histx   = new TH1D*  [Nhist];
  TH1D**    histy   = new TH1D*  [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1", 4, 1 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis) << std::endl;
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.26 ); // Mbc sideband
  }
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j ) << std::endl;
    chain[j]->GetCut()->Display(0);
    chain[j]->MakeTree();
  }

  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  std::cout << " ************************ make tmphist *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    //tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[Nchain]->GetChange()), xbin,xmin,xmax, ybin,ymin,ymax );
    tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[Nchain]->GetChange()), xbin,xbins[fl_mode_ll], ybin,ymin,ymax );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] ); 
    tmphist[j]->Sumw2();
    tmphist[j]->Scale( 1/scale_event_sig );
    // +++++++ display ++++++++++++++++++++++++++++++++++
    std::cout << std::setw(10) << Form("<tmphist%d>",j);
    std::cout << std::setw(12) << tmphist[j]->GetEntries() << " events : " 
	      << add_cut[j] << std::endl;
  }

  // +++++++ make hist2  ++++++++++++++++++++++++++++++++++
  std::cout << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    //hist2[i] = new TH2D( Form("hist%d",i),Form("hist%d",i), xbin,xmin,xmax, ybin,ymin,ymax ); // uniform binning
    hist2[i] = new TH2D( Form("hist%d",i),Form("hist%d",i), xbin,xbins[fl_mode_ll], ybin,ymin,ymax ); // non-uniform binning
    ((TGaxis*)hist2[i]->GetXaxis())->SetMaxDigits(3);
    ((TGaxis*)hist2[i]->GetYaxis())->SetMaxDigits(3);
    hist2[i]->GetXaxis()->CenterTitle();
    hist2[i]->GetYaxis()->CenterTitle();
    hist2[i]->SetXTitle( xlabel );
    hist2[i]->SetYTitle( ylabel );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist2[i]->Add( tmphist[j], (Double_t) add[i][j] );
    }
    // +++++++ display ++++++++++++++++++++++++++++++++++
    std::cout << std::setw(10) << Form("<hist%d>",i);
    //std::cout << std::setw(12) << hist2[i]->GetEntries() << " events, Added-Files( ";
    std::cout << std::setw(12) << hist2[i]->Integral() << " events, Added-Files( ";
    for( Int_t j=0; j<Nchain; j++ ) if( add[i][j] ) std::cout << add[i][j]*j << ",";
    std::cout << " )" << std::endl; 
  }

  // +++++++ make hist[xy] ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nhist; i++ ){
    histx[i] = hist2[i]->ProjectionX();
    histy[i] = hist2[i]->ProjectionY();
    histx[i]->Scale( 1/histx[i]->Integral() );
    histy[i]->Scale( 1/histy[i]->Integral() );
    Deco( histx[i], 0, i+1, i+1 );
  }

  // +++++++ make hist1 ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nhist; i++ ){
    //hist1[i] = new TH1D( Form("afb%d",i),Form("afb%d",i), xbin,xmin,xmax ); // uniform binning
    hist1[i] = new TH1D( Form("afb%d",i),Form("afb%d",i), xbin,xbins[fl_mode_ll]     ); // non-uniform binning
    Deco( hist1[i], 0, i+1, i+1 );
    for( Int_t k=0; k<xbin; k++ ){
      Double_t theta_pE = 0;
      Double_t theta_p = hist2[i]->IntegralAndError( k+1,k+2, 6, 10, theta_pE );
      Double_t theta_mE = 0;
      Double_t theta_m = hist2[i]->IntegralAndError( k+1,k+2, 1,  5, theta_mE );
      Double_t afb = (theta_p-theta_m)/(theta_p+theta_m);
      Double_t afbE = 2/(theta_p+theta_m)/(theta_p+theta_m)*sqrt(theta_m*theta_m*theta_pE*theta_pE + theta_p*theta_p*theta_mE*theta_mE);
      std::cout << "bin = " << k+1 << ": "
		<< "Np = " << theta_p
		<< " +- "  << theta_p << ", "
		<< "Nm = " << theta_m
		<< " +- "  << theta_m << " -> "
		<< "AFB = " << afb
		<< " +- " << afbE
		<< std::endl;
      hist1[i]->SetBinContent( k+1, afb );
      hist1[i]->SetBinError  ( k+1, afbE );
    }
  }
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  Double_t allymax_x = GetYMax( Nhist, histx);
  Double_t allymax_y = GetYMax( Nhist, histy);
  for( Int_t i=0; i<Nhist; i++ ){
    histx[i]->SetMaximum( 1.1*allymax_x );
    histy[i]->SetMaximum( 1.1*allymax_y );
  }
  c1->Draw();

  c1->cd(1);
  hist2[0]->Draw( "COLZ" );
  c1->cd(2);
  hist1[0]->Draw();
  c1->cd(3);
  histx[0]->Draw();
  c1->cd(4);
  histy[0]->Draw();
  
  // +++++++ save +++++++++++++++++++++++++++++++++++++++
  c1->Update();

  if( flag_save ) c1->Print( Form("pic/2d_q2_theta_sig_lep%d_set%s.eps", fl_mode_ll, setname) );
  std::cout << "finish" << std::flush;  
  if( fl_appRun ) app.Run();
  
  delete[] chain;
  delete[] tmphist;
  delete[] hist1;
  delete[] hist2;
  delete[] histx;
  delete[] histy;
  delete   c1;
  
  return 0;
}

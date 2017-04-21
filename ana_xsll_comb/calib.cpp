#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Util/Stat.h"
#include "../Set/Nominal_cut_selection.h"
#include "../Set/Branch.h"

#include "draws_.h"

#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>
#include <TString.h>
#include <TFile.h>

Int_t main( Int_t argc, Char_t** argv ){
  if( argc!=3 ) std::cerr << "wrong input" << std::endl
			  << " Usage : ./calib (char*)infile (char*)outdir" << std::endl,abort();
  const Char_t* infile = argv[1];
  const Char_t* outdir = argv[2];
  
  std::string basename = gSystem->BaseName( infile );  

  if( basename.find(".root") == std::string::npos ){
    std::cerr << "[infile] " << std::setw(80) << infile
	      << " -> wrong file name" << std::endl;
    abort();
  }

  Int_t fl_sbr = -1; // 0(sigMC), 1(gMC), 2(RD)
  if     ( basename.find("sigMC") != std::string::npos ) fl_sbr = 0;
  else if( basename.find("gMC")   != std::string::npos ) fl_sbr = 1;
  else if( basename.find("CC")    != std::string::npos ) fl_sbr = 1;
  else if( basename.find("RD")    != std::string::npos ) fl_sbr = 2;
  else std::cerr << "[ABORT] wrong file name : " << infile << std::endl, abort();

  basename.erase( basename.rfind(".root") );
  // --------------------------------------------------------------------  
  const Char_t* tname_B       = "h511";
  const Char_t* tname_Gen     = "h12";
  const Double_t scale_factor_chi2_sig[2] = {1.14, 1.14}; // [mm,ee]
  const Double_t scale_factor_chi2_bkg    = 1.04;
  const Double_t scale_factor_dzll_sig[2] = {1.28, 1.28}; // [mm,ee]
  const Double_t scale_factor_dzll_bkg    = 1.17;
  const Double_t scale_factor_de_mean[2]  = {0.0009, 0.0019}; // [mm,ee]
  const Double_t scale_factor_de_sigma[2] = {1.15,   1.13};   // [mm,ee]
  // --------------------------------------------------------------------  
  TChain* chain_B   = new TChain( tname_B   );
  TChain* chain_Gen = new TChain( tname_Gen );
  chain_B  ->Add( infile );
  chain_Gen->Add( infile );
  // --------------------------------------------------------------------  
  Float_t rm_l, kfbdgf;
  Float_t kfbcl,      kfbchi,      dzll3d,      de;
  Float_t kfbclorg,   kfbchiorg,   dzll3dorg,   deorg;
  Float_t kfbclcalib, kfbchicalib, dzll3dcalib, decalib;
  
  chain_B->SetBranchAddress( "rm_l",   &rm_l   );
  chain_B->SetBranchAddress( "kfbdgf", &kfbdgf );
  chain_B->SetBranchAddress( "kfbcl",  &kfbcl  );
  chain_B->SetBranchAddress( "kfbchi", &kfbchi );
  chain_B->SetBranchAddress( "dzll3d", &dzll3d );
  chain_B->SetBranchAddress( "de",     &de     );

  TTree* newtree_B = chain_B->CloneTree(0);
  newtree_B->Branch("kfbclorg",    &kfbclorg,    "kfbclorg/F"    );
  newtree_B->Branch("kfbchiorg",   &kfbchiorg,   "kfbchiorg/F"   );
  newtree_B->Branch("dzll3dorg",   &dzll3dorg,   "dzll3dorg/F"   );
  newtree_B->Branch("deorg",       &deorg,       "deorg/F"       );
  newtree_B->Branch("kfbclcalib",  &kfbclcalib,  "kfbclcalib/F"  );
  newtree_B->Branch("kfbchicalib", &kfbchicalib, "kfbchicalib/F" );
  newtree_B->Branch("dzll3dcalib", &dzll3dcalib, "dzll3dcalib/F" );
  newtree_B->Branch("decalib",     &decalib,     "decalib/F"     );

  Int_t  nevt      = 0;
  while( chain_B->GetEntry(nevt, 0) ) {
    Int_t fl_mode_ll = rm_l ? 1 : 0;
    
    // original values
    kfbclorg  = kfbcl;
    kfbchiorg = kfbchi;
    dzll3dorg = dzll3d;
    deorg     = de;

    // calibrated values by signal scale factor
    if( fl_sbr!=2 ){ // for MC
      dzll3dcalib = scale_factor_dzll_sig[fl_mode_ll] * dzll3dorg;    
      kfbchicalib = scale_factor_chi2_sig[fl_mode_ll] * kfbchiorg;
      kfbclcalib  = chisq2cl( (Int_t)kfbdgf, kfbchicalib );
      decalib = scale_factor_de_sigma[fl_mode_ll] * deorg - scale_factor_de_mean[fl_mode_ll];
    }else{ // for real-data
      dzll3dcalib = dzll3dorg;
      kfbchicalib = kfbchiorg;
      kfbclcalib  = kfbclorg;
      decalib     = deorg;
    }

    // defalut values
    if( fl_sbr==0 || fl_sbr==1 ){
      if( fl_sbr==0 ){ // sig
	dzll3d = scale_factor_dzll_sig[fl_mode_ll] * dzll3dorg;    
	kfbchi = scale_factor_chi2_sig[fl_mode_ll] * kfbchiorg;
      }else if( fl_sbr==1 ){ // bkg
	dzll3d = scale_factor_dzll_bkg * dzll3dorg;    
	kfbchi = scale_factor_chi2_bkg * kfbchiorg;
      }
      kfbcl = chisq2cl( (Int_t)kfbdgf, kfbchi );
    }
    if( fl_sbr==0 ) de = decalib;

    newtree_B->Fill();
    nevt++;
  }

  std::cout << "[infile] " << std::setw(80) << infile
	    << "(" << tname_B    << ", " << std::setw(6) << chain_B  ->GetEntries() << ")"
	    << std::endl;
  TTree* newtree_Gen = chain_Gen->CloneTree();
  // --------------------------------------------------------------------
  const Char_t* outfile = Form( "%s%s_calib.root",outdir,basename.c_str() );
  TFile*        rootf   = new TFile( outfile, "RECREATE" );
  newtree_B->Write();
  if( strstr(infile, "_m1m_") != NULL || strstr(infile, "_m9999m_") != NULL ){
    if( fl_sbr==0 ) newtree_Gen->Write();
  }
  rootf->Close();

  delete chain_B;
  delete newtree_B;
  delete newtree_Gen;
  delete rootf;

  return 0;
}

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
#include <TF1.h>

Double_t Eval_func( Double_t fvar, TF1* func ){
  Double_t area = func->Integral( func->GetXmin(), func->GetXmax() );
  Double_t var  = func->Eval( fvar );
  return var/area;  
}


Int_t main( Int_t argc, Char_t** argv ){
  if( argc!=3 ) std::cerr << "wrong input" << std::endl
			  << " Usage : ./cal_delh (char*)infile (char*)outdir" << std::endl,abort();
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
  const Char_t* fl_type       = "tot";

  // --------------------------------------------------------------------  
  TFile file_de_lep1( "~/ewp/ana/ana_xsll_comb/NB_lep_calib/pdf/pdf_de_lep1_s03-5_setK-U.root" );
  TFile file_de_lep0( "~/ewp/ana/ana_xsll_comb/NB_lep_calib/pdf/pdf_de_lep0_s03-5_setK-U.root" );
  if( file_de_lep1.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_de_lep1.GetName() << std::endl, abort();
  if( file_de_lep0.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_de_lep0.GetName() << std::endl, abort();
  TF1* func_sig_de_lep1         = (TF1*)file_de_lep1.Get     (      "pdf_func_de_lep1_sig"             );
  TF1* func_sig_de_lep1_pi0     = (TF1*)file_de_lep1.Get     (      "pdf_func_de_lep1_pi0_sig"         );
  TF1* func_sig_de_lep0         = (TF1*)file_de_lep0.Get     (      "pdf_func_de_lep0_sig"             );
  TF1* func_sig_de_lep0_pi0     = (TF1*)file_de_lep0.Get     (      "pdf_func_de_lep0_pi0_sig"         );
  TF1* func_bkg_de_lep1         = (TF1*)file_de_lep1.Get     ( Form("pdf_func_de_lep1_%s_bkg",fl_type) );
  TF1* func_bkg_de_lep0         = (TF1*)file_de_lep0.Get     ( Form("pdf_func_de_lep0_%s_bkg",fl_type) );
  if( func_sig_de_lep1         == NULL ) std::cerr << "[ABORT] can not find func  1 : " << std::endl, abort();
  if( func_sig_de_lep1_pi0     == NULL ) std::cerr << "[ABORT] can not find func  2 : " << std::endl, abort();
  if( func_sig_de_lep0         == NULL ) std::cerr << "[ABORT] can not find func  3 : " << std::endl, abort();
  if( func_sig_de_lep0_pi0     == NULL ) std::cerr << "[ABORT] can not find func  4 : " << std::endl, abort();
  if( func_bkg_de_lep1         == NULL ) std::cerr << "[ABORT] can not find func  5 : " << std::endl, abort();
  if( func_bkg_de_lep0         == NULL ) std::cerr << "[ABORT] can not find func  8 : " << std::endl, abort();

  // --------------------------------------------------------------------  
  const Char_t* outfile = Form( "%s%s_delh.root",outdir,basename.c_str() );
  TFile*        rootf   = new TFile( outfile, "RECREATE" );

  // --------------------------------------------------------------------  
  TChain* chain_B = new TChain( tname_B   );
  chain_B ->Add( infile );
  
  // --------------------------------------------------------------------  
  Float_t rm_l, rm_xs, de, delh;
  
  chain_B->SetBranchAddress( "rm_l",   &rm_l   );
  chain_B->SetBranchAddress( "rm_xs",  &rm_xs  );
  chain_B->SetBranchAddress( "de",     &de     );

  TTree* newtree_B = chain_B->CloneTree(0);
  newtree_B->Branch("delh", &delh, "delh/F" );

  // --------------------------------------------------------------------  

  Int_t nevt = 0;
  while( chain_B->GetEntry(nevt, 0) ) {
    Int_t fl_mode_ll = rm_l ? 1 : 0;
    Double_t ls = 0;
    Double_t lb = 0;
    if( fl_mode_ll==1 ){ // ee
      if( de<-0.10 || de>0.05 ){
	ls = 0;
	lb = 0;
      }else{
	if( rm_xs > 999 ) ls = Eval_func( de, func_sig_de_lep1_pi0 ); // w/  pi0
	else              ls = Eval_func( de, func_sig_de_lep1     ); // w/o pi0
	lb = Eval_func( de, func_bkg_de_lep1 ); // bkg(tot)
      }
    }else if( fl_mode_ll==0 ){ // mu
      if( de<-0.05 || de>0.05 ){
	ls = 0;
	lb = 0;
      }else{
	if( rm_xs > 999 ) ls = Eval_func( de, func_sig_de_lep0_pi0 ); // w/  pi0
	else              ls = Eval_func( de, func_sig_de_lep0     ); // w/o pi0
	lb = Eval_func( de, func_bkg_de_lep1 ); // bkg(tot)
      }
    }
    
    if( ls==0 && lb==0 ) delh = -999;
    else                 delh = ls/(ls+lb);

    nevt++;
    if( (rm_l==1 && de>-0.10 && de<0.05) || (rm_l==0 && de>-0.05 && de<0.05) ){
      de = delh;
      newtree_B->Fill();
    }else continue;
  }
  
  std::cout << "[infile] " << std::setw(80) << infile
	    << "(" << tname_B    << ", " << std::setw(6) << chain_B  ->GetEntries() << ")"
	    << std::endl;

  // --------------------------------------------------------------------

  newtree_B->Write();
  rootf->Close();

  return 0;
}



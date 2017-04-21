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
#include "../Set_fake_lepton/Nominal_cut_selection.h"
#include "../Set_fake_lepton/Branch.h"
#include "../Set_fake_lepton/makeCut.h"

#include "draws_.h"

#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>
#include <TString.h>
#include <TFile.h>
#include <TPaveText.h>

Int_t main( Int_t argc, Char_t** argv ){
  if( argc!=4 ){
    std::cerr << "wrong input" << std::endl
	      << " Usage : ./merge_cut2 (char*)indir (char*)outdir (char*)exp"
	      << std::endl;
    abort();
  }

  const Char_t* indir  = argv[1];
  const Char_t* outdir = argv[2];
  const Char_t* exp    = argv[3];
  // --------------------------------------------------------------------
  const Bool_t  fl_message = !true;
  const Char_t* tname_B    = "h511";
  TFile* file_fl  = new TFile( "fakerate_dstr_data/fake_rate.root" );
  TH2D** hist_fl  = new TH2D*[2];
  TH2D** hist_flE = new TH2D*[2];
  hist_fl [0] =  (TH2D*)file_fl->Get( "hist2_fake" ); // mm
  hist_fl [1] =  (TH2D*)file_fl->Get( "hist3_fake" ); // ee
  hist_flE[0] =  (TH2D*)file_fl->Get( "hist2_err"  ); // mm
  hist_flE[1] =  (TH2D*)file_fl->Get( "hist3_err"  ); // ee
  // --------------------------------------------------------------------
  TChain* chain_B = new TChain( tname_B );
  
  std::stringstream sTmp;
  sTmp << indir << "/RD_e0" << exp << "*.root";

  Int_t nfile_B = chain_B->Add( sTmp.str().c_str() );
  if( !nfile_B ) std::cerr << "[ABORT] NO file : " << sTmp.str().c_str() << std::endl, abort();

  // --------------------------------------------------------------------  
  TCut cut_lrnb = "1";
  // --------------------------------------------------------------------
  const Char_t* outfile = Form( "%s/RD_e0%s_caseB_double_weight.root",outdir,exp );
  TFile* rootf  = new TFile( outfile, "RECREATE" );
  TTree* tree_B = new TTree();
  tree_B = chain_B->CopyTree( cut_lrnb );
  // ------------------------------------------------------------------------
  Float_t rm_l;                             // Mode
  Float_t epp, emp, mpp, mmp;               // Momentum
  Float_t lpc, lmc;                         // Direction
  Float_t lp_eid, lm_eid, lp_muid, lm_muid; // Lepton-PID
  Float_t fl_p[2], fl_c[2];                 // Momentum and Direction of fake-lepton [l+, l-]
  Float_t weight, weightE;                  // Weight ( New Branch )
  
  tree_B->SetBranchAddress( "rm_l",    &rm_l    );
  tree_B->SetBranchAddress( "epp",     &epp     );
  tree_B->SetBranchAddress( "emp",     &emp     );
  tree_B->SetBranchAddress( "mpp",     &mpp     );
  tree_B->SetBranchAddress( "mmp",     &mmp     );
  tree_B->SetBranchAddress( "lpc",     &lpc     );
  tree_B->SetBranchAddress( "lmc",     &lmc     );
  tree_B->SetBranchAddress( "lp_eid",  &lp_eid  );
  tree_B->SetBranchAddress( "lm_eid",  &lm_eid  );
  tree_B->SetBranchAddress( "lp_muid", &lp_muid );
  tree_B->SetBranchAddress( "lm_muid", &lm_muid );

  tree_B->Branch( "weight",  &weight,  "weight/F"  );
  tree_B->Branch( "weightE", &weightE, "weightE/F" );

  // ------------------------------------------------------------------------

  TTree* newtree_B = tree_B->CloneTree( 0 );

  Int_t nev = 0;
  while (tree_B->GetEntry(nev, 0)) {
    nev++;
    if( rm_l == 1 ){ // electron mdoe
	fl_p[0] = epp;
	fl_c[0] = lpc;
	fl_p[1] = emp;
	fl_c[1] = lmc;
    }else if( rm_l == 0 ){ // muon mode
	fl_p[0] = mpp;
	fl_c[0] = lpc;
	fl_p[1] = mmp;
	fl_c[1] = lmc;
    }else{
      std::cerr << "[ABORT] : Invalid modes(rm_l) : " << rm_l << std::endl, abort();
    }

    Int_t   bin_x[2];
    Int_t   bin_y[2];
    Float_t tmp_weight [2];
    Float_t tmp_weightE[2];
    for( Int_t i=0; i<2; i++ ){
      bin_x[i]       =  hist_fl[(Int_t)rm_l]->FindBin((Double_t)fl_c[i]);
      bin_y[i]       = (hist_fl[(Int_t)rm_l]->FindBin((Double_t)fl_c[i], (Double_t)fl_p[i]))/(hist_fl[(Int_t)rm_l]->GetNbinsX()+2);
      tmp_weight [i] = hist_fl [(Int_t)rm_l]->GetBinContent(bin_x[i], bin_y[i]);
      tmp_weightE[i] = hist_flE[(Int_t)rm_l]->GetBinContent(bin_x[i], bin_y[i]);
      if( bin_x[i]==8 && bin_y[i]==1 ){
	tmp_weight[i]  = (hist_fl [(Int_t)rm_l]->GetBinContent(7,1) + hist_fl [(Int_t)rm_l]->GetBinContent(9,1) )/2.0;
	tmp_weightE[i] = sqrt( hist_flE[(Int_t)rm_l]->GetBinContent(7,1)*hist_flE[(Int_t)rm_l]->GetBinContent(7,1)
			+hist_flE[(Int_t)rm_l]->GetBinContent(9,1)*hist_flE[(Int_t)rm_l]->GetBinContent(9,1) )/2.0;
      }
      if( bin_x[i]==3 && bin_y[i]==1 ){
	tmp_weight[i]  = (hist_fl [(Int_t)rm_l]->GetBinContent(2,1) + hist_fl [(Int_t)rm_l]->GetBinContent(4,1) )/2.0;
	tmp_weightE[i] = sqrt( hist_flE[(Int_t)rm_l]->GetBinContent(2,1)*hist_flE[(Int_t)rm_l]->GetBinContent(2,1)
			+hist_flE[(Int_t)rm_l]->GetBinContent(4,1)*hist_flE[(Int_t)rm_l]->GetBinContent(4,1) )/2.0;
      }
    }
    weight = tmp_weight[0] * tmp_weight[1];
    if( weight==0 ) weightE = 0;
    else            weightE = weight * sqrt( tmp_weightE[0]*tmp_weightE[0]/tmp_weight[0]/tmp_weight[0] + tmp_weightE[1]*tmp_weightE[1]/tmp_weight[1]/tmp_weight[1] );
    
    if( fl_message || (weight==0 && weightE==0) ) std::cout << "[CHECK] "
							    << std::setw( 7) << std::right << nev            << " : "
							    << std::setw( 2) << std::right << rm_l           << " : ("
							    << std::setw(12) << std::right << fl_c[0]        << ", "
							    << std::setw(12) << std::right << fl_p[0]        << ") -> ("
							    << std::setw( 3) << std::right << bin_x[0]       << ", "
							    << std::setw( 3) << std::right << bin_y[0]       << ") => "
							    << std::setw(12) << std::right << tmp_weight[0]  << " +- "
							    << std::setw(12) << std::right << tmp_weightE[0] << ", "
							    << std::setw(12) << std::right << fl_c[1]        << ", "
							    << std::setw(12) << std::right << fl_p[1]        << ") -> ("
							    << std::setw( 3) << std::right << bin_x[1]       << ", "
							    << std::setw( 3) << std::right << bin_y[1]       << ") => "
							    << std::setw(12) << std::right << tmp_weight[1]  << " +- "
							    << std::setw(12) << std::right << tmp_weightE[1] << std::endl;
    newtree_B->Fill();
  }

  std::cout << "exp0" << std::setw(3) << std::left  << exp
	    << ", "   << std::setw(3) << std::right << nfile_B << " files(" << tname_B    << ", " << std::setw(9) << std::right << chain_B->GetEntries()
	    << " -> "
	    << std::setw(9) << std::right << newtree_B  ->GetEntries()
	    << "[" << std::setw(2) << std::right << int(100*((double)newtree_B->GetEntries()/chain_B->GetEntries())) << "%]"
	    << ")"
	    << std::endl;
  // --------------------------------------------------------------------
  
  newtree_B->Write();
  //rootf    ->Close();
  
  delete chain_B;
  delete tree_B;
  delete newtree_B;
  delete rootf;
  delete file_fl;

  return 0;
}

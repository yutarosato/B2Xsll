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

#include "draws_.h"

#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>
#include <TString.h>
#include <TFile.h>

Int_t main( Int_t argc, Char_t** argv ){
  if( argc!=3 ){
    std::cerr << "wrong input" << std::endl
	      << " Usage : ./through_cut2 (char*)infile (char*)outdir"
	      << std::endl;
    abort();
  }

  const Char_t* infile = argv[1];
  const Char_t* outdir = argv[2];

  std::string basename = gSystem->BaseName( infile );
  
  if( basename.find(".root") == std::string::npos ){
    std::cerr << "[infile] " << std::setw(80) << infile
	      << " -> wrong file name" << std::endl;
    abort();
  }
  basename.erase( basename.rfind(".root") );
  // --------------------------------------------------------------------
  const Char_t* tname_B   = "h511";
  // --------------------------------------------------------------------
  TChain* chain_B = new TChain( tname_B );
  Int_t  nfile_B  = chain_B->Add( infile );
  // --------------------------------------------------------------------  
  
  MCut_array* cut_ep= new MCut_array( branch_table() );
  make_cut_single_ep( cut_ep, tname_B ); // single lepton(m+) veto
  cut_ep->Set( "Mbc",  0 );
  cut_ep->Set( "de",   0 );
  cut_ep->Set( "pi1_eid", 0 );
  cut_ep->Set( "pi2_eid", 0 );
  cut_ep->Set( "pi3_eid", 0 );
  cut_ep->Set( "pi4_eid", 0 );
  //cut_ep->Display(1);

  MCut_array* cut_em = new MCut_array( branch_table() );
  make_cut_single_em( cut_em, tname_B ); // single lepton(m-) veto
  cut_em->Set( "Mbc",  0 );
  cut_em->Set( "de",   0 );
  cut_em->Set( "pi1_muid", 0 );
  cut_em->Set( "pi2_muid", 0 );
  cut_em->Set( "pi3_muid", 0 );
  cut_em->Set( "pi4_muid", 0 );
  //cut_em->Display(1);

  MCut_array* cut_mp= new MCut_array( branch_table() );
  make_cut_single_mp( cut_mp, tname_B ); // single lepton(m+) veto
  cut_mp->Set( "Mbc",  0 );
  cut_mp->Set( "de",   0 );
  cut_mp->Set( "pi1_muid", 0 );
  cut_mp->Set( "pi2_muid", 0 );
  cut_mp->Set( "pi3_muid", 0 );
  cut_mp->Set( "pi4_muid", 0 );
  //cut_mp->Display(1);

  MCut_array* cut_mm = new MCut_array( branch_table() );
  make_cut_single_mm( cut_mm, tname_B ); // single lepton(m-) veto
  cut_mm->Set( "Mbc",  0 );
  cut_mm->Set( "de",   0 );
  cut_mm->Set( "pi1_muid", 0 );
  cut_mm->Set( "pi2_muid", 0 );
  cut_mm->Set( "pi3_muid", 0 );
  cut_mm->Set( "pi4_muid", 0 );
  //cut_mm->Display(1);

  TCut cut_kfbcl  = "kfbcl>1.0e-18";
  TCut cut_pi_sel_ee = "pi1_eid >0.80 || pi2_eid >0.80 || pi3_eid >0.80 || pi4_eid >0.80 || pi1p<0.40 || pi2p<0.40 || pi3p<0.40 || pi4p<0.40";
  TCut cut_pi_sel_mm = "pi1_muid>0.97 || pi2_muid>0.97 || pi3_muid>0.97 || pi4_muid>0.97 || pi1p<0.80 || pi2p<0.80 || pi3p<0.80 || pi4p<0.80";

  // --------------------------------------------------------------------
  TTree* tree_B_ep = new TTree();
  TTree* tree_B_em = new TTree();
  TTree* tree_B_mp = new TTree();
  TTree* tree_B_mm = new TTree();
 
  tree_B_ep = chain_B->CopyTree( cut_ep->Output() && cut_kfbcl && cut_pi_sel_ee );
  tree_B_em = chain_B->CopyTree( cut_em->Output() && cut_kfbcl && cut_pi_sel_ee );
  tree_B_mp = chain_B->CopyTree( cut_mp->Output() && cut_kfbcl && cut_pi_sel_mm );
  tree_B_mm = chain_B->CopyTree( cut_mm->Output() && cut_kfbcl && cut_pi_sel_mm );

  Int_t entry = 0;
  entry += chain_B->GetEntries( cut_ep->Output()  && cut_kfbcl && cut_pi_sel_ee );
  entry += chain_B->GetEntries( cut_em->Output()  && cut_kfbcl && cut_pi_sel_ee );
  entry += chain_B->GetEntries( cut_mp->Output()  && cut_kfbcl && cut_pi_sel_mm );
  entry += chain_B->GetEntries( cut_mm->Output()  && cut_kfbcl && cut_pi_sel_mm );
  
  TTree* tree_B = new TTree();

  if( entry ){
    std::cout << "[single-lepton-veto(ee or mm)] ";
    tree_B = tree_B_ep->CloneTree(0);
    tree_B = tree_B_em->CloneTree(0);
    tree_B = tree_B_mp->CloneTree(0);
    tree_B = tree_B_mm->CloneTree(0);
    tree_B->CopyEntries( tree_B_ep );
    tree_B->CopyEntries( tree_B_em );
    tree_B->CopyEntries( tree_B_mp );
    tree_B->CopyEntries( tree_B_mm );
  }else{
    std::cerr << "[WARNING] empty entry : " << infile << " -> skip" << std::endl;
    return 0;
  }
  TTree* newtree_B = tree_B->CloneTree();
  std::cout << "[infile] "
	    << std::setw(80) << infile << "(" << tname_B    << ", "
	    << std::setw( 6) << chain_B  ->GetEntries() << " -> "
	    << std::setw( 6) << newtree_B->GetEntries()
	    << "[" << int(100*((double)newtree_B->GetEntries()/chain_B->GetEntries())) << "%]" << ")"
    // << " [check] : " << tree_B_e->GetEntries() << " + " << tree_B_m->GetEntries();
	    << std::endl;
  // --------------------------------------------------------------------
  const Char_t* outfile = Form( "%s%s_single_cut2.root",outdir,basename.c_str() );
  TFile* rootf = new TFile( outfile, "RECREATE" );
  
  newtree_B->Write();
  rootf    ->Close();
  
  delete chain_B;
  delete newtree_B;
  //delete rootf;

  return 0;
}

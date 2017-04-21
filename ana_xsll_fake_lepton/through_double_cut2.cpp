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
	      << " Usage : ./merge_cut2 (char*)infile (char*)outdir"
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
  Int_t nfile_B = chain_B->Add( infile );    

  // --------------------------------------------------------------------  
  MCut_array* cut = new MCut_array( branch_table() );
  make_cut_double_ee( cut, tname_B ); // double lepton(ee) veto
  cut->Set( "Mbc",  0 );
  cut->Set( "de",   0 );
  //cut->Set( "xs_m", 0 );
  //cut->Set(    441, 0 );
  //cut->Set(    443, 0 );
  //cut->Set( 100443, 0 );
  TCut cut_ee = cut->Output();
  make_cut_double_mm( cut, tname_B ); // double lepton(mm) veto
  cut->Set( "Mbc",  0 );
  cut->Set( "de",   0 );
  //cut->Set( "xs_m", 0 );
  //cut->Set(    441, 0 );
  //cut->Set(    443, 0 );
  //cut->Set( 100443, 0 );
  TCut cut_mm = cut->Output();

  TCut cut_kfbcl = "kfbcl>1.0e-18";
  // --------------------------------------------------------------------
  TTree* tree_B = new TTree();
  std::cout << "[double-lepton-veto(ee or mm)] ";
  tree_B = chain_B->CopyTree( (cut_ee || cut_mm) && cut_kfbcl );

  TTree* newtree_B = tree_B->CloneTree();
  std::cout << "[infile] "
	    << std::setw(80) << infile << "(" << tname_B    << ", "
	    << std::setw( 6) << chain_B  ->GetEntries() << " -> "
	    << std::setw( 6) << newtree_B->GetEntries()
	    << "[" << int(100*((double)newtree_B->GetEntries()/chain_B->GetEntries())) << "%]" << ")"
	    << std::endl;
  // --------------------------------------------------------------------
  const Char_t* outfile = Form( "%s%s_double_cut2.root",outdir,basename.c_str() );
  TFile* rootf = new TFile( outfile, "RECREATE" );
  
  newtree_B->Write();
  rootf->Close();
  
  delete chain_B;
  delete cut;
  delete tree_B;
  delete newtree_B;
  delete rootf;

  return 0;
}

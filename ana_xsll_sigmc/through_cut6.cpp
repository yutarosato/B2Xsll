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
			  << " Usage : ./through_cut (char*)infile (char*)outdir" << std::endl,abort();
  
  const Char_t* infile = argv[1];
  const Char_t* outdir = argv[2];
  
  std::string basename = gSystem->BaseName( infile );
  std::string dirname  = gSystem->DirName ( infile );

  if( basename.find(".root") == std::string::npos ){
    std::cerr << "[infile] " << std::setw(80) << infile
	      << " -> wrong file name" << std::endl;
    abort();
  }
  basename.erase( basename.rfind(".root") );
  // --------------------------------------------------------------------  
  const Char_t* tname_B   = "h511";
  const Char_t* tname_Gen = "h12";
  // --------------------------------------------------------------------  
  TChain* chain_B   = new TChain( tname_B   );
  TChain* chain_Gen = new TChain( tname_Gen );
  chain_B  ->Add( infile );
  chain_Gen->Add( infile );
  // --------------------------------------------------------------------  
  MCut_array* cut1 = new MCut_array( branch_table() );
  make_cut_ee( cut1, tname_B );
  cut1->Set( "Mbc",  0 );
  cut1->Set( "de",   0 );
  cut1->Set( "gm_fl_xs*gm_m_xs", 1, 0.0, 0.0, 1.00 );  // added in cut6 !!
  //cut1->Set( "xs_m", 0 );
  //cut1->Set(    441, 0 );
  //cut1->Set(    443, 0 );
  //cut1->Set( 100443, 0 );
  TCut cut_ee = cut1->Output();
  MCut_array* cut2 = new MCut_array( branch_table() );
  make_cut_mm( cut2, tname_B );
  cut2->Set( "Mbc",  0 );
  cut2->Set( "de",   0 );
  cut1->Set( "gm_fl_xs*gm_m_xs", 1, 0.0, 0.0, 1.00 );  // added in cut6 !!
  //cut2->Set( "xs_m", 0 );
  //cut2->Set(    441, 0 );
  //cut2->Set(    443, 0 );
  //cut2->Set( 100443, 0 );
  TCut cut_mm = cut2->Output();

  MCut_array* cut1_gen = new MCut_array( branch_table() );
  make_cut_ee( cut1_gen, tname_Gen );
  TCut cut_ee_gen = cut1_gen->Output();
  MCut_array* cut2_gen = new MCut_array( branch_table() );
  make_cut_mm( cut2_gen, tname_Gen );
  TCut cut_mm_gen = cut2_gen->Output();

  TCut cut_kfbcl = "kfbcl>1.0e-18";
  // --------------------------------------------------------------------
  TTree* tree_B = new TTree;
  std::cout << "[ee or mm] ", tree_B = chain_B->CopyTree( (cut_ee||cut_mm) && cut_kfbcl );
  TTree* newtree_B = tree_B->CloneTree();
  std::cout << "[infile] " << std::setw(80) << infile
	    << "(" << tname_B    << ", " << std::setw(6) << chain_B  ->GetEntries()
	    << " -> "
	    << std::setw(6) << newtree_B  ->GetEntries()
	    << "[" << int(100*((double)newtree_B->GetEntries()/chain_B->GetEntries())) << "%]" << ")"
	    << std::endl;
  TTree* tree_Gen = new TTree;
  tree_Gen = chain_Gen->CopyTree( (cut_ee_gen||cut_mm_gen) );
  TTree* newtree_Gen = tree_Gen->CloneTree();
  //TTree* newtree_Gen = chain_Gen->CloneTree();
  // --------------------------------------------------------------------
  const Char_t* outfile = Form( "%s%s_cut6.root",outdir,basename.c_str() );
  TFile*        rootf   = new TFile( outfile, "RECREATE" );
  newtree_B  ->Write();
  if( strstr(infile, "_m1m_") != NULL || strstr(infile, "_m9999m_") != NULL ) newtree_Gen->Write();
  rootf->Close();
  delete chain_B;
  delete cut1;
  delete cut2;
  delete tree_B;
  delete newtree_B;
  delete newtree_Gen;
  delete rootf;

  return 0;
}

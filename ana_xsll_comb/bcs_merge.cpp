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

Int_t main( Int_t argc, Char_t** argv ){ // exp setname
  if( argc!=5 ){
    std::cerr << "wrong input" << std::endl
	      << " Usage : ./bcs_merge (char*)exp (char*)setname (char*)indir (char*)outdir"
	      << std::endl;
    abort();
  }
  const Char_t* expno      = argv[1];
  const Char_t* fl_setname = argv[2];
  // --------------------------------------------------------------------  
  const Char_t* tname_B = "h511";
  const Char_t* indir   = argv[3];
  const Char_t* outdir  = argv[4];
  // --------------------------------------------------------------------  

  std::stringstream sTmp;
  sTmp << indir << "sigMC_" << expno << "*_set" << fl_setname << "_*.root";
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;


  TChain* chain_B   = new TChain( tname_B );
  Int_t nfile = chain_B->Add(infile);
  if( nfile != 12 ) std::cerr << infile << " --> [ABORT] wrong nfile :" << nfile << std::endl, abort();
  //if( nfile != 16 ) std::cerr << infile << " --> [ABORT] wrong nfile :" << nfile << std::endl, abort(); // for xsspin
  //if( nfile != 4 ) std::cerr << infile << " --> [ABORT] wrong nfile :" << nfile << std::endl, abort(); // for xsjpsi sigMC
  // -------------< DISPLAY > -------------------------------------------
  std::cout << "exp" << std::setw(3) << std::right << expno                 << ", "
	    << "set" << std::setw(3) << std::right << fl_setname            << ", "
	             << std::setw(3) << std::right << nfile                 << " files, "
	             << std::setw(8) << std::right << chain_B->GetEntries() << " events" << std::endl;
  // --------------------------------------------------------------------
  TTree* tree_B = new TTree();
  tree_B = chain_B->CloneTree();
  if( tree_B->GetEntries() != chain_B->GetEntries() ) std::cerr << infile << " --> [ABORT] wrong entries :"
								<< tree_B->GetEntries()
								<< " ?=? "
								<< chain_B->GetEntries()
								<< std::endl, abort();
  
  const Char_t* outfile = Form( "%ssigMC_%s_m9999m_caseB_set%s_bcsmerge.root",outdir,expno, fl_setname );
  //const Char_t* outfile = Form( "%ssigMC_%s_m9999m_caseB_c10flip_set%s_bcsmerge.root",outdir,expno, fl_setname ); // c10-flip

  TFile* rootf = new TFile( outfile, "RECREATE" );
  
  tree_B->Write();
  rootf->Close();

  delete chain_B;
  delete tree_B;
  delete rootf;

  return 0;
}

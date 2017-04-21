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
  if( argc!=6 ){
    std::cerr << "wrong input" << std::endl
	      << " Usage : ./merge_cut2 (char*)indir (char*)outdir (char*)exp (int)stream (char*)type"
	      << std::endl;
    abort();
  }

  const Char_t* indir  = argv[1];
  const Char_t* outdir = argv[2];
  const Char_t* exp    = argv[3];
  const Int_t   stream = atoi( argv[4] );
  const Char_t* type   = argv[5];
  // --------------------------------------------------------------------
  const Bool_t  fl_message = !true;
  const Char_t* tname_B    = "h511";
  Int_t         fl_rd      = -1; // 1(rd)  0(gmc)
  TChain* chain_B = new TChain( tname_B );
  
  std::stringstream sTmp;
  if( strstr(type, "rd")!=NULL ){
    fl_rd = 1; // RD
    sTmp << indir << "/RD_e0" << exp << "*.root";
  }else{
    fl_rd = 0; // gMC
    sTmp << indir << "/gMC_" << type << "_e0" << exp << "*s0" << stream << "*.root";
  }

  Int_t nfile_B = chain_B->Add( sTmp.str().c_str() );
  if( !nfile_B ) std::cerr << "[ABORT] NO file : " << sTmp.str().c_str() << std::endl, abort();
  sTmp.str("");
  sTmp.clear();

  // --------------------------------------------------------------------  
  const Int_t fl_q2 = 0;
  TCut cut_lrnb = (char*)makeCut_LRNB_2d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.91, 0.39, 0.93, 0.91, 0.86, 0.56, 0.92, 0.87).c_str(); // 2d NB_lep(bcs=bb) // after calibration(dzll<190)
  // --------------------------------------------------------------------
  if     ( fl_rd==1 ) sTmp << Form( "%s/RD_e0%s_caseB_double_mergecut2.root",         outdir,      exp         ); // RD
  else if( fl_rd==0 ) sTmp << Form( "%s/gMC_%s_e0%s_s0%d_caseB_double_mergecut2.root",outdir,type, exp, stream ); // gMC
 
  Char_t outfile[1024];
  strcpy( outfile, (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  TFile* rootf  = new TFile( outfile, "RECREATE" );
  TTree* tree_B = new TTree();
  tree_B = chain_B->CopyTree( cut_lrnb );
  TTree* newtree_B = tree_B->CloneTree();


  std::cout << std::setw(7) << std::right << type
	    << ", s0"   << stream
	    << ", exp0" << std::setw(3) << std::left  << exp
	    << ", "     << std::setw(3) << std::right << nfile_B << " files(" << tname_B    << ", " << std::setw(9) << std::right << chain_B->GetEntries()
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

  return 0;
}

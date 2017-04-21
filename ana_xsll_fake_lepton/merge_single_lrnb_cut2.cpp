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
#include <TH2D.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TString.h>
#include <TArrow.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TPaletteAxis.h>

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

  TChain* chain_B = new TChain( tname_B );
  
  std::stringstream sTmp;
  sTmp << indir << "/gMC_" << type << "_e0" << exp << "*s0" << stream << "*.root";

  Int_t nfile_B = chain_B->Add( sTmp.str().c_str() );
  if( !nfile_B ) std::cerr << "[ABORT] NO file : " << sTmp.str().c_str() << std::endl, abort();

  // --------------------------------------------------------------------  
  const Int_t fl_q2 = 0;
  TCut cut_lrnb = (char*)makeCut_LRNB_2d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.91, 0.39, 0.93, 0.91, 0.86, 0.56, 0.92, 0.87).c_str(); // 2d NB_lep(bcs=bb) // after calibration(dzll<190)
  // --------------------------------------------------------------------
  const Char_t* outfile = Form( "%s/gMC_%s_e0%s_s0%d_caseB_single_mergecut2.root",outdir,type, exp, stream );
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

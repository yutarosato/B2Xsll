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
#include "../Set/makeCut.h"

#include "draws_.h"

#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>
#include <TString.h>
#include <TFile.h>

Int_t main( Int_t argc, Char_t** argv ){

  if( argc!=4 ) std::cerr << "wrong input" << std::endl
			  << " Usage : ./through_cut (char*)infile (char*)outdir (int)fl_mode_ll" << std::endl,abort();

  const Char_t* infile     = argv[1];
  const Char_t* outdir     = argv[2];
  const Int_t   fl_mode_ll = atoi(argv[3]);
  
  std::string basename = gSystem->BaseName( infile );  

  if( basename.find(".root") == std::string::npos ){
    std::cerr << "[infile] " << std::setw(80) << infile
	      << " -> wrong file name" << std::endl;
    abort();
  }
  basename.erase( basename.rfind(".root") );
  // --------------------------------------------------------------------  
  const Char_t* tname_B = "h511";
  TChain*       chain_B = new TChain( tname_B );
  chain_B->Add( infile );
  // --------------------------------------------------------------------  
  const Int_t fl_q2 = 0;
  TCut* cut_lrnb = new TCut[2];
  //cut_lrnb[0] = (char*)makeCut_LRNB_2d("nb_lep0_orgksfw_vtxcl_fmiss1_%s", fl_q2, 0.85, 0.55, 0.93, 0.88).c_str(); // for mm // before calib
  //cut_lrnb[1] = (char*)makeCut_LRNB_2d("nb_lep1_orgksfw_vtxcl_fmiss1_%s", fl_q2, 0.91, 0.37, 0.94, 0.92).c_str(); // for ee // before calib
  TCut cut_cc   = "1";
  cut_lrnb[0] = (char*)makeCut_LRNB_2d("nb_lep0_orgksfw_vtxcl_fmiss1_%s", fl_q2, 0.86, 0.56, 0.92, 0.87).c_str(); // for mm // after calib
  cut_lrnb[1] = (char*)makeCut_LRNB_2d("nb_lep1_orgksfw_vtxcl_fmiss1_%s", fl_q2, 0.91, 0.39, 0.93, 0.91).c_str(); // for ee // after calib
  //TCut cut_cc   = "!(lpgt==3 && lmgt==3 && (dzll3dcalib>0.0190 || kfbclcalib<1.0e-18) )"; // charmonium events
  // --------------------------------------------------------------------
  TTree* tree_B = new TTree();
  std::cout << "[lrnb cut] ", tree_B = chain_B->CopyTree( cut_lrnb[fl_mode_ll] && cut_cc );
  TTree* newtree_B;
  if( chain_B->GetEntries(cut_lrnb[fl_mode_ll] && cut_cc) ) newtree_B = tree_B ->CloneTree();
  else                                newtree_B = chain_B->CloneTree(0,"newtree");
  std::cout << "[infile] " << std::setw(80) << infile
	    << "(" << tname_B    << ", " << std::setw(6) << chain_B  ->GetEntries()
	    << " -> "
	    << std::setw(6) << newtree_B  ->GetEntries()
	    << "[" << int(100*((double)newtree_B->GetEntries()/chain_B->GetEntries())) << "%]" << ")"
	    << std::endl;
  // --------------------------------------------------------------------
  const Char_t* outfile = Form( "%s%s_cut_lrnb_emu%d.root",outdir,basename.c_str(),fl_mode_ll );
  TFile*        rootf   = new TFile( outfile, "RECREATE" );
  newtree_B->Write();
  rootf->Close();

  delete   chain_B;
  delete   tree_B;
  delete   rootf;
  delete[] cut_lrnb;

  return 0;
}

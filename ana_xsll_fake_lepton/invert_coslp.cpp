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
  else if( basename.find("RD")    != std::string::npos ) fl_sbr = 2;
  else std::cerr << "[ABORT] wrong file name : " << infile << std::endl, abort();

  basename.erase( basename.rfind(".root") );
  // --------------------------------------------------------------------  
  const Char_t* tname_B       = "h511";
  // --------------------------------------------------------------------  
  TChain* chain_B   = new TChain( tname_B   );
  chain_B  ->Add( infile );
  // --------------------------------------------------------------------  
  Float_t coslp, recbfl;
  
  chain_B->SetBranchAddress( "coslp",  &coslp  );
  chain_B->SetBranchAddress( "recbfl", &recbfl );
  
  TTree* newtree_B = chain_B->CloneTree(0);

  Int_t  nevt      = 0;
  while( chain_B->GetEntry(nevt, 0) ) {
    if( fabs(recbfl)==1 ) coslp = -coslp;
    newtree_B->Fill();
    nevt++;
  }

  std::cout << "[infile] " << std::setw(80) << infile
	    << "(" << tname_B    << ", " << std::setw(6) << chain_B  ->GetEntries() << ")"
	    << std::endl;

  // --------------------------------------------------------------------
  const Char_t* outfile = Form( "%s%s_inv.root",outdir,basename.c_str() );
  TFile*        rootf   = new TFile( outfile, "RECREATE" );
  newtree_B->Write();
  rootf->Close();

  delete chain_B;
  delete newtree_B;
  delete rootf;

  return 0;
}

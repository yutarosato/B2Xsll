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
  if( argc!=6 ){
    std::cerr << "wrong input" << std::endl
	      << " Usage : ./merge_cut3 (char*)indir (char*)outdir (char*)exp (int)stream (char*)type"
	      << std::endl;
    abort();
  }

  const Char_t* indir  = argv[1];
  const Char_t* outdir = argv[2];
  const Char_t* exp    = argv[3];
  const Int_t   stream = atoi( argv[4] );
  const Char_t* type   = argv[5];
  // --------------------------------------------------------------------
  const Char_t* tname_B   = "h511";
  // --------------------------------------------------------------------
  Int_t fl_rd  = -1; // 1(rd)  0(gmc)
  Int_t fl_emu = -1; // 1(emu) 0(ee or mm)
  TChain* chain_B = new TChain( tname_B );
  
  std::stringstream sTmp;
  if( strstr(type, "rd")!=NULL ){
    fl_rd = 1; // RD
    sTmp << indir << "/" << type << "/9999/RD_e0" << exp << "*.root";
  }else{
    fl_rd = 0; // gMC
    sTmp << indir << "/" << type << "/9999/gMC_" << type << "_e0" << exp << "*s0" << stream << "*.root";
    //sTmp << indir << "/gMC_" << type << "_e0" << exp << "*s0" << stream << "*.root"; // tmppppp
  }
  Int_t nfile_B = chain_B->Add( sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();

  // --------------------------------------------------------------------  
  MCut_array* cut1 = new MCut_array( branch_table() );
  make_cut_ee( cut1, tname_B );
  cut1->Set( "Mbc",     0 );
  cut1->Set( "de",      0 );
  //cut1->Set( "xs_m",    0 );
  cut1->Set(    441,    0 );
  cut1->SetFunc(443,    1 );
  cut1->Set( "cc_morg", 0 );
  cut1->Set( 100443,    0 );
  cut1->Set( "dzll3d",      0         ); // for calib
  cut1->Set( "dzll3dcalib", 1, 0.0190 ); // for calib
  TCut cut_ee = cut1->Output();

  MCut_array* cut2 = new MCut_array( branch_table() );
  make_cut_mm( cut2, tname_B );
  cut2->Set( "Mbc",     0 );
  cut2->Set( "de",      0 );
  //cut2->Set( "xs_m",    0 );
  cut2->Set(    441,    0 );
  cut2->SetFunc(443,    1 );
  cut2->Set( "cc_morg", 0 );
  cut2->Set( 100443,    0 );
  cut2->Set( "dzll3d",      0         ); // for calib
  cut2->Set( "dzll3dcalib", 1, 0.0190 ); // for calib
  TCut cut_mm = cut2->Output();
  
  MCut_array* cut3 = new MCut_array( branch_table() );
  make_cut_emu( cut3, tname_B );
  cut3->Set( "Mbc",     0 );
  cut3->Set( "de",      0 );
  //cut2->Set( "xs_m",    0 );
  cut3->Set(    441,    0 );
  cut3->SetFunc(443,    1 );
  cut3->Set( "cc_morg", 0 );
  cut3->Set( 100443,    0 );
  cut3->Set( "dzll3d",      0         ); // for calib
  cut3->Set( "dzll3dcalib", 1, 0.0190 ); // for calib
  TCut cut_emu = cut3->Output();

  //TCut cut_kfbcl = "kfbcl>1.0e-18";
  TCut cut_kfbcl = "kfbclcalib>1.0e-18"; // for calib
  // --------------------------------------------------------------------
  TTree* tree_B     = new TTree;
  if( strstr(indir, "emu")!=NULL ){
    fl_emu = 1; //emu
    //std::cout << "[  emu   ] ";
    tree_B = chain_B->CopyTree( cut_emu && cut_kfbcl );
  }else{
    fl_emu = 0; // ee or mm
    //std::cout << "[ee or mm] ";
    tree_B = chain_B->CopyTree( (cut_ee || cut_mm) && cut_kfbcl );
      }

  TTree* newtree_B = tree_B->CloneTree();


  // ------------------- display ----------------------------------------
  std::cout << Form( "fl_rd(%d), fl_emu(%d) : ", fl_rd, fl_emu)
	    << std::setw(7) << std::right << type
	    << ", s0"       << stream
	    << ", exp0"     << std::setw(3) << std::left  << exp
	    << ", "         << std::setw(3) << std::right << nfile_B << " files(" << tname_B    << ", " << std::setw(9) << std::right << chain_B->GetEntries()
	    << " -> "
	    << std::setw(9) << std::right << newtree_B  ->GetEntries()
	    << "[" << std::setw(2) << std::right << int(100*((double)newtree_B->GetEntries()/chain_B->GetEntries())) << "%]"
	    << ")"
	    << std::endl;
  // --------------------------------------------------------------------
  if     ( fl_rd && fl_emu ) sTmp << Form( "%s/RD_e0%s_caseB_emu_mergecut3.root",         outdir,       exp         ); // [RD,    emu    ]
  else if( fl_rd           ) sTmp << Form( "%s/RD_e0%s_caseB_mergecut3.root",             outdir,       exp         ); // [RD,  ee or mm ]
  else if(          fl_emu ) sTmp << Form( "%s/gMC_%s_e0%s_s0%d_caseB_emu_mergecut3.root",outdir, type, exp, stream ); // [gMC,   emu    ]
  else                       sTmp << Form( "%s/gMC_%s_e0%s_s0%d_caseB_mergecut3.root",    outdir, type, exp, stream ); // [gMC, ee or mm ]

  Char_t outfile[1024];
  strcpy( outfile, (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  TFile* rootf = new TFile( outfile, "RECREATE" );


  newtree_B  ->Write();
  rootf->Close();
  
  delete chain_B;
  delete cut1;
  delete cut2;
  delete cut3;
  delete tree_B;
  delete newtree_B;
  delete rootf;

  return 0;
}

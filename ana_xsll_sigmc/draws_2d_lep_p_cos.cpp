#include <iostream>

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

#include <vector>
#include <stdlib.h>
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TArrow.h>
#include <TFile.h>


Int_t main( Int_t argc, Char_t** argv ){
  using namespace sigmc;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( !(argc==3 || argc==4) ) std::cerr << "wrong input" << std::endl
					  << " Usage : ./draws_2d (int)fl_mode_ll (char*)setname [(int)fl_appRun]" << std::endl
					  << "[fl_mode_ll] : 1(e), 0(mu)"           << std::endl, abort();

  Int_t   fl_mode_ll = atoi( argv[1] ); // 1(e), 0(mu)
  Char_t* setname    = argv[2]; 
  Int_t   fl_appRun  = 1;
  if( argc==4 ) fl_appRun = atoi( argv[3] );

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 1; // fixed
  const Int_t Nhist         = 1; // fixed
  const Int_t nfile[Nchain] = {0};
  
  std::stringstream sTmp;
  sTmp << indir << "sigMC_*_m9999m_caseB_set[" << setname << "]";
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;

  TCut* add_cut = new TCut[Nchain];
  add_cut[0] = "self==1";
  
  const Bool_t  flag_save  = true;

  const Char_t*  tname    = "h511";
  const Char_t*  axis[2]  = {"mpp:lpc","epp:lpc"}; // mu, e
  const Int_t    ybin     =    10; // 10 or 15
  const Double_t ymin     =  0.40;
  const Double_t ymax     =  6.40;
  const Char_t*  ylabel   = "P [GeV]";
  const Double_t offset   =   0.0;
  const Int_t    xbin     =    10;
  const Double_t xmin     =  -1.0;
  const Double_t xmax     =   1.0;
  const Char_t*  xlabel   = "cos#theta";
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH2D**    tmphist = new TH2D*  [Nchain];
  TH2D**    hist    = new TH2D*  [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",1 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++

  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis[fl_mode_ll]) << std::endl;
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile, tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  //
  //chain[0]->GetCut()->Set(    441, 0 );
  //chain[0]->GetCut()->Set(    443, 0 );
  //chain[0]->GetCut()->Set( 100443, 0 );
  //
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j ) << std::endl;
    chain[j]->GetCut()->Display(0);
    chain[j]->MakeTree();
  }
  
  //chain[0]->GetTree()->Scan("cc_m:cc_morg:coslp:lppt:lmpt:mpp:mmp:epp:emp:lpc:lmc:lpmoid:lmmoid", "cc_m*cc_m<1.0 && (coslp>0.6||coslp<-0.6)" );
  //chain[0]->GetTree()->Draw( "bc:lmc", "self==1 && (coslp>0.6 || coslp<-0.6) && cc_m*cc_m<1.0" );
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin,offset+xmin,offset+xmax, ybin,ymin,ymax );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis[fl_mode_ll], add_cut[j] );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nchain; i++ ){
    hist[i] = new TH2D( Form("lep%d",fl_mode_ll), Form("lep%d",fl_mode_ll), xbin,offset+xmin,offset+xmax, ybin,ymin,ymax );
    ((TGaxis*)hist[i]->GetXaxis())->SetMaxDigits(3);
    ((TGaxis*)hist[i]->GetYaxis())->SetMaxDigits(3);
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    hist[i]->SetXTitle( xlabel );
    hist[i]->SetYTitle( ylabel );
    hist[i]->Add( tmphist[i] );
    // +++++++ display ++++++++++++++++++++++++++++++++++
    std::cout << std::setw(10) << Form("<hist%d>",i);
    std::cout << std::setw(9)  << hist[i]->GetEntries() << " events"
	      << std::endl; 
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  for( Int_t i=0; i<Nhist; i++ ) hist[i]->Draw( "COLZ" );
    // +++++++ tlegend ++++++++++++++++++++++++++++++++++
  c1->Update();
  if( flag_save ){
    c1->Print(     Form("pic/2d_lep_p_cos_lep%d_set%s.eps",  fl_mode_ll, setname) );
    TFile outfile( Form("pic/2d_lep_p_cos_lep%d_set%s.root", fl_mode_ll, setname), "RECREATE" );
    for( Int_t i=0; i<Nhist; i++ ) hist[i]->Write();
    outfile.Close();
  }
        
  std::cout << "finish" << std::flush;  
  if( fl_appRun ) app.Run();
  
  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete   c1;
    
  return 0;
}

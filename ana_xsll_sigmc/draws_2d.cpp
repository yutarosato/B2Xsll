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


Int_t main( Int_t argc, Char_t** argv ){
  using namespace sigmc;
  TApplication app( "app", &argc, argv );
  Style();

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 1;
  const Int_t nfile[Nchain] = {0};
  const Int_t fl_mode_ll    = 1;
  
  std::stringstream sTmp;
  sTmp << indir << "sigMC_*_caseB";
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;

  TCut* add_cut = new TCut[Nchain];
  add_cut[0] = "";
  const Bool_t  flag_save  = true;
  const Char_t* outfile    = "pic/2d_Mbc_deltaE.eps";
  
  using namespace Mbc_deltaE;
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  MChain**  chain   = new MChain*[Nchain];
  TH2D**    tmphist = new TH2D*  [Nchain];
  TCanvas*  c1      = Canvas( "c1","c1",2 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile, tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  //
  for( Int_t j=0; j<Nchain; j++ ){
    chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
    chain[j]->GetCut()->Set( "de",  1, -0.20, 0.0, 0.20 );
  }
  //
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j ) << std::endl;
    chain[j]->GetCut()->Display(0);
    chain[j]->MakeTree();
  }

  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin,offset+xmin,offset+xmax, ybin,ymin,ymax );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  
  TH2D* hist = new TH2D( Form("%s",axis),Form("%s",axis), xbin,offset+xmin,offset+xmax, ybin,ymin,ymax );
  ((TGaxis*)hist->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)hist->GetYaxis())->SetMaxDigits(3);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->SetXTitle( xlabel );
  hist->SetYTitle( ylabel );
  for( Int_t j=0; j<Nchain; j++ ){
    hist->Add( tmphist[j] );
  }
  
  // +++++++ display ++++++++++++++++++++++++++++++++++
  
  std::cout << std::setw(10) << Form("<hist>" );
  std::cout << std::setw(9)  << hist->GetEntries() << " events( "
	    << std::setw(7)  << hist->Integral()   << " ) "
	    << std::endl;
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  
  c1->Draw();
  c1->cd(1);
  hist->Draw( "COLZ" );
  
  // +++++++ tlegend ++++++++++++++++++++++++++++++++++
  c1->cd(2);
  TPaveText* box    = new TPaveText( 0.0,0.4,1.0,1.0 );
  box->SetTextSize(0.015);
  for( Int_t j=0; j<Nchain; j++ ){
    box->AddText( Form("tmphist%d : Mode(%s : %s)",j,chain[j]->GetModeName(),gSystem->BaseName(chain[j]->GetFileName())) );
    box->AddText( Form("   %s",gSystem->DirName(chain[j]->GetFileName())) );
    box->AddText( Form("   %s",tmphist[j]->GetTitle()) );
  }
  
  box->Draw();
  
  c1->Update();
  if( flag_save ) c1->Print( outfile );
  std::cout << std::flush;  
  app.Run();
  
  delete[] chain;
  delete[] tmphist;
  delete   hist;
  delete   c1;
    
  return 0;
}

#include <iostream>
#include <math.h>

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
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TF1.h>
#include <TArrow.h>

Int_t main( Int_t argc, Char_t** argv ){
  using namespace gmc;  
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==3 || argc==4) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_Mbc (char*)stream (double)used_nstream [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Char_t*  stream       = argv[1];
  Double_t used_nstream = atof(argv[2]); 
  Int_t    fl_appRun    = 1;
  if( argc==4 ) fl_appRun = atoi( argv[3] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 2*nbkgtype;
  const Int_t Nhist         = 3*nbkgtype;
  const Int_t nfile[Nchain] = {0};
  const Int_t fl_q2         = 0;
  const Int_t fl_mode_ll[Nchain] = {1,1,1,1, // 1(e)
				    0,0,0,0, // 0(mu)
  };
  const Double_t scale_event_bkg = used_nstream; // gmc : N -> N/alpha
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    sTmp << indir << "/gMC_" << bkgtype[i%nbkgtype] << "_*_s0[" << stream << "]";
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1},       // uds    (ee)
    {1,1},     // charm  (ee)
    {1,1,1},   // mixed  (ee)
    {1,1,1,1}, // charged(ee)
    {0,0,0,0,1},       // uds    (mm)
    {0,0,0,0,1,1},     // charm  (mm)
    {0,0,0,0,1,1,1},   // mixed  (mm)
    {0,0,0,0,1,1,1,1}, // charged(mm)
    {1,0,0,0,1},         // uds    (ee+mm)
    {1,1,0,0,1,1},       // charm  (ee+mm)
    {1,1,1,0,1,1,1},     // mixed  (ee+mm)
    {1,1,1,1,1,1,1,1},   // charged(ee+mm)
  };

  Char_t* LRNB_cut = new Char_t[2048];
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.91, 0.36, 0.94, 0.92, 0.87, 0.57, 0.92, 0.89).c_str() ); // 2d NB_lep(bcs=bb)
  LRNB_cut = "1"; // LR or NB cut should be applied 
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    //sTmp << makeCut_5body_veto().c_str();
    //sTmp << " && ";
    //sTmp << makeCut_unflavor_veto().c_str();
    //sTmp << " && ";
    sTmp << LRNB_cut;
    strcpy( add_cut[i], (Char_t*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  using namespace Mbc_bkg_wide;
  //using namespace Mbc_comb;
  //using namespace Mbc;

  const Bool_t flag_save  = true; // outfile.eps and outfile.root

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1", 3, 1 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll[j] )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 );
    //chain[j]->GetCut()->Set( "de",  1, -0.20, 0.0, -0.15 );
  }
  //
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j ) << std::endl;
    chain[j]->GetCut()->Display(0);
    chain[j]->MakeTree();
    std::cout << "add_cut : " << add_cut[j] << std::endl;
  }
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin,offset+xmin,offset+xmax );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
    tmphist[j]->Scale( 1/scale_event_bkg );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    Deco( hist[i], 3, col_fil[i%nbkgtype], col_fil[i%nbkgtype] );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }

  // +++++++ display ++++++++++++++++++++++++++++++++++
  std::cout << "=================================" << std::endl
	    << "[Scale Factor]"                    << std::endl
	    << "bkg : " << scale_event_bkg         << std::endl
	    << "=================================" << std::endl;
  Double_t entry_sig[Nhist] = {0}; // # of events in signal box region
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
    entry_sig[i] = hist[i]->Integral( hist[i]->FindBin(5.27+0.0000001), xbin );
    TString sTmp = "Added-Files( ";
    for( Int_t j=0; j<Nchain; j++ ) if( add[i][j] ) sTmp += j, sTmp += ",";
    sTmp += ")[";
    sTmp += entry_all;
    sTmp += " events ( canvas : ";
    sTmp += entry_canvas;
    sTmp += " / under : ";
    sTmp += entry_under;
    sTmp += " / over  : ";
    sTmp += entry_over;
    sTmp += " / sig  : ";
    sTmp += entry_sig[i];
    sTmp += "]";
    
    hist[i]->SetTitle( sTmp.Data() );
    std::cout << sTmp.Data() << std::endl;
  }

  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
	    << std::setw(12) << std::right << "event-type"
	    << " || "
	    << std::setw(12) << std::right << "ee" 
	    << std::setw(12) << std::right << "mm" 
	    << std::setw(12) << std::right << "ee+mm" << std::endl
	    << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  for( Int_t i=0; i<Nchain/2; i++ ){
    if     ( i== 0 ) std::cout << std::setw(12) << std::right << "uds";
    else if( i== 1 ) std::cout << std::setw(12) << std::right << "charm";
    else if( i== 2 ) std::cout << std::setw(12) << std::right << "mixed";
    else if( i== 3 ) std::cout << std::setw(12) << std::right << "charged";
    std::cout << " || "
	      << std::setw(12) << std::right << tmphist[         i]->Integral( tmphist[i]->FindBin(5.27+0.00000001),tmphist[i]->FindBin(5.29-0.000000001) ) // ee
	      << std::setw(12) << std::right << tmphist[nbkgtype+i]->Integral( tmphist[i]->FindBin(5.27+0.00000001),tmphist[i]->FindBin(5.29-0.000000001) ) // mm
	      << std::setw(12) << std::right
	      << tmphist[i]->Integral( tmphist[i]->FindBin(5.27+0.00000001),tmphist[i]->FindBin(5.29-0.000000001) ) + tmphist[nbkgtype+i]->Integral( tmphist[i]->FindBin(5.27+0.00000001),tmphist[i]->FindBin(5.29-0.000000001) ) // ee+mm
	      << std::endl;
  }
  std::cout << std::setw(12) << std::right << "qq"
	    << " || "
	    << std::setw(12) << std::right << tmphist[0]->Integral( tmphist[0]->FindBin(5.27+0.00000001),tmphist[0]->FindBin(5.29-0.000000001) ) + tmphist[1]->Integral( tmphist[1]->FindBin(5.27+0.00000001),tmphist[1]->FindBin(5.29-0.000000001) )
	    << std::setw(12) << std::right << tmphist[4]->Integral( tmphist[4]->FindBin(5.27+0.00000001),tmphist[4]->FindBin(5.29-0.000000001) ) + tmphist[5]->Integral( tmphist[5]->FindBin(5.27+0.00000001),tmphist[5]->FindBin(5.29-0.000000001) )
	    << std::setw(12) << std::right
	    << tmphist[0]->Integral( tmphist[0]->FindBin(5.27+0.00000001),tmphist[0]->FindBin(5.29-0.000000001) ) + tmphist[1]->Integral( tmphist[1]->FindBin(5.27+0.00000001),tmphist[1]->FindBin(5.29-0.000000001) ) + tmphist[4]->Integral( tmphist[4]->FindBin(5.27+0.00000001),tmphist[4]->FindBin(5.29-0.000000001) ) + tmphist[5]->Integral( tmphist[5]->FindBin(5.27+0.00000001),tmphist[5]->FindBin(5.29-0.000000001) )
	    << std::endl;
  std::cout << std::setw(12) << std::right << "bb"
	    << " || "
	    << std::setw(12) << std::right << tmphist[2]->Integral( tmphist[2]->FindBin(5.27+0.00000001),tmphist[2]->FindBin(5.29-0.000000001) ) + tmphist[3]->Integral( tmphist[3]->FindBin(5.27+0.00000001),tmphist[3]->FindBin(5.29-0.000000001) )
	    << std::setw(12) << std::right << tmphist[6]->Integral( tmphist[6]->FindBin(5.27+0.00000001),tmphist[6]->FindBin(5.29-0.000000001) ) + tmphist[7]->Integral( tmphist[7]->FindBin(5.27+0.00000001),tmphist[7]->FindBin(5.29-0.000000001) )
	    << std::setw(12) << std::right
	    << tmphist[2]->Integral( tmphist[2]->FindBin(5.27+0.00000001),tmphist[2]->FindBin(5.29-0.000000001) ) + tmphist[3]->Integral( tmphist[3]->FindBin(5.27+0.00000001),tmphist[3]->FindBin(5.29-0.000000001) ) + tmphist[6]->Integral( tmphist[6]->FindBin(5.27+0.00000001),tmphist[6]->FindBin(5.29-0.000000001) ) + tmphist[7]->Integral( tmphist[7]->FindBin(5.27+0.00000001),tmphist[7]->FindBin(5.29-0.000000001) )
	    << std::endl;
  std::cout << std::setw(12) << std::right << "total"
	    << " || "
	    << std::setw(12) << std::right << tmphist[0]->Integral( tmphist[0]->FindBin(5.27+0.00000001),tmphist[0]->FindBin(5.29-0.000000001) ) + tmphist[1]->Integral( tmphist[1]->FindBin(5.27+0.00000001),tmphist[1]->FindBin(5.29-0.000000001) ) + tmphist[2]->Integral( tmphist[2]->FindBin(5.27+0.00000001),tmphist[2]->FindBin(5.29-0.000000001) ) + tmphist[3]->Integral( tmphist[3]->FindBin(5.27+0.00000001),tmphist[3]->FindBin(5.29-0.000000001) )
	    << std::setw(12) << std::right << tmphist[4]->Integral( tmphist[4]->FindBin(5.27+0.00000001),tmphist[4]->FindBin(5.29-0.000000001) ) + tmphist[5]->Integral( tmphist[5]->FindBin(5.27+0.00000001),tmphist[5]->FindBin(5.29-0.000000001) ) + tmphist[6]->Integral( tmphist[6]->FindBin(5.27+0.00000001),tmphist[6]->FindBin(5.29-0.000000001) ) + tmphist[7]->Integral( tmphist[7]->FindBin(5.27+0.00000001),tmphist[7]->FindBin(5.29-0.000000001) )
	    << std::setw(12) << std::right
	    << tmphist[0]->Integral( tmphist[0]->FindBin(5.27+0.00000001),tmphist[0]->FindBin(5.29-0.000000001) ) + tmphist[1]->Integral( tmphist[1]->FindBin(5.27+0.00000001),tmphist[1]->FindBin(5.29-0.000000001) ) + tmphist[2]->Integral( tmphist[2]->FindBin(5.27+0.00000001),tmphist[2]->FindBin(5.29-0.000000001) ) + tmphist[3]->Integral( tmphist[3]->FindBin(5.27+0.00000001),tmphist[3]->FindBin(5.29-0.000000001) ) + tmphist[4]->Integral( tmphist[4]->FindBin(5.27+0.00000001),tmphist[4]->FindBin(5.29-0.000000001) ) + tmphist[5]->Integral( tmphist[5]->FindBin(5.27+0.00000001),tmphist[5]->FindBin(5.29-0.000000001) ) + tmphist[6]->Integral( tmphist[6]->FindBin(5.27+0.00000001),tmphist[6]->FindBin(5.29-0.000000001) ) + tmphist[7]->Integral( tmphist[7]->FindBin(5.27+0.00000001),tmphist[7]->FindBin(5.29-0.000000001) )
	    << std::endl;
  
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

  
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();

  TH2D* waku = Waku( Nhist, hist, xlabel );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s (stream=%s)", axis, stream) );

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend** legend1 = new TLegend*[3];
  for( Int_t i=0; i<3; i++ ){
    legend1[i] = new TLegend( 0.75,0.75,0.99,0.99 );
    if     ( i== 0 ) legend1[i]->SetHeader("ee"       );
    else if( i== 1 ) legend1[i]->SetHeader("#mu#mu"   );
    else if( i== 2 ) legend1[i]->SetHeader("ee+#mu#mu");
    legend1[i]->AddEntry( hist[4*i+3],"charged", "P" );
    legend1[i]->AddEntry( hist[4*i+2]," mixed",  "P" );
    legend1[i]->AddEntry( hist[4*i+1]," charm",  "P" );
    legend1[i]->AddEntry( hist[4*i+0],"  uds",   "P" );
    c1->cd(i+1);
    waku->Draw();
    legend1[i]->Draw();
    for(Int_t j=nbkgtype-1; j>=0; j-- ) hist[nbkgtype*i+j]->Draw("same");
  }
  
  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_type_s0%s.eps",  axis, stream) );
    c1->Print( Form("pic/%s_type_s0%s.root", axis, stream) );
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] func;
  delete[] add_cut;
  delete   c1;
    
  return 0;
}


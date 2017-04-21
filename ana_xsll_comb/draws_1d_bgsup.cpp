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
#include <TFile.h>
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
  using namespace sig_bkg;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage   : ./draws (char*)stream (char*)setname (int)fl_axis [(int)fl_appRun]" << std::endl
					<< " Example : ./draws          01            A-B          0" << std::endl
					<< "[fl_axis] 0(evis) 1(mmiss) 2(deroe) 3(kfbchi) 4(dzll) 5(cos-theta_B) 6(Fmiss) 7(ksfw) 8(Fmiss_qq) 9(Fmiss_bb)"
					<< std::endl, abort();
  Char_t* stream     = argv[1];
  Char_t* setname    = argv[2];
  Int_t   fl_axis    = atoi( argv[3] );
  Int_t   fl_appRun  = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 10;    // sig(e,mu), uds(e,mu), charm(e,mu), mixed(e,mu), charged(e,mu)
  const Int_t Nhist         = 4;     // sig, bkg(qq), bkg(bb), bkg(total)
  const Int_t fl_sb[Nchain] = {1,1}; // 1(sig), 0(bkg)
  const Int_t nfile[Nchain] = {0};
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if( fl_sb[i] ) sTmp << indir[1] << "sigMC_*_m9999m_caseB_set[" << setname << "]";         // sig
    else           sTmp << indir[0] << "/gMC_" << bkgtype[(i-2)/2] << "_*_s0[" << stream << "]";  // bkg
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1,1},                 // sig
    {0,0,1,1,1,1},         // bkg(qq)
    {0,0,0,0,0,0,1,1,1,1}, // bkg(bb)
    {0,0,1,1,1,1,1,1,1,1}, // bkg(total)
  };
  const Int_t fl_mode_ll[Nchain] = {1,0,1,0,1,0,1,0,1,0};

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];
  strcpy( add_cut[0], makeCut_category(1).c_str() ); // true
  strcpy( add_cut[1], makeCut_category(1).c_str() ); // true
  add_cut[2] = "";
  add_cut[3] = "";
  add_cut[4] = "";
  add_cut[5] = "";
  add_cut[6] = "lpgt!=3 || lmgt!=3";
  add_cut[7] = "lpgt!=3 || lmgt!=3";
  add_cut[8] = "lpgt!=3 || lmgt!=3";
  add_cut[9] = "lpgt!=3 || lmgt!=3";
  
  //using namespace Mbc;
  //using namespace deltaE;
  //using namespace M_Xs;
  //using namespace M_ll;
  using namespace bgsup;

  const Bool_t  flag_fit     = true;
  const Int_t   sel_fun[2][8] = {{0,0,0,1092,3,0, 10,10}, // bkg
				 {0,0,0, 350,5,2,410,10}, // sig
  };
  
  const Bool_t flag_save  = true; // outfile.eps and outfile.root
  const Bool_t flag_scale = true;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",2 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis[fl_axis]) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll[j] )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
    //chain[j]->GetCut()->Set( "de",  1, -0.20, 0.0, 0.20 );
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
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis[fl_axis], add_cut[j] );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  hist[0] = new TH1D( Form("pdf_hist_%s_sig",    fname[fl_axis]), Form("pdf_hist_%s_set%s_sig",    fname[fl_axis],setname), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[1] = new TH1D( Form("pdf_hist_%s_qq_bkg", fname[fl_axis]), Form("pdf_hist_%s_qq_s0%s_bkg",  fname[fl_axis],stream),  xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[2] = new TH1D( Form("pdf_hist_%s_bb_bkg", fname[fl_axis]), Form("pdf_hist_%s_bb_s0%s_bkg",  fname[fl_axis],stream),  xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[3] = new TH1D( Form("pdf_hist_%s_tot_bkg",fname[fl_axis]), Form("pdf_hist_%s_tot_s0%s_bkg", fname[fl_axis],stream),  xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  for( Int_t i=0; i<Nhist; i++ ){
    //hist[i] = new TH1D( Form("pdf_hist_%s_%s_s0%s_bkg",fname[fl_axis],event_type,stream),Form("pdf_hist_%s_%s_s0%s_bkg",fname[fl_axis],event_type,stream), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
    Deco( hist[i], 2, i+1, i+1);
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }

  // +++++++ display ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin[fl_axis]+1);
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
    sTmp += "]";
    
    hist[i]->SetTitle( sTmp.Data() );
    std::cout << sTmp.Data() << std::endl;
    if( flag_scale ){
      hist[i]->Sumw2();
      hist[i]->Scale( 1/(hist[i]->GetEntries()) );
    }
  }

  // +++++++ make fit-func ++++++++++++++++++++++++++++++++++
  if(flag_fit){
    func[0] = new TF1( Form("pdf_func_%s_sig",     fname[fl_axis]), make_func(sel_fun[1][fl_axis]),offset[fl_axis]+xmin_fit[fl_axis],offset[fl_axis]+xmax_fit[fl_axis],n_fitfunc_par(sel_fun[1][fl_axis]) );
    func[1] = new TF1( Form("pdf_func_%s_qq_bkg",  fname[fl_axis]), make_func(sel_fun[0][fl_axis]),offset[fl_axis]+xmin_fit[fl_axis],offset[fl_axis]+xmax_fit[fl_axis],n_fitfunc_par(sel_fun[0][fl_axis]) );
    func[2] = new TF1( Form("pdf_func_%s_bb_bkg",  fname[fl_axis]), make_func(sel_fun[0][fl_axis]),offset[fl_axis]+xmin_fit[fl_axis],offset[fl_axis]+xmax_fit[fl_axis],n_fitfunc_par(sel_fun[0][fl_axis]) );
    func[3] = new TF1( Form("pdf_func_%s_tot_bkg", fname[fl_axis]), make_func(sel_fun[0][fl_axis]),offset[fl_axis]+xmin_fit[fl_axis],offset[fl_axis]+xmax_fit[fl_axis],n_fitfunc_par(sel_fun[0][fl_axis]) );
    for( Int_t i=0; i<Nhist; i++ ){
      //func[i] = new TF1( Form("pdf_func_%s_%s_s0%s_bkg",fname[fl_axis],event_type,stream), make_func(sel_fun[fl_axis]),offset[fl_axis]+xmin_fit[fl_axis],offset[fl_axis]+xmax_fit[fl_axis],n_fitfunc_par(sel_fun[fl_axis]) );
      func[i]->SetLineColor(i+1);
      if( sel_fun[fl_sb[i]][fl_axis]==1092 ){
	func[i]->SetParNames( "area1","alpha1","n1","area2","alpha2", "n2" );
	func[i]->SetParameters(  0.019,  0.36,  1.9, 0.11,    2.0,    3.1  );
      }else{
	func_set_parameters(sel_fun[fl_sb[i]][fl_axis], func[i], hist[i], xbin[fl_axis], offset[fl_axis]+xmin[fl_axis], offset[fl_axis]+xmax[fl_axis]);
      }
    }
  }
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);
  TH2D* waku = Waku( Nhist, hist, xlabel[fl_axis] );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s", xlabel[fl_axis]) );
  waku->Draw();

  if( flag_fit ){
    for(Int_t i=0; i<Nhist; i++ ){
      hist[i]->Fit(func[i],"R","PE0same");
      std::cout << "chi2/NDF = " << func[i]->GetChisquare() << " / " << func[i]->GetNDF()
		<< " = " << func[i]->GetChisquare()/func[i]->GetNDF()
		<< std::endl;
    }
  }else{
    for( Int_t i=0; i<Nhist; i++ ) hist[i]->Draw( "same" );
  }
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[0],"sig",      "L" );
  legend1->AddEntry( hist[1],"bkg(qq )", "L" );
  legend1->AddEntry( hist[2],"bkg(bb )", "L" );
  legend1->AddEntry( hist[3],"bkg(tot)", "L" );
  legend1->Draw();
  // +++++++ tlegend ++++++++++++++++++++++++++++++++++
  c1->cd(2);
  TPaveText* box    = new TPaveText( 0.0,0.4,1.0,1.0 );
  box->SetTextSize(0.015);
  TLegend*   legend = new TLegend  ( 0.0,0.0,1.0,0.4 );
  legend->SetTextSize(0.015);
  for( Int_t j=0; j<Nchain; j++ ){
    box->AddText( Form("tmphist%d : Mode(%s : %s : %d files)", j,chain[j]->GetModeName(),gSystem->BaseName(chain[j]->GetFileName()),chain[j]->GetNfile()) );
    box->AddText( Form("   %s", add_cut[j]) );
  }

  for( Int_t i=0; i<Nhist; i++ ){
    legend->AddEntry( hist[i],"","PL" );
  }
  box->Draw();
  legend->Draw();

  c1->Update();
  if( flag_save ){
    c1->Print(     Form("pic/pdf_%s_s0%s_set%s.eps",  fname[fl_axis], stream, setname)             );
    TFile outfile( Form("pic/pdf_%s_s0%s_set%s.root", fname[fl_axis], stream, setname), "RECREATE" );
    c1->Write();
    for( Int_t i=0; i<Nhist; i++ ) hist[i]->Write();
    if( flag_fit ){
      for( Int_t i=0; i<Nhist; i++ ) func[i]->Write();
    }
    outfile.Close();
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] func;
  delete[] add_cut;
  delete   c1;
  delete   legend;
    
  return 0;
}

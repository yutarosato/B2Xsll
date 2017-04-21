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
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TF1.h>
#include <TArrow.h>
#include <TPaveStats.h>
Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==7 || argc==8) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_diagonal (char*)stream (char*)setname (double)used_nstream (double)used_nset (char*)tag (char*)brname [(int)fl_appRun]" << std::endl
					<< "[  stream  ] 0,1,2,3,4,5,01,0-5" << std::endl
					<< "[  setname ] A,B,..,U, AB, A-U"  << std::endl
					<< std::endl, abort();
  Char_t*  stream       = argv[1];
  Char_t*  setname      = argv[2];
  Double_t used_nstream = atof(argv[3]);
  Double_t used_nset    = atof(argv[4]);
  Char_t*  tag          = argv[5];
  Char_t*  brname       = argv[6];
  Int_t    fl_appRun    = 1;
  if( argc==8 ) fl_appRun = atoi( argv[7] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Int_t fl_lrnb = 0; // 0(LR), 1(NB)
  if(      strncmp(brname, "lr", 2) == 0 ) fl_lrnb=0; // LR
  else if( strncmp(brname, "nb", 2) == 0 ) fl_lrnb=1; // NB
  else std::cerr << "[ABORT] Wrong branch-name : " << brname << std::endl, abort();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t    Nchain             = 16; // fixed !! [4 category]x[2 lepton-type]x[2 Xs-region], [ sigMC, sigMC(t,t), sigMC(except-t,t), gMC ]-> [low-xs(ee)] : [high-xs(ee)] : [low-xs(mm)] : [high-xs(mm)]
  const Int_t    Nhist              =  4; // fixed !! [ sig/bkg  ]                x[2 Xs-region], [ sig, bkg ] -> [low-xs] : [high-xs]
  const Int_t    nfile[Nchain]      = {0};
  const Int_t    fl_mode_ll[Nchain] = { 1,1,1,1,1,1,1,1,   // 1(ee)
					0,0,0,0,0,0,0,0 }; // 0(mm)
  const Char_t*  indir_sig[2]       = { Form("hbk_bcs/sig/sig%s_522_merge/", tag), Form("NB/hbk/hbk_%s_522_bcs_merge/",tag) };
  const Char_t*  indir_bkg[2]       = { Form("hbk_bcs/bkg/bkg%s_522/",       tag), Form("NB/hbk/hbk_%s_522_bcs/",      tag) };
  const Char_t*  tail               = "*.root";
  const Double_t scale_event_sig    = sigmc_amount*(used_nset/(Double_t)nset); // sigmc : N -> N/alpha
  const Double_t scale_event_bkg    = used_nstream;                            //   gmc : N -> N/alpha 
  const Bool_t   fl_message         = !true;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if(      i%4 < 2.5 ) sTmp << indir_sig[fl_lrnb] << "/sigMC_*_set[" << setname << "]"; // sig
    else if( i%4 > 2.5 ) sTmp << indir_bkg[fl_lrnb] << "/gMC_*_s0["    << stream  << "]"; // bkg
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0}, // sig(low-xs)
    {0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0}, // bkg(low-xs)
    {0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0}, // sig(high-xs)
    {0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1}, // bkg(high-xs)
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( int i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[1024];

  for( int i=0; i<Nchain; i++ ) {
    if(      i%4==0 ) sTmp << "self==1";                              // true
    else if( i%4==1 ) sTmp << "self!=1 && " << makeCut_q2fl( 1,1   ); // false  with  q2=true, fl=true
    else if( i%4==2 ) sTmp << "self!=1 && " << makeCut_q2fl( 1,1,1 ); // false except q2=true, fl=true
    else if( i%4==3 ) sTmp << "1";
    
    if(      (i/4)%2==0 ) sTmp << " && xs_m < 1.1"; // low-Xs
    else if( (i/4)%2==1 ) sTmp << " && xs_m > 1.1"; // high-Xs
    
    strcpy( add_cut[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  const Char_t*  tname     = "h511";
  const Char_t*  axis      = brname;
  const Double_t offset    =  0.0;
  const Int_t    xbin  [2] = {100,     200 }; // {LR, NB}
  const Double_t xmin  [2] = { 0.0,   -1.0 }; // {LR, NB}
  const Double_t xmax  [2] = { 1.0,    1.0 }; // {LR, NB}
  const Char_t*  xlabel[2] = { "LR",  "NB" }; // {LR, NB}
  
  const Bool_t flag_save  = true; // outfile.eps and outfile.root
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",4 );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
    //chain[j]->GetCut()->Set( "de",  1, -0.50, 0.0, 0.50 );
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
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin[fl_lrnb],offset+xmin[fl_lrnb],offset+xmax[fl_lrnb] );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
    //////////////////////////////// correction events with LR=1 or NB=1
    if( chain[j]->GetTree()->GetEntries( Form("(%s)&&(%s>1)", add_cut[j],axis) ) ){
      std::cerr << "[ABORT] Strange Evnets having LR>1 or NB>1  Exist !" << std::endl;
      std::cerr << Form( "(%s)&&(%s>1)", add_cut[j],axis ) << std::endl;
      chain[j]->GetTree()->SetScanField(0);
      chain[j]->GetTree()->Scan( "*", Form("(%s)&&(%s>1)", add_cut[j],axis ) ); 
      abort();
    }
    Double_t tmp_entry = tmphist[j]->GetEntries();
    Double_t tmp_bin   = tmphist[j]->GetBinContent( xbin[fl_lrnb]   );
    Double_t tmp_over  = tmphist[j]->GetBinContent( xbin[fl_lrnb]+1 );
    tmphist[j]->SetBinContent( xbin[fl_lrnb],   tmp_bin+tmp_over );
    tmphist[j]->SetBinContent( xbin[fl_lrnb]+1, 0                );
    tmphist[j]->SetEntries( tmp_entry );
    ////////////////////////////////
    tmphist[j]->Sumw2();
    if(      j%4 < 2.5 ) tmphist[j]->Scale( 1/scale_event_sig );
    else if( j%4 > 2.5 ) tmphist[j]->Scale( 1/scale_event_bkg );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    if(      i==0 ) hist[i] = new TH1D( Form("hist_%s_lowxs_sig",  axis), Form("hist_%s_lowxs_sig",  axis), xbin[fl_lrnb],offset+xmin[fl_lrnb],offset+xmax[fl_lrnb] );
    else if( i==1 ) hist[i] = new TH1D( Form("hist_%s_lowxs_bkg",  axis), Form("hist_%s_lowxs_bkg",  axis), xbin[fl_lrnb],offset+xmin[fl_lrnb],offset+xmax[fl_lrnb] );
    else if( i==2 ) hist[i] = new TH1D( Form("hist_%s_highxs_sig", axis), Form("hist_%s_highxs_sig", axis), xbin[fl_lrnb],offset+xmin[fl_lrnb],offset+xmax[fl_lrnb] );
    else if( i==3 ) hist[i] = new TH1D( Form("hist_%s_highxs_bkg", axis), Form("hist_%s_highxs_bkg", axis), xbin[fl_lrnb],offset+xmin[fl_lrnb],offset+xmax[fl_lrnb] );
    Deco( hist[i], 2, i+1, i+1 );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }
  // +++++++ calculation +++++++++++++++++++++
  TGraph**       graph_eff  = new TGraph*      [Nhist/2]; // low-xs, high-xs
  TGraphErrors** graph_sign = new TGraphErrors*[Nhist/2]; // low-xs, high-xs
  for( Int_t i=0; i<Nhist/2; i++ ){
    graph_eff[i]  = new TGraph();
    graph_sign[i] = new TGraphErrors();
    Deco( graph_eff[i],  0, i+1, i+1 );
    Deco( graph_sign[i], 0, i+1, i+1 );
  }
  graph_eff [0]->SetName( "eff_lowxs"   );
  graph_eff [1]->SetName( "eff_highxs"  );
  graph_sign[0]->SetName( "sign_lowxs"  );
  graph_sign[1]->SetName( "sign_highxs" );
    
  Double_t entry_sig[Nhist/2] = { hist[0]->Integral(), hist[2]->Integral() }; // { low-xs, high-xs }
  Double_t entry_bkg[Nhist/2] = { hist[1]->Integral(), hist[3]->Integral() }; // { low-xs, high-xs }

  std::cout << "=================================" << std::endl
	    << "[Scale Factor]"                    << std::endl
    	    << "sig : " << scale_event_sig
	    << "("      << used_nset << " set)"    << std::endl
	    << "bkg : " << scale_event_bkg         << std::endl
	    << "=================================" << std::endl;
  if( fl_message ){
    std::cout << "===================================================" << std::endl
	      << "nsig( low-xs) = " << entry_sig[0] << std::endl
	      << "nsig(high-xs) = " << entry_sig[1] << std::endl
	      << "nbkg( low-xs) = " << entry_bkg[0] << std::endl
	      << "nbkg(high-xs) = " << entry_bkg[1] << std::endl
	      << "===================================================" << std::endl;
  }

  Double_t max_sign[Nhist/2] = {0}; // maximum significance             [low-xs, high-xs]
  Double_t max_sig [Nhist/2] = {0}; // # of sig at maximum significance [low-xs, high-xs]
  Double_t max_bkg [Nhist/2] = {0}; // # of bkg at maximum significance [low-xs, high-xs]
  Double_t max_lr  [Nhist/2] = {0}; // LR value at maximum significance [low-xs, high-xs]
  for( Int_t k=0; k<xbin[fl_lrnb]; k++ ){
    for( Int_t i=0; i<Nhist/2; i++ ){
      Double_t nsig_err   = 0;
      Double_t nbkg_err   = 0;
      Double_t nsig       = hist[2*i  ]->IntegralAndError( k+1, xbin[fl_lrnb], nsig_err );
      Double_t nbkg       = hist[2*i+1]->IntegralAndError( k+1, xbin[fl_lrnb], nbkg_err );
      Double_t sign       = nsig==0 ? 0 : nsig / sqrt(nsig+nbkg);
      Double_t round_sig  = nsig==0 ? 0 : sign/nsig - sign*sign*sign/nsig/nsig/2.0;
      Double_t round_bkg  = nsig==0 ? 0 : -sign*sign*sign/nsig/nsig/2.0;
      Double_t sign_error = sqrt( round_sig*round_sig*nsig_err*nsig_err + round_bkg*round_bkg*nbkg_err*nbkg_err );
      graph_eff [i]->SetPoint     ( k, nsig/entry_sig[i],           1.0-nbkg/entry_bkg[i]     );
      graph_sign[i]->SetPoint     ( k, hist[i]->GetBinLowEdge(k+1), sign                      );
      graph_sign[i]->SetPointError( k,                           0, sign_error                );
      if( sign > max_sign[i] ){
	max_sign[i] = sign;
	max_sig [i] = nsig;
	max_bkg [i] = nbkg;
	max_lr  [i] = hist[2*i]->GetBinLowEdge(k+1);
      }
      if( fl_message ){
	Double_t tmp_x, tmp_y;
	std::cout << "Bin = " << k+1 << ", i = " << i << std::endl;
	std::cout << " S = " << std::setw(7) << std::right << nsig << ", B = "   << std::setw(7) << std::right << nbkg << std::endl;
	graph_eff [i]->GetPoint( k, tmp_x, tmp_y ); std::cout << "[ eff(sig) = " << std::setw(8) << std::right << tmp_x << ", eff(bkg) = " << std::setw(8) << std::right << tmp_y << " ] ";
	graph_sign[i]->GetPoint( k, tmp_x, tmp_y ); std::cout << "[ LH = "       << std::setw(8) << std::right << tmp_x << ", sign = "     << std::setw(8) << std::right << tmp_y << " ] " << std::endl;
      }
    }
  }

  // +++++++ display ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin[fl_lrnb]+1);
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
    hist[i]->Scale( 1/hist[i]->Integral() );
  }

  std::cout << "[ Maximum Significane ]"  << std::endl;
  for( Int_t i=0; i<Nhist/2; i++ ){
    if(      i==0 ) std::cout << " ( low-xs) ";
    else if( i==1 ) std::cout << " (high-xs) ";
    std::cout << std::setw(10) << std::right << max_sign[i] << " ( S = "
	      << std::setw( 7) << std::right << max_sig [i] << ", B = "
	      << std::setw( 7) << std::right << max_bkg [i] << ", " << xlabel[fl_lrnb] << " = "
	      << std::setw( 5) << std::right << max_lr  [i] << " )"
	      << std::endl;
  }

  std::cout << " ( total ) "
	    << std::setw(10) << std::right << (max_sig[0] + max_sig[1])/sqrt(max_sig[0]+max_sig[1]+max_bkg[0]+max_bkg[1]) << " ( S = "
	    << std::setw( 7) << std::right <<  max_sig[0] + max_sig[1] << ", B = "
	    << std::setw( 7) << std::right <<  max_bkg[0] + max_bkg[1] << " )"
	    << std::endl;
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);
  TH2D* waku = Waku( Nhist, hist, xlabel[fl_lrnb] );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s", xlabel[fl_lrnb]) );
  waku->Draw();

  for( Int_t i=Nhist-1; i>=0; i-- ) hist[i]->Draw( "same" );
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[0],"sig(Low-X_{s})",  "L" );
  legend1->AddEntry( hist[1],"bkg(Low-X_{s})",  "L" );
  legend1->AddEntry( hist[2],"sig(High-X_{s})", "L" );
  legend1->AddEntry( hist[3],"bkg(High-X_{s})", "L" );
  legend1->Draw();
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->cd(2);
  TH2D* waku_eff  = new TH2D( "Performance of Likelihood Ratio Cut", "Performance of Likelihood Ratio Cut", 2,0,1,2,0,1  );
  waku_eff->GetXaxis()->CenterTitle();
  waku_eff->GetYaxis()->CenterTitle();
  waku_eff->SetXTitle( "Sig. Eff."           );
  waku_eff->SetYTitle( "Bkg. Rejection Eff." );
  waku_eff->Draw();
  for( Int_t i=0; i<Nhist/2; i++ ) graph_eff[i]->Draw("Psame");

  // +++++++ tlegend2 ++++++++++++++++++++++++++++++++++
  TLegend* legend2 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend2->AddEntry( graph_eff[0],"Low-X_{s}",  "L" );
  legend2->AddEntry( graph_eff[1],"High-X_{s}", "L" );
  legend2->Draw();

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->cd(3);
  TH2D* waku_sign = new TH2D( "Significance", "Significance", 2,xmin[fl_lrnb],xmax[fl_lrnb],2,0,GetYMax(Nhist/2, graph_sign) );
  waku_sign->GetXaxis()->CenterTitle();
  waku_sign->GetYaxis()->CenterTitle();
  waku_sign->SetXTitle( xlabel[fl_lrnb] );
  waku_sign->SetYTitle( "Significance"  );
  waku_sign->Draw();
  for( Int_t i=0; i<Nhist/2; i++ ) graph_sign[i]->Draw("Psame");
  // +++++++ tlegend3 ++++++++++++++++++++++++++++++++++
  TLegend* legend3 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend3->AddEntry( graph_sign[0],"Low-X_{s}",  "L" );
  legend3->AddEntry( graph_sign[1],"High-X_{s}", "L" );
  legend3->Draw();
  // +++++++ tlegend ++++++++++++++++++++++++++++++++++
  c1->cd(4);
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
    c1->Print(     Form("pic/1d_%s_s0%s_set%s.eps",  axis,stream,setname)             );
    TFile outfile( Form("pic/1d_%s_s0%s_set%s.root", axis,stream,setname), "RECREATE" );
    c1->Write();
    for( Int_t i=0; i<Nhist;   i++ ) hist[i]->Write();
    for( Int_t i=0; i<Nhist/2; i++ ) graph_eff[i]->Write(), graph_sign[i]->Write();
    outfile.Close();
  }

  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] add_cut;
  delete legend;
  delete legend1;
  delete legend2;
  delete legend3;
  delete box;
  delete waku;
  delete waku_eff;
  delete waku_sign;
  delete c1;
  return 0;
}

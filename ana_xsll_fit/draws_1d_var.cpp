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

using namespace sig_gmc_rd_comp;
using namespace bgsup_var;

const Bool_t flag_k4pi     = !true; // 1(veto  K4pi    modes)
const Bool_t flag_unflavor = !true; // 1(veto unflavor modes)
const Bool_t flag_ccpi0    = !true; // 1(veto ccpi0 peak    )
const Bool_t flag_save     = true; // outfile.eps and outfile.root
const Bool_t flag_norm     = true;

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==6 || argc==7) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_var (int)fl_axis (char*)stream (char*)setname (double)used_nstream (double)used_nset[(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Int_t    fl_axis      = atoi(argv[1]);
  Char_t*  stream       = argv[2];
  Char_t*  setname      = argv[3];
  Double_t used_nstream = atof(argv[4]); 
  Double_t used_nset    = atof(argv[5]);
  Int_t   fl_appRun      = 1;
  if( argc==7 ) fl_appRun = atoi( argv[6] );

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain             = 6*2; // [gmc,gmc(emu),gmc(jpsi),sigmc,rd(emu),rd(jpsi)] x [ee,mm]
  const Int_t Ncategory          = 7;   // [gmc(qq,bb),gmc(emu),gmc(jpsi),sigmc,rd(emu),rd(jpsi)]
  const Int_t Nplot              = 8;   // [bkg(qq,bb,qq+bb),bkg(emu),bkg(jpsi),sigmc,rd(emu),rd(jpsi)]
  const Int_t Ntmp               = Ncategory*2; // x[ee,mm      ]
  const Int_t Nhist              = Nplot    *3; // x[ee,mm,ee+mm]
  const Int_t fl_sb[Ntmp]        = {0,0,1,2,3, 4, 5,  // 0[bkg], 1[bkg(emu)], 2[bkg(jpsi)], 3[sig],  4[rd(emu)],  5[rd(jpsi)] for ee
				    6,6,7,8,9,10,11}; // 6[bkg], 7[bkg(emu)], 8[bkg(jpsi)], 9[sig], 10[rd(emu)], 11[rd(jpsi)] for mm
  const Int_t fl_mode_ll[Nchain] = {1,2,1,1,2,1,  // 1(ee),2(emu)
				    0,2,0,0,2,0}; // 0(mm)

  const Double_t scale_event_sig = sigmc_amount*(used_nset/(Double_t)nset); // sigmc : N -> N/alpha
  const Double_t scale_event_bkg = used_nstream;                            // gmc   : N -> N/alpha
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t* LRNB_cut = new Char_t[4096];
  LRNB_cut = "1";
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if     ( i==0 || i== 6 ) sTmp << indir[0] << "gMC_*s0["    << stream  << "]_"; // bkg
    else if( i==1 || i== 7 ) sTmp << indir[1] << "gMC_*s0["    << stream  << "]_"; // bkg(emu)
    else if( i==2 || i== 8 ) sTmp << indir[2] << "gMC_*s0["    << stream  << "]_"; // bkg(jpsi)
    //else if( i==2 || i== 8 ) sTmp << indir[2] << "gMC_*e031*s0["    << stream  << "]_"; // bkg(jpsi) // small sample
    else if( i==3 || i== 9 ) sTmp << indir[3] << "sigMC_*set[" << setname << "]_"; // sig
    else if( i==4 || i==10 ) sTmp << indir[4] << "RD_";                            // rd(emu)
    else if( i==5 || i==11 ) sTmp << indir[5] << "RD_";                            // rd(jpsi)
    //else if( i==5 || i==11 ) sTmp << indir[5] << "RD_*e031";                            // rd(jpsi) // small sample
    else                   std::cerr << "[ABORT] Wrong fl_sb : " << fl_sb[i] << std::endl, abort();
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Ntmp] ={
    {1},             // gmc(qq),    ee
    {0,1},           // gmc(bb),    ee
    {1,1},           // gmc(qq+bb), ee
    {0,0,1},         // gmc(emu),   ee
    {0,0,0,1},       // gmc(jpsi),  ee
    {0,0,0,0,1},     // sigmc,      ee
    {0,0,0,0,0,1},   // rd(emu),    ee
    {0,0,0,0,0,0,1}, // rd(jpsi),   ee
    {0,0,0,0,0,0,0,1},             // gmc(qq),    mm
    {0,0,0,0,0,0,0,0,1},           // gmc(bb),    mm
    {0,0,0,0,0,0,0,1,1},           // gmc(qq+bb), mm
    {0,0,0,0,0,0,0,0,0,1},         // gmc(emu),   mm
    {0,0,0,0,0,0,0,0,0,0,1},       // gmc(jpsi),  mm
    {0,0,0,0,0,0,0,0,0,0,0,1},     // sigmc,      mm
    {0,0,0,0,0,0,0,0,0,0,0,0,1},   // rd(emu),    mm
    {0,0,0,0,0,0,0,0,0,0,0,0,0,1}, // rd(jpsi),   mm
    {1,0,0,0,0,0,0,1},             // gmc(qq),    ee+mm
    {0,1,0,0,0,0,0,0,1},           // gmc(bb),    ee+mm
    {1,1,0,0,0,0,0,1,1},           // gmc(qq+bb), ee+mm
    {0,0,1,0,0,0,0,0,0,1},         // gmc(emu),   ee+mm
    {0,0,0,1,0,0,0,0,0,0,1},       // gmc(jpsi),  ee+mm
    {0,0,0,0,1,0,0,0,0,0,0,1},     // sigmc,      ee+mm
    {0,0,0,0,0,1,0,0,0,0,0,0,1},   // rd(emu),    ee+mm
    {0,0,0,0,0,0,1,0,0,0,0,0,0,1}, // rd(jpsi),   ee+mm
  };

  Char_t** add_cut = new Char_t*[Ntmp];
  for( Int_t i=0; i<Ntmp; i++ ){
    add_cut[i] = new Char_t[4096];
    sTmp << LRNB_cut;
    if     ( i%Ncategory==4 ) sTmp << "&& self==1";   // sigmc
    else if( i%Ncategory==0 ) sTmp << "&& genbfl==0"; // gmc(qq)
    else if( i%Ncategory==1 ) sTmp << "&& genbfl!=0"; // gmc(bb)
    /*
    sTmp << "&& ( " // K or K*
	 << Form( " (rm_xs==  1) || " ) // K+
	 << Form( " (rm_xs== 10) || " ) // Ks
	 << Form( " (rm_xs==101  && abs(xs_m-%f)<0.050) || ", PDGmass::kstr0 ) // K+pi-
	 << Form( " (rm_xs==110  && abs(xs_m-%f)<0.050) || ", PDGmass::kstrp ) // Kspi-
	 << Form( " (rm_xs==1001 && abs(xs_m-%f)<0.050) || ", PDGmass::kstrp ) // K+pi0
	 << Form( " (rm_xs==1010 && abs(xs_m-%f)<0.050) || ", PDGmass::kstr0 ) // Kspi0
	 << "0 )";
    */
    
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    if( flag_ccpi0    ) sTmp << " && " << cut_ccpi0;
    strcpy( add_cut[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Ntmp];
  TH1D**    hist    = new TH1D*  [Nhist];
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis[fl_axis]) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), 0, tail );
    nominal_cut_selection( chain[j], fl_mode_ll[j] )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.26 );
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.29 );
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
  }
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make tmphist *************************************" << std::endl;
  for( Int_t j=0; j<Ntmp; j++ ){
    std::cout << Form( "add_cut%d(chain-tree%d) : ", j, fl_sb[j])  << add_cut[j] << std::endl;
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[fl_sb[j]]->GetChange()), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
    if     ( fl_axis==2 && strcmp(axis[fl_axis],"kfbchi/kfbdgf")==0 && j%Ncategory==3 ) chain[fl_sb[j]]->GetTree()->Project( Form("tmphist%d",j), "kfbchicalib/kfbdgf", add_cut[j] ), std::cout << Form("calib %d", j) << std::endl;
    else if( fl_axis==3 && strcmp(axis[fl_axis],"kfbcl"        )==0 && j%Ncategory==3 ) chain[fl_sb[j]]->GetTree()->Project( Form("tmphist%d",j), "kfbclcalib",         add_cut[j] ), std::cout << Form("calib %d", j) << std::endl;
    else if( fl_axis==4 && strcmp(axis[fl_axis],"10000*dzll3d" )==0 && j%Ncategory==3 ) chain[fl_sb[j]]->GetTree()->Project( Form("tmphist%d",j), "10000*dzll3dcalib",  add_cut[j] ), std::cout << Form("calib %d", j) << std::endl;
    else if( fl_axis==6 && strcmp(axis[fl_axis],"de"           )==0 && j%Ncategory==3 ) chain[fl_sb[j]]->GetTree()->Project( Form("tmphist%d",j), "decalib",            add_cut[j] ), std::cout << Form("calib %d", j) << std::endl;
    else                                                                                chain[fl_sb[j]]->GetTree()->Project( Form("tmphist%d",j), axis[fl_axis], add_cut[j] );
    
    tmphist[j]->Sumw2();
  }
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
    
    for( Int_t j=0; j<Ntmp; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }

  if( flag_norm ){
    for( Int_t i=0; i<Nhist; i++ ){
      if( fl_axis==6 ){
	if     ( i== 0 || i== 1 ) hist[i]->Scale( 1/hist[ 2]->Integral() );
	else if( i== 8 || i== 9 ) hist[i]->Scale( 1/hist[10]->Integral() );
	else if( i==16 || i==17 ) hist[i]->Scale( 1/hist[18]->Integral() );
	else                      hist[i]->Scale( 1/hist[ i]->Integral() );
      }else{
	if     ( i== 0 || i== 1 ) hist[i]->Scale( 1/hist[ 2]->GetEntries() );
	else if( i== 8 || i== 9 ) hist[i]->Scale( 1/hist[10]->GetEntries() );
	else if( i==16 || i==17 ) hist[i]->Scale( 1/hist[18]->GetEntries() );
	else                      hist[i]->Scale( 1/hist[ i]->GetEntries() );
      }
    }
  }
  
  for( Int_t i=0; i<3; i++ ){
    Deco( hist[Nplot*i+0], 1, 8, 8 ); // gmc(qq)
    Deco( hist[Nplot*i+1], 1, 4, 4 ); // gmc(bb)
    Deco( hist[Nplot*i+2], 1, 3, 3 ); // gmc(qq+bb)
    Deco( hist[Nplot*i+3], 1, 3, 3 ); // gmc(emu)
    Deco( hist[Nplot*i+4], 1, 2, 2 ); // gmc(jpsi)
    Deco( hist[Nplot*i+5], 1, 2, 2 ); // sigmc
    Deco( hist[Nplot*i+6], 1, 3, 3 ); // rd(emu)
    Deco( hist[Nplot*i+7], 1, 2, 2 ); // rd(jpsi)
    hist[Nplot*i+0]->SetLineStyle(2); // gmc(qq)
    hist[Nplot*i+1]->SetLineStyle(2); // gmc(bb)
  }
  
  // +++++++ display ++++++++++++++++++++++++++++++++++
  std::cout << "=================================" << std::endl
	    << "[Scale Factor]"                    << std::endl
	    << "sig : " << scale_event_sig         << std::endl
	    << "bkg : " << scale_event_bkg         << std::endl
	    << "=================================" << std::endl;
  Double_t entry_sig[Nhist] = {0}; // # of events in signal box region
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin[fl_axis]+1);
    for( Int_t j=hist[i]->FindBin(5.27+0.000001); j<=xbin[fl_axis]; j++ ) entry_sig[i] += hist[i]->GetBinContent(j);

    for( Int_t j=0; j<Ntmp; j++ ) if( add[i][j] ) std::cout << j << ",";
    std::cout <<  ")[ "
	      << std::setw(12) << std::right << entry_all
	      <<  " events ( canvas : "
	      << std::setw(12) << std::right << std::setprecision(3) << entry_canvas
	      << " / under : "
	      << std::setw(12) << std::right << std::setprecision(3) << entry_under
	      << " / over  : "
	      << std::setw(12) << std::right << std::setprecision(3) << entry_over
	      << " / sig  : "
	      << std::setw(12) << std::right << std::setprecision(3) << entry_sig[i]
	      << "]"
	      << std::endl;

  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  TCanvas* c1 = Canvas( "c1","c1",3, 4 );
  c1->Draw();

  TH2D** waku = new TH2D*[3];

  if( fl_axis==3 ){ // for kfbcl
    waku[0] = Waku( Nhist, hist, xlabel[fl_axis], "ee",    "ee",        2 );
    waku[1] = Waku( Nhist, hist, xlabel[fl_axis], "mm",    "#mu#mu",    2 );
    waku[2] = Waku( Nhist, hist, xlabel[fl_axis], "ee;mm", "ee+#mu#mu", 2 );
  }else{
    waku[0] = Waku( Nhist, hist, xlabel[fl_axis], "ee",    "ee"        );
    waku[1] = Waku( Nhist, hist, xlabel[fl_axis], "mm",    "#mu#mu"    );
    waku[2] = Waku( Nhist, hist, xlabel[fl_axis], "ee+mm", "ee+#mu#mu" );
  }

  TPaveText** box = new TPaveText*[2];
  box[0] = new TPaveText( 0.05, 0.87, 0.70, 0.93, "BRNDC" ); box[0]->AddText("gmc(emu,jpsi[dot])-gmc(qq+bb,sig[hist])");
  box[1] = new TPaveText( 0.05, 0.87, 0.70, 0.93, "BRNDC" ); box[1]->AddText( "rd(emu,jpsi[dot])-gmc(emu,jpsi[hist])" );
  
  for(Int_t i=0; i<3; i++ ){

    c1->cd(i+1);
    if( fl_axis==3 ) gPad->SetLogy(); // for CL(B-vtx)
    waku[i]->Draw();
    box[0]->Draw();
    hist[Nplot*i+2]->DrawCopy("hist same"); // gmc(qq+bb)
    hist[Nplot*i+3]->DrawCopy("same");      // gmc(emu)
    hist[Nplot*i+4]->DrawCopy("same");      // gmc(jpsi)
    hist[Nplot*i+5]->DrawCopy("hist same"); // sigmc

    c1->cd(i+4);
    if( fl_axis==3 ) gPad->SetLogy(); // for CL(B-vtx)
    waku[i]->Draw();
    box[1]->Draw();
    hist[Nplot*i+3]->DrawCopy("hist same"); // gmc(emu)
    hist[Nplot*i+4]->DrawCopy("hist same"); // gmc(jpsi)
    hist[Nplot*i+6]->DrawCopy("same"); // rd(emu)
    hist[Nplot*i+7]->DrawCopy("same"); // rd(jpsi)

    c1->cd(i+7);
    if( fl_axis==3 ) gPad->SetLogy(); // for CL(B-vtx)
    waku[i]->Draw();
    hist[Nplot*i+0]->DrawCopy("hist same"); // gmc(qq)
    hist[Nplot*i+1]->DrawCopy("hist same"); // gmc(bb)
    hist[Nplot*i+2]->DrawCopy("hist same"); // gmc(qq+bb)
    hist[Nplot*i+5]->DrawCopy("hist same"); // sigmc

    c1->cd(i+10);
    if( fl_axis==3 ) gPad->SetLogy(); // for CL(B-vtx)
    waku[i]->Draw();
    hist[Nplot*i+0]->SetLineStyle(1);
    hist[Nplot*i+1]->SetLineStyle(1);
    hist[Nplot*i+0]->DrawNormalized("hist same"); // gmc(qq)
    hist[Nplot*i+1]->DrawNormalized("hist same"); // gmc(bb)
    hist[Nplot*i+5]->DrawCopy("hist same"); // sigmc
  }

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[0], "gmc(qq)",    "L" );
  legend1->AddEntry( hist[1], "gmc(bb)",    "L" );
  legend1->AddEntry( hist[2], "gmc(qq+bb)", "L" );
  legend1->AddEntry( hist[3], "gmc(emu)",   "L" );
  legend1->AddEntry( hist[4], "gmc(jpsi)",  "L" );
  legend1->AddEntry( hist[5], "sigmc",      "L" );
  legend1->AddEntry( hist[6], "rd(emu)",    "L" );
  legend1->AddEntry( hist[7], "rd(jpsi)",   "L" );
  legend1->Draw();

  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_set%s_s0%s.eps",  fname[fl_axis], setname, stream) );
    c1->Print( Form("pic/%s_set%s_s0%s.root", fname[fl_axis], setname, stream) );
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] add_cut;
  delete   c1;
    
  return 0;
}


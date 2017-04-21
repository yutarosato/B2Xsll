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
#include <TFile.h>


//using namespace sig_gmc_rd_cut2;
using namespace sig_gmc_rd_cut2_all;
//using namespace Mbc_comb;
using namespace Mbc_bkg;

const Bool_t flag_norm     = true;
const Bool_t flag_k4pi     = true; // 1(veto  K4pi    modes)
const Bool_t flag_unflavor = true; // 1(veto unflavor modes)
const Bool_t flag_ccpi0    = true; // 1(veto ccpi0 peak    )
const Bool_t flag_save     = true; // outfile.eps and outfile.root

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==2 || argc==3) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_Mbc (char*) setname [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Char_t*  setname      = argv[1];
  Int_t    fl_appRun    = 1;
  if( argc==3 ) fl_appRun = atoi( argv[2] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t fl_5body_veto    = !true;
  const Int_t Nchain           = 2; // [ee(sigmc),mm(sigmc)]
  const Int_t Nctgry           = 5; // [true, false with (t,t), false with (t,f), false with(f,t), false with(f,f)]
  const Int_t Nplot            = 4; // [true, false with (t,t), false with (t,f), other false]
  const Int_t Ntmp             = 2*Nctgry; // [ee,mm      ]
  const Int_t Nhist            = 3*Nplot;  // [ee,mm,ee+mm]
  const Int_t nfile[Nchain]    = {0};
  const Int_t fl_sb[Nchain]    = {1,1}; // 0(bkg), 1(sig)
  const Int_t fl_mode_ll[Ntmp] = {1,1,1,1,1,
				  0,0,0,0,0}; // 1(e), 0(mu)
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t* LRNB_cut = new Char_t[4096];
  LRNB_cut = "1";
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<3; i++ ){
    infile[i] = new Char_t[1024];
    //if     ( fl_sb[i]==0 ) sTmp << indir[fl_sb[i]] << "gMC_*_s0[" << stream << "]";     // bkg
    //else if( fl_sb[i]==1 ) sTmp << indir[fl_sb[i]] << "sigMC_*_set[" << setname << "]"; // sig
    //else if( fl_sb[i]==2 ) sTmp << indir[fl_sb[i]] << "RD_";                            // rd
    if( fl_sb[i]==1 ) sTmp << indir[fl_sb[i]] << "sigMC_*_set[" << setname << "]"; // sig
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Ntmp] ={
    {1},         // true
    {0,1},       // (t,t)
    {0,0,1},     // (t,f)
    {0,0,0,1,1}, // other false
    {0,0,0,0,0,1},         // true
    {0,0,0,0,0,0,1},       // (t,t)
    {0,0,0,0,0,0,0,1},     // (t,f)
    {0,0,0,0,0,0,0,0,1,1}, // other
    {1,0,0,0,0,1},           // true
    {0,1,0,0,0,0,1},         // (t,t)
    {0,0,1,0,0,0,0,1},       // (t,f)
    {0,0,0,1,1,0,0,0,1,1},   // other
  };
  
  Char_t** add_cut = new Char_t*[Ntmp];
  for( Int_t i=0; i<Ntmp; i++ ) add_cut[i] = new Char_t[4096];
  for( Int_t i=0; i<Ntmp; i++ ){
    if     ( i%Nctgry==0 ) sTmp << LRNB_cut << " && self==1";   // true
    else if( i%Nctgry==1 ) sTmp << LRNB_cut << " && self!=1 &&" // false with (t,t)
				<< makeCut_q2fl( 1, 1 );
    else if( i%Nctgry==2 ) sTmp << LRNB_cut << " && self!=1 &&" // false with (t,f)
				<< makeCut_q2fl( 1, 0 );
    else if( i%Nctgry==3 ) sTmp << LRNB_cut << " && self!=1 &&" // false with (f,t)
				<< makeCut_q2fl( 0, 1 );
    else if( i%Nctgry==4 ) sTmp << LRNB_cut << " && self!=1 &&" // false with (f,f)
				<< makeCut_q2fl( 0, 0 );
    
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
  TCanvas*  c1      = Canvas( "c1","c1", 2, 2 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[fl_sb[j]], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], j )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 );
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
	    << " ************************ make hist(tmp) *************************************" << std::endl;
  for( Int_t j=0; j<Ntmp; j++ ){
    std::cout << Form( "<tmphist %d > add_cut : ", j ) << add_cut[j] << std::endl;
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[fl_mode_ll[j]]->GetChange()), xbin,offset+xmin,offset+xmax );
    chain[fl_mode_ll[j]]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );

  }
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist(plot) *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    Deco( hist[i], 1, i/Nplot+2, i/Nplot+2 );
    for( Int_t j=0; j<Ntmp; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
    hist[i]->Sumw2();
    if( flag_norm ) hist[i]->Scale( 1/hist[i]->Integral() );
    if     ( i%Nplot==0 ) hist[i]->SetTitle( "true"          );
    else if( i%Nplot==1 ) hist[i]->SetTitle( "(q2,fl)=(t,t)" );
    else if( i%Nplot==2 ) hist[i]->SetTitle( "(q2,fl)=(t,f)" );
    else if( i%Nplot==3 ) hist[i]->SetTitle( "(q2,fl)=(f,?)" );
  }


  // +++++++ display ++++++++++++++++++++++++++++++++++
  Double_t entry_sig[Nhist] = {0}; // # of events in signal box region
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
    for( Int_t j=hist[i]->FindBin(5.27+0.000001); j<=xbin; j++ ) entry_sig[i] += hist[i]->GetBinContent(j);
    TString sTmp = "Added-Files( ";
    for( Int_t j=0; j<Ntmp; j++ ) if( add[i][j] ) std::cout << j <<  ",";
    std::cout <<  ")["
	      << entry_all
	      << " events ( canvas : "
	      << entry_canvas
	      << " / under : "
	      << entry_under
	      << " / over  : "
	      << entry_over
	      << " / sig  : "
	      << entry_sig[i]
	      << "]"
	      << std::endl;
  }

  std::cout << std::endl;

  Double_t entry_sig_each[Ntmp] = {0}; // # of events in signal box region (tmphist)
  for( Int_t i=0; i<Ntmp; i++ ){
    std::cout << Form("<tmphist %d> ",i);
    for( Int_t j=tmphist[i]->FindBin(5.27+0.000000001); j<=xbin; j++ ) entry_sig_each[i] += tmphist[i]->GetBinContent(j);
    std::cout << tmphist[i]->Integral() << " events(canvas), "
	      << entry_sig_each[i]      << " events(signal-box)" << std::endl;
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  for(Int_t i=0; i<Nhist; i++ ){
    c1->cd(i%Nplot+1);
    if     ( i/Nplot==0 ) hist[i]->Draw("P"); // ee
    else if( i/Nplot==1 ) hist[i]->Draw("Psame"); // mm
    //else if( i/Nplot==2 ) hist[i]->Draw("Psame"); // ee+mm
  }
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[0+0*Nplot],"ee",        "P" );
  legend1->AddEntry( hist[0+1*Nplot],"#mu#mu",    "P" );
  legend1->AddEntry( hist[0+2*Nplot],"ee+#mu#mu", "P" );
  legend1->Draw();

  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_self_cf_tot_set%s.eps",  axis, setname) );
    TFile outfile( Form("pic/%s_self_cf_tot_set%s.root", axis, setname), "RECREATE" );
    c1->Write();
    for( Int_t i=0; i<Nhist; i++ ) hist[i]->Write();
    outfile.Close();
  }


  // +++++++++++++++ Log for Roofit +++++++++++++++++++++
  std::cout << " +++++++++++++++ Log for Roofit +++++++++++++++++++++ "  << std::endl
	    << "  r_tt[0] = new RooRealVar( \"r_tt0\", \"r_tt0\", " << std::setw(10) << std::right << hist[ 1]->GetEntries() / hist[0]->GetEntries() << " ); "
	    << "  r_tf[0] = new RooRealVar( \"r_tf0\", \"r_tf0\", " << std::setw(10) << std::right << hist[ 2]->GetEntries() / hist[0]->GetEntries() << " ); "
	    << "  r_f [0] = new RooRealVar( \"r_f0\",  \"r_f0\",  " << std::setw(10) << std::right << hist[ 3]->GetEntries() / hist[0]->GetEntries() << " );"  << std::endl
	    << "  r_tt[1] = new RooRealVar( \"r_tt1\", \"r_tt1\", " << std::setw(10) << std::right << hist[ 5]->GetEntries() / hist[4]->GetEntries() << " ); "
	    << "  r_tf[1] = new RooRealVar( \"r_tf1\", \"r_tf1\", " << std::setw(10) << std::right << hist[ 6]->GetEntries() / hist[4]->GetEntries() << " ); "
	    << "  r_f [1] = new RooRealVar( \"r_f1\",  \"r_f1\",  " << std::setw(10) << std::right << hist[ 7]->GetEntries() / hist[4]->GetEntries() << " );"  << std::endl
	    << "  r_tt[2] = new RooRealVar( \"r_tt2\", \"r_tt2\", " << std::setw(10) << std::right << hist[ 9]->GetEntries() / hist[8]->GetEntries() << " ); "
	    << "  r_tf[2] = new RooRealVar( \"r_tf2\", \"r_tf2\", " << std::setw(10) << std::right << hist[10]->GetEntries() / hist[8]->GetEntries() << " ); "
	    << "  r_f [2] = new RooRealVar( \"r_f2\",  \"r_f2\",  " << std::setw(10) << std::right << hist[11]->GetEntries() / hist[8]->GetEntries() << " );"  << std::endl;
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] add_cut;
  delete   c1;
    
  return 0;
}

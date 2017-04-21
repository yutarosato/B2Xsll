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


using namespace sig_gmc_rd_cut2;
//using namespace sig_gmc_rd_cut2_cc;
//using namespace Mbc_comb;
using namespace Mbc_bkg;

const Bool_t flag_scale    = true;
const Bool_t flag_scan     = true;
const Bool_t flag_k4pi     = true; // 1(veto  K4pi    modes)
const Bool_t flag_unflavor = true; // 1(veto unflavor modes)
const Bool_t flag_ccpi0    = true; // 1(veto ccpi0 peak    )
const Bool_t flag_save     = true; // outfile.eps and outfile.root

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==3 || argc==4) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_Mbc (char*)stream (int)used_nstream [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Char_t*  stream       = argv[1];
  Double_t used_nstream = atof(argv[2]); 
  Int_t    fl_appRun   = 1;
  if( argc==4 ) fl_appRun = atoi( argv[3] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t fl_5body_veto    = !true;
  const Int_t Nchain           =  2;   // [ee(sigmc),mm(sigmc)]
  const Int_t Nctgry           = 14;   // [seven q^2] x [two coslp]
  const Int_t Nplot            = 14+1; // [seven q^2] x [two coslp] + total
  const Int_t Ntmp             =  2*Nctgry; // [ee,mm      ]
  const Int_t Nhist            =  3*Nplot;  // [ee,mm,ee+mm]
  const Int_t nfile[Nchain]    = {0};
  const Int_t fl_sb[Nchain]    = {0,0}; // 0(bkg), 1(sig)
  const Int_t fl_mode_ll[Ntmp] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,
				  0,0,0,0,0,0,0,0,0,0,0,0,0,0,}; // 1(e), 0(mu)
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t* LRNB_cut = new Char_t[4096];
  LRNB_cut = "1";
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<3; i++ ){
    infile[i] = new Char_t[1024];
    if     ( fl_sb[i]==0 ) sTmp << indir[fl_sb[i]] << "gMC_*_s0[" << stream << "]";     // bkg
    //if     ( fl_sb[i]==0 ) sTmp << indir[fl_sb[i]] << "CC_";                              // cc mc // tmppppp
    //else if( fl_sb[i]==1 ) sTmp << indir[fl_sb[i]] << "sigMC_*_set[" << setname << "]"; // sig
    //else if( fl_sb[i]==2 ) sTmp << indir[fl_sb[i]] << "RD_";                            // rd
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Ntmp] ={
    {1},
    {0,1},
    {0,0,1},
    {0,0,0,1},
    {0,0,0,0,1},
    {0,0,0,0,0,1},
    {0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {1,1,1,1,1,1,1,1,1,1,1,1,1,1},

    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1},

    {1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
  };
  
  Char_t** add_cut = new Char_t*[Ntmp];
  for( Int_t i=0; i<Ntmp; i++ ) add_cut[i] = new Char_t[4096];
  for( Int_t i=0; i<Ntmp; i++ ){
    sTmp << "  (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999))"; // charmonium events
    if     ( i%Nctgry== 0 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 1,  1 ).c_str();
    else if( i%Nctgry== 1 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 1, -1 ).c_str();
    else if( i%Nctgry== 2 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 2,  1 ).c_str();
    else if( i%Nctgry== 3 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 2, -1 ).c_str();
    else if( i%Nctgry== 4 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 3,  1 ).c_str();
    else if( i%Nctgry== 5 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 3, -1 ).c_str();
    else if( i%Nctgry== 6 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 5,  1 ).c_str();
    else if( i%Nctgry== 7 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 5, -1 ).c_str();
    else if( i%Nctgry== 8 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 7,  1 ).c_str();
    else if( i%Nctgry== 9 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 7, -1 ).c_str();
    else if( i%Nctgry==10 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 8,  1 ).c_str();
    else if( i%Nctgry==11 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 8, -1 ).c_str();
    else if( i%Nctgry==12 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 9,  1 ).c_str();
    else if( i%Nctgry==13 ) sTmp << " && " << LRNB_cut << " && " << makeCut_q2(fl_mode_ll[i], 9, -1 ).c_str();
    
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
  TCanvas*  c1      = Canvas( "c1","c1", 4, 4 );
  TCanvas*  c2      = Canvas( "c2","c2", 1, 1 );
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
    if( flag_scan ){
      std::cout << "+++++++++++++++++++++++++++++++++++++++++++++< SCAN tmphist" << j << " : "
		<< chain[fl_mode_ll[j]]->GetTree()->GetEntries( Form("5.27<%s && 5.29>%s && %s", axis, axis, add_cut[j]) )
		<< " entries >++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		<< std::endl;
      chain[fl_mode_ll[j]]->GetTree()->SetScanField(0);
      chain[fl_mode_ll[j]]->GetTree()->Scan(
					    "exprun:event:rm_xs:Mbc:cc_m-cc_morg:cc_m:cc_morg:de:lpgt:lmgt:lporg:lmorg:lpself:lmself:lpselfid:lmselfid:lpmoid:lmmoid:korg:pi1org:pi2org:pi3org:pi4org:pi0self:pi0gam1e:pi0gam2e:gm_ccng:gm_ccge:cc_ng:cc_g1_se:cc_g2_se:heg_self:leg_self:heg_moid:leg_moid:cc_mheg:cc_mleg:cc_morgh:cc_morgl:gb1_semi:gb2_semi:gb1nd:gb2nd:gb1d1_se:gb1d2_se:gb2d1_se:gb2d2_se:gm_bg1:gm_bg2:gm_b1:gm_b2:gm_l1:gm_l2:gm_nu1:gm_nu2:rest:rest_sw:rest2:rest2_sw:dntrk:ntrk",
					    Form("5.27<%s && 5.29>%s && %s", axis, axis, add_cut[j]) );
      std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }
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
    if( flag_scale ) hist[i]->Scale( 1/used_nstream );
    if     ( i%Nplot== 0 ) hist[i]->SetTitle( "(q2,coslp)=(1,+)" );
    else if( i%Nplot== 1 ) hist[i]->SetTitle( "(q2,coslp)=(1,-)" );
    else if( i%Nplot== 2 ) hist[i]->SetTitle( "(q2,coslp)=(2,+)" );
    else if( i%Nplot== 3 ) hist[i]->SetTitle( "(q2,coslp)=(2,-)" );
    else if( i%Nplot== 4 ) hist[i]->SetTitle( "(q2,coslp)=(3,+)" );
    else if( i%Nplot== 5 ) hist[i]->SetTitle( "(q2,coslp)=(3,-)" );
    else if( i%Nplot== 6 ) hist[i]->SetTitle( "(q2,coslp)=(5,+)" );
    else if( i%Nplot== 7 ) hist[i]->SetTitle( "(q2,coslp)=(5,-)" );
    else if( i%Nplot== 8 ) hist[i]->SetTitle( "(q2,coslp)=(7,+)" );
    else if( i%Nplot== 9 ) hist[i]->SetTitle( "(q2,coslp)=(7,-)" );
    else if( i%Nplot==10 ) hist[i]->SetTitle( "(q2,coslp)=(8,+)" );
    else if( i%Nplot==11 ) hist[i]->SetTitle( "(q2,coslp)=(8,-)" );
    else if( i%Nplot==12 ) hist[i]->SetTitle( "(q2,coslp)=(9,+)" );
    else if( i%Nplot==13 ) hist[i]->SetTitle( "(q2,coslp)=(9,-)" );
    else if( i%Nplot==14 ) hist[i]->SetTitle( "Total"            );
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
  c2->Draw();
  for(Int_t i=0; i<Nhist; i++ ){
    if( i%Nplot==Nplot-1 ){
      c2->cd();
      if     ( i/Nplot==0 ) hist[i]->Draw("P");     // ee
      else if( i/Nplot==1 ) hist[i]->Draw("Psame"); // mm
      //else if( i/Nplot==2 ) hist[i]->Draw("Psame"); // ee+mm
    }else{
      c1->cd(i%Nplot+1);
      if     ( i/Nplot==0 ) hist[i]->Draw("P");     // ee
      else if( i/Nplot==1 ) hist[i]->Draw("Psame"); // mm
      //else if( i/Nplot==2 ) hist[i]->Draw("Psame"); // ee+mm
    }
  }
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[0+0*Nplot],"ee",        "P" );
  legend1->AddEntry( hist[0+1*Nplot],"#mu#mu",    "P" );
  legend1->AddEntry( hist[0+2*Nplot],"ee+#mu#mu", "P" );
  legend1->Draw();

  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_peak_cc_s0%s_c1.eps",  axis, stream) );
    c2->Print( Form("pic/%s_peak_cc_s0%s_c2.eps",  axis, stream) );
    TFile outfile( Form("pic/%s_peak_cc_s0%s.root", axis, stream), "RECREATE" );
    c1->Write();
    c2->Write();
    for( Int_t i=0; i<Nhist; i++ ) hist[i]->Write();
    outfile.Close();
  }

  // +++++++++++++++ Log for Roofit +++++++++++++++++++++
  std::cout << " +++++++++++++++ Log for Roofit +++++++++++++++++++++ "  << std::endl
	    << "  RooRealVar*** r_peak_cc = new RooRealVar**[Nroohist];" << std::endl
	    << "  for( Int_t i=0; i<Nroohist; i++ ) r_peak_cc[i] = new RooRealVar*[Nbin_afb+1];" << std::endl;
  for( Int_t i=0; i<3; i++ ){
    for( Int_t j=0; j<Nplot; j++ ){
      Int_t k = j + i*Nplot;
      std::cout << "  r_peak_cc"
		<< "[" << std::setw(2) << std::right << i << "]"
		<< "[" << std::setw(2) << std::right << j << "] "
		<< "= new RooRealVar( "
		<< "\"r_peak_cc" << k << "\", "
		<< "\"r_peak_cc" << k << "\", "
		<< hist[k]->Integral()
		<< " );" << std::endl;
    }
  }
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  
  return 0;
}

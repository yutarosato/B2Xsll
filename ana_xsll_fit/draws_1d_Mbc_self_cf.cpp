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
  const Int_t Nbin_afb = 14;
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
  
  Char_t*** add_cut = new Char_t**[Nbin_afb];
  for( Int_t i=0; i<Nbin_afb; i++ ){
    add_cut[i] = new Char_t*[Ntmp];
    for( Int_t j=0; j<Ntmp; j++ ){
      Int_t k = j+i*Ntmp;
      add_cut[i][j] = new Char_t[4096];
      
      //for( Int_t i=0; i<Ntmp; i++ ){
      if     ( j%Nctgry==0 ) sTmp << LRNB_cut << " && self==1";   // true
      else if( j%Nctgry==1 ) sTmp << LRNB_cut << " && self!=1 &&" // false with (t,t)
				  << makeCut_q2fl( 1, 1 );
      else if( j%Nctgry==2 ) sTmp << LRNB_cut << " && self!=1 &&" // false with (t,f)
				  << makeCut_q2fl( 1, 0 );
      else if( j%Nctgry==3 ) sTmp << LRNB_cut << " && self!=1 &&" // false with (f,t)
				  << makeCut_q2fl( 0, 1 );
      else if( j%Nctgry==4 ) sTmp << LRNB_cut << " && self!=1 &&" // false with (f,f)
				  << makeCut_q2fl( 0, 0 );
      

      if     ( i ==  0 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 1,  1 ).c_str();
      else if( i ==  1 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 1, -1 ).c_str();
      else if( i ==  2 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 2,  1 ).c_str();
      else if( i ==  3 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 2, -1 ).c_str();
      else if( i ==  4 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 3,  1 ).c_str();
      else if( i ==  5 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 3, -1 ).c_str();
      else if( i ==  6 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 5,  1 ).c_str();
      else if( i ==  7 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 5, -1 ).c_str();
      else if( i ==  8 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 7,  1 ).c_str();
      else if( i ==  9 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 7, -1 ).c_str();
      else if( i == 10 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 8,  1 ).c_str();
      else if( i == 11 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 8, -1 ).c_str();
      else if( i == 12 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 9,  1 ).c_str();
      else if( i == 13 ) sTmp << " && " << makeCut_q2( fl_mode_ll[j], 9, -1 ).c_str();
      
      if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
      if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
      if( flag_ccpi0    ) sTmp << " && " << cut_ccpi0;
      strcpy( add_cut[i][j], (char*)sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
    }
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D***   tmphist = new TH1D** [Nbin_afb]; //[Ntmp ]
  TH1D***   hist    = new TH1D** [Nbin_afb]; //[Nhist]
  for( Int_t i=0; i<Nbin_afb; i++ ){
    tmphist[i] = new TH1D*[Ntmp];
    hist   [i] = new TH1D*[Nhist];
  }
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
  for( Int_t i=0; i<Nbin_afb; i++ ){
    for( Int_t j=0; j<Ntmp; j++ ){
      Int_t k = j+i*Ntmp;
      std::cout << Form( "<tmphist %d > add_cut : ", k ) << add_cut[i][j] << std::endl;
      tmphist[i][j] = new TH1D( Form("tmphist%d",k), Form("%s",chain[fl_mode_ll[j]]->GetChange()), xbin,offset+xmin,offset+xmax );
      chain[fl_mode_ll[j]]->GetTree()->Project( Form("tmphist%d",k), axis, add_cut[i][j] );
    }
  }
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist(plot) *************************************" << std::endl;
  for( Int_t i=0; i<Nbin_afb; i++ ){
    for( Int_t j=0; j<Nhist; j++ ){
      Int_t k = j+i*Nhist;
      hist[i][j] = new TH1D( Form("hist%d",k),Form("hist%d",k), xbin,offset+xmin,offset+xmax );
      Deco( hist[i][j], 1, j/Nplot+2, j/Nplot+2 );
      for( Int_t m=0; m<Ntmp; m++ ){
	if( add[j][m] ) hist[i][j]->Add( tmphist[i][m] );
      }
      hist[i][j]->Sumw2();
      sTmp << i+1 << " q^{2}- cos#theta bin, ";
      if     ( j%Nplot==0 ) sTmp << "true";
      else if( j%Nplot==1 ) sTmp << "(q2,fl)=(t,t)";
      else if( j%Nplot==2 ) sTmp << "(q2,fl)=(t,f)";
      else if( j%Nplot==3 ) sTmp <<  "(q2,fl)=(f,?)";
      hist[i][j]->SetTitle( sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
    }
  }
  
  
  // +++++++ display ++++++++++++++++++++++++++++++++++
  Double_t entry_sig   [Nbin_afb][Nhist] = {0}; // # of events in signal box region
  Double_t entry_canvas[Nbin_afb][Nhist] = {0}; // # of events in     fit    region
  for( Int_t i=0; i<Nbin_afb; i++ ){
    for( Int_t j=0; j<Nhist; j++ ){
      Int_t k = j+i*Nhist;
      std::cout << Form("<hist [%d][%d] > ",i,j);
      Double_t entry_all    = hist[i][j]->GetEntries();
      entry_canvas[i][j]    = hist[i][j]->Integral();
      Double_t entry_under  = hist[i][j]->GetBinContent(0);
      Double_t entry_over   = hist[i][j]->GetBinContent(xbin+1);
      for( Int_t m=hist[i][j]->FindBin(5.27+0.000001); m<=xbin; m++ ) entry_sig[i][j] += hist[i][j]->GetBinContent(j);
      TString sTmp = "Added-Files( ";
      for( Int_t m=0; m<Ntmp; m++ ) if( add[j][m] ) std::cout << k <<  ",";
      std::cout <<  ")["
		<< entry_all
		<< " events ( canvas : "
		<< entry_canvas[i][j]
		<< " / under : "
		<< entry_under
		<< " / over  : "
		<< entry_over
		<< " / sig  : "
		<< entry_sig[i][j]
		<< "]"
		<< std::endl;
      if( flag_norm ) hist[i][j]->Scale( 1/hist[i][j]->Integral() );
    }
  }
  std::cout << std::endl;
  
  Double_t entry_sig_each[Nbin_afb][Ntmp] = {0}; // # of events in signal box region (tmphist)
  for( Int_t i=0; i<Nbin_afb; i++ ){
    for( Int_t j=0; j<Ntmp; j++ ){
      std::cout << Form("<tmphist [%d][%d]> ",i,j);
      for( Int_t m=tmphist[i][j]->FindBin(5.27+0.000000001); m<=xbin; m++ ) entry_sig_each[i][j] += tmphist[i][j]->GetBinContent(j);
      std::cout << tmphist[i][j]->Integral() << " events(canvas), "
		<< entry_sig_each[i][j]      << " events(signal-box)" << std::endl;
    }
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  TCanvas** c1 = new TCanvas*[Nplot]; // true, (t,t), (t,f), (f,?)
  for( Int_t i=0; i<Nplot; i++ ){
    c1[i] = Canvas( Form("c1_%d",i+1), Form("c1_%d",i+1), 4, 4 );
    c1[i]->Draw();
  }

  for( Int_t i=0; i<Nbin_afb; i++ ){
    for(Int_t j=0; j<Nhist; j++ ){
      c1[j%Nplot]->cd(i+1);
      if     ( j/Nplot==0 ) hist[i][j]->Draw("P"); // ee
      else if( j/Nplot==1 ) hist[i][j]->Draw("Psame"); // mm
      //else if( j/Nplot==2 ) hist[i][j]->Draw("Psame"); // ee+mm
    }
  }

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[0][0+0*Nplot],"ee",        "P" );
  legend1->AddEntry( hist[0][0+1*Nplot],"#mu#mu",    "P" );
  legend1->AddEntry( hist[0][0+2*Nplot],"ee+#mu#mu", "P" );
  legend1->Draw();

  for( Int_t i=0; i<Nplot; i++ ) c1[i]->Update();

  if( flag_save ){
    for( Int_t i=0; i<Nplot; i++ ) c1[i]->Print( Form("pic/%s_self_cf_set%s_%d.eps",  axis, setname, i+1) );
    TFile outfile( Form("pic/%s_self_cf_set%s.root", axis, setname), "RECREATE" );
    for( Int_t i=0; i<Nplot; i++ ) c1[i]->Write();
    for( Int_t i=0; i<Nbin_afb; i++ ){
      for( Int_t j=0; j<Nhist; j++ ){
	hist[i][j]->Write();
      }
    }
    outfile.Close();
  }


  // +++++++++++++++ Log for Roofit +++++++++++++++++++++
  Double_t total_evt_true = 0;
  Double_t total_evt_tt   = 0;
  Double_t total_evt_tf   = 0;
  Double_t total_evt_f    = 0;
  
  for( Int_t i=0; i<3; i++ ){
    for( Int_t j=0; j<Nbin_afb; j++ ){
      total_evt_true += entry_canvas[j][i*Nplot+0];
      total_evt_tt   += entry_canvas[j][i*Nplot+1];
      total_evt_tf   += entry_canvas[j][i*Nplot+2];
      total_evt_f    += entry_canvas[j][i*Nplot+3];
    }
  }
  std::cout << "TOTAL (true) : " << total_evt_true << std::endl
	    << "TOTAL (t,t ) : " << total_evt_tt   << std::endl
	    << "TOTAL (t,f ) : " << total_evt_tf   << std::endl
	    << "TOTAL (f,? ) : " << total_evt_f    << std::endl;

  std::cout << " +++++++++++++++ Log for Roofit +++++++++++++++++++++ "  << std::endl
	    << "  RooRealVar*** r_tt = new RooRealVar**[Nroohist];"      << std::endl
	    << "  RooRealVar*** r_tf = new RooRealVar**[Nroohist];"      << std::endl
	    << "  RooRealVar*** r_f  = new RooRealVar**[Nroohist];"      << std::endl
	    << "  for( Int_t i=0; i<Nroohist; i++ ){"     << std::endl
	    << "    r_tt[i] = new RooRealVar*[Nbin_afb];" << std::endl
	    << "    r_tf[i] = new RooRealVar*[Nbin_afb];" << std::endl
	    << "    r_f[i]  = new RooRealVar*[Nbin_afb];" << std::endl
	    << "  }"                                      << std::endl;

  for( Int_t i=0; i<3; i++ ){
    for( Int_t j=0; j<Nbin_afb/2; j++ ){
      Int_t k = 2*j+i*Nbin_afb;
      Float_t ratio_tt_f = entry_canvas[2*j+0][i*Nplot+1] / entry_canvas[2*j+0][i*Nplot+0]; // Forward (t,t)/Forward (true)
      Float_t ratio_tt_b = entry_canvas[2*j+1][i*Nplot+1] / entry_canvas[2*j+1][i*Nplot+0]; // Backward(t,t)/Backward(true)
      Float_t ratio_tf_f = entry_canvas[2*j+0][i*Nplot+2] / entry_canvas[2*j+1][i*Nplot+0]; // Forward (t,f)/Backward(true)
      Float_t ratio_tf_b = entry_canvas[2*j+1][i*Nplot+2] / entry_canvas[2*j+0][i*Nplot+0]; // Backward(t,f)/Forward (true)
      Float_t ratio_f_f  = entry_canvas[2*j+0][i*Nplot+3] / total_evt_true;                 // Forward (f,?)/Total   (true)
      Float_t ratio_f_b  = entry_canvas[2*j+1][i*Nplot+3] / total_evt_true;                 // Backward(f,?)/Total   (true) 
      
      std::cout << "  r_tt[" << std::setw(2) << std::right << i << "][" << std::setw(2) << std::right << 2*j+0 << "] = new RooRealVar( " << std::setw(11) << std::left << Form("\"r_tt%d\",",k+0) << std::setw(11) << std::left << Form(" \"r_tt%d\",",2*j+0) << std::setw(10) << std::right << ratio_tt_f << "); "
		<< "  r_tf[" << std::setw(2) << std::right << i << "][" << std::setw(2) << std::right << 2*j+0 << "] = new RooRealVar( " << std::setw(11) << std::left << Form("\"r_tf%d\",",k+0) << std::setw(11) << std::left << Form(" \"r_tf%d\",",2*j+0) << std::setw(10) << std::right << ratio_tf_f << "); "
		<< "  r_f [" << std::setw(2) << std::right << i << "][" << std::setw(2) << std::right << 2*j+0 << "] = new RooRealVar( " << std::setw(11) << std::left << Form("\"r_f%d\",", k+0) << std::setw(11) << std::left << Form(" \"r_f%d\",", 2*j+0) << std::setw(10) << std::right << ratio_f_f  << "); "
		<< std::endl
		<< "  r_tt[" << std::setw(2) << std::right << i << "][" << std::setw(2) << std::right << 2*j+1 << "] = new RooRealVar( " << std::setw(11) << std::left << Form("\"r_tt%d\",",k+1) << std::setw(11) << std::left << Form(" \"r_tt%d\",",2*j+1) << std::setw(10) << std::right << ratio_tt_b << "); "
		<< "  r_tf[" << std::setw(2) << std::right << i << "][" << std::setw(2) << std::right << 2*j+1 << "] = new RooRealVar( " << std::setw(11) << std::left << Form("\"r_tf%d\",",k+1) << std::setw(11) << std::left << Form(" \"r_tf%d\",",2*j+1) << std::setw(10) << std::right << ratio_tf_b << "); " 
		<< "  r_f [" << std::setw(2) << std::right << i << "][" << std::setw(2) << std::right << 2*j+1 << "] = new RooRealVar( " << std::setw(11) << std::left << Form("\"r_f%d\",", k+1) << std::setw(11) << std::left << Form(" \"r_f%d\",", 2*j+1) << std::setw(10) << std::right << ratio_f_b  << "); "
		<< std::endl; 
    }
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  
  return 0;
}

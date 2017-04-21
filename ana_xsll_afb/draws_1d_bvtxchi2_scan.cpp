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
  //using namespace sig_gmc_rd_cut2_beforebgsup;
  using namespace sig_gmc_rd_cut3_beforebgsup;
  //using namespace sig_gmc_rd_emu_beforebgsup;
  using namespace bgsup_var_org;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_var (int)fl_factor (char*)stream (double)used_nstream [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Int_t    fl_axis       = 2; // chi2(B-vtx), fixed
  Int_t    fl_factor     = atoi(argv[1]); // 0.001*fl_factor
  Char_t*  stream        = argv[2];
  Double_t used_nstream  = atof(argv[3]); 
  Int_t   fl_appRun      = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Bool_t flag_calib    = true;
  const Bool_t flag_k4pi     = !true; // 1(veto  K4pi    modes)
  const Bool_t flag_unflavor = !true; // 1(veto unflavor modes)
  const Bool_t flag_save     = true; // outfile.eps and outfile.root
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t    Nchain             = 2*2; // [gmc,rd] x [ee,mm]
  const Int_t    Nhist              = 2*3; // [gmc,rd] x [ee,mm,ee+mm]
  const Int_t    nfile[Nchain]      = {0};
  const Int_t    fl_sb[Nchain]      = {0,2,0,2}; // 0(bkg), 1(sig), 2(rd)
  const Int_t    fl_mode_ll[Nchain] = {1,1,  // 1(ee)
				       //0,0}; // 0(mm)
				       1,1}; // 0(mm) // tmppppp
  
  const Double_t scale_event_bkg = used_nstream; // gmc : N -> N/alpha
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t* LRNB_cut = new Char_t[4096];
  LRNB_cut = "1";
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if     ( fl_sb[i]==0 ) sTmp << indir[fl_sb[i]] << "gMC_*_s0[" << stream << "]";     // bkg
    else if( fl_sb[i]==1 ) sTmp << indir[fl_sb[i]] << "sigMC_";                         // sig
    else if( fl_sb[i]==2 ) sTmp << indir[fl_sb[i]] << "RD_";                            // rd
    else                   std::cerr << "[ABORT] Wrong fl_sb : " << fl_sb[i] << std::endl, abort();
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1},       // gmc, ee
    {0,1},     // rd,  ee
    {0,0,1},   // gmc, mm
    {0,0,0,1}, // rd,  mm
    {1,0,1,0}, // gmc, ee+mm
    {0,1,0,1}, // rd,  ee+mm
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    sTmp << LRNB_cut;
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    /*
      sTmp << "&& ( "//  with   pi0 modes
		      << Form( "(rm_xs==1001 && abs(xs_m-%f)<0.050) || ", PDGmass::kstrp ) // K+pi0
		      << Form( "(rm_xs==1010 && abs(xs_m-%f)<0.050)",     PDGmass::kstr0 ) // Kspi0
		      << " ) ";
    */
    /*
    sTmp << "&& ( " // without pi0 modes (K, Ks, Kpi+, Kspi+)
	 << Form( " (rm_xs==  1) || " ) // K+
	 << Form( " (rm_xs== 10) || " ) // Ks
	 << Form( " (rm_xs==101 && abs(xs_m-%f)<0.050) || ", PDGmass::kstr0 ) // K+pi-
	 << Form( " (rm_xs==110 && abs(xs_m-%f)<0.050) ",    PDGmass::kstrp ) // Kspi-
	 << " )";
    */
    strcpy( add_cut[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",3, 1 );
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
    if     ( fl_sb[j] )      chain[j]->GetCut()->Set( "dzll3d",  1, 0.0200        );// rd
    //if     ( fl_sb[j] )      chain[j]->GetCut()->Set( "dzll3dorg",  1, 0.0200        );// rd
    else                     chain[j]->GetCut()->Set( "dzll3dorg",  1, 0.01811594203 );// mm(mc)
    //if     ( fl_sb[j] )      chain[j]->GetCut()->Set( "dzll3dorg",  1, 0.0150        );// rd
    //else if( fl_mode_ll[j] ) chain[j]->GetCut()->Set( "dzll3dorg",  1, 0.01134644478 );// ee(mc)
    //else                     chain[j]->GetCut()->Set( "dzll3dorg",  1, 0.01213592233 );// mm(mc)
    //if     ( fl_sb[j] )      chain[j]->GetCut()->Set( "dzll3dorg",  1, 0.0200        );// rd
    //else if( fl_mode_ll[j] ) chain[j]->GetCut()->Set( "dzll3dorg",  1, 0.01512859304 );// ee(mc)
    //else                     chain[j]->GetCut()->Set( "dzll3dorg",  1, 0.01618122977 );// mm(mc)
    chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.26 );
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
    std::cout << "add_cut : " << add_cut[j] << std::endl;
  }
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
    //if( flag_calib && fl_sb[j]!=2 ) chain[j]->GetTree()->Project( Form("tmphist%d",j), Form("(1+0.001*%d)*%s",  fl_factor, axis[fl_axis]), add_cut[j] ); // chi2(b-vtx)
    //else                            chain[j]->GetTree()->Project( Form("tmphist%d",j), axis[fl_axis],                                      add_cut[j] );
    if( flag_calib && fl_sb[j]!=2 ){
      if( fl_sb[j]==2 ) chain[j]->GetTree()->Project( Form("tmphist%d",j), Form("(1+0.001*%d)*%s",  fl_factor, axisrd[fl_axis]), add_cut[j] ); // chi2(b-vtx), rd
      else              chain[j]->GetTree()->Project( Form("tmphist%d",j), Form("(1+0.001*%d)*%s",  fl_factor, axis  [fl_axis]), add_cut[j] ); // chi2(b-vtx), mc
    }else{
      if( fl_sb[j]==2 ) chain[j]->GetTree()->Project( Form("tmphist%d",j), axisrd[fl_axis], add_cut[j] ); // rd
      else              chain[j]->GetTree()->Project( Form("tmphist%d",j), axis  [fl_axis], add_cut[j] ); // mc
    }
    
    tmphist[j]->Sumw2();
    //if( fl_sb[j]==0 ) tmphist[j]->Scale( 1/scale_event_bkg );
  }
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
    Deco( hist[i], 2, i%2+1, i%2+1 );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
    //hist[i]->Scale( 1/hist[i]->GetEntries() );
    hist[i]->Scale( 1/hist[i]->Integral() );
  }
  // +++++++ display ++++++++++++++++++++++++++++++++++
  std::cout << "=================================" << std::endl
	    << "[Scale Factor]"                    << std::endl
	    << "bkg : " << scale_event_bkg         << std::endl
	    << "=================================" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin[fl_axis]+1);

    for( Int_t j=0; j<Nchain; j++ ) if( add[i][j] ) std::cout << j << ",";
    std::cout <<  ")[ "
	      << std::setw(12) << std::right << entry_all
	      <<  " events ( canvas : "
	      << std::setw(12) << std::right << std::setprecision(3) << entry_canvas
	      << " / under : "
	      << std::setw(12) << std::right << std::setprecision(3) << entry_under
	      << " / over  : "
	      << std::setw(12) << std::right << std::setprecision(3) << entry_over
	      << "]"
	      << std::endl;

  }

  // +++++++ calculation of chi2 for scale factor ++++++++++++++++++++
  Double_t diff_chi2[3] = {0}; // [ee,mm,ee+mm]
  Double_t diff_ndf[3]  = {0}; // [ee,mm,ee+mm]
  for( Int_t k=0; k<xbin[fl_axis]; k++ ){
    Double_t tmp_delta1 = hist[1]->GetBinContent(k+1) - hist[0]->GetBinContent(k+1);
    Double_t tmp_error1 = sqrt( hist[1]->GetBinError(k+1)*hist[1]->GetBinError(k+1) + hist[0]->GetBinError(k+1)*hist[0]->GetBinError(k+1) );
    Double_t tmp_delta2 = hist[3]->GetBinContent(k+1) - hist[2]->GetBinContent(k+1);
    Double_t tmp_error2 = sqrt( hist[3]->GetBinError(k+1)*hist[3]->GetBinError(k+1) + hist[2]->GetBinError(k+1)*hist[2]->GetBinError(k+1) );
    Double_t tmp_delta3 = hist[5]->GetBinContent(k+1) - hist[4]->GetBinContent(k+1);
    Double_t tmp_error3 = sqrt( hist[5]->GetBinError(k+1)*hist[5]->GetBinError(k+1) + hist[4]->GetBinError(k+1)*hist[4]->GetBinError(k+1) );
    if( tmp_error1 ){
      diff_chi2[0] += tmp_delta1*tmp_delta1/tmp_error1/tmp_error1;
      diff_ndf[0]++;
    }
    if( tmp_error2 ){
      diff_chi2[1] += tmp_delta2*tmp_delta2/tmp_error2/tmp_error2;
      diff_ndf[1]++;
    }
    if( tmp_error3 ){
      diff_chi2[2] += tmp_delta3*tmp_delta3/tmp_error3/tmp_error3;
      diff_ndf[2]++;
    }
  }
  std::cout << "Difference_chi2(  ee ) : " << fl_factor << " : " << diff_chi2[0]/diff_ndf[0] << " : " << diff_chi2[0] << " / " << diff_ndf[0] << std::endl;
  std::cout << "Difference_chi2(  mm ) : " << fl_factor << " : " << diff_chi2[1]/diff_ndf[1] << " : " << diff_chi2[1] << " / " << diff_ndf[1] << std::endl;
  std::cout << "Difference_chi2(ee+mm) : " << fl_factor << " : " << diff_chi2[2]/diff_ndf[2] << " : " << diff_chi2[2] << " / " << diff_ndf[2] << std::endl;
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();

  TH2D** waku = new TH2D*[3];

  if( fl_axis==3 ){
    waku[0] = Waku( Nhist, hist, xlabel[fl_axis], Form("ee",   xlabel[fl_axis]), "ee",        2 );
    waku[1] = Waku( Nhist, hist, xlabel[fl_axis], Form("mm",   xlabel[fl_axis]), "#mu#mu",    2 );
    waku[2] = Waku( Nhist, hist, xlabel[fl_axis], Form("eemm", xlabel[fl_axis]), "ee+#mu#mu", 2 );
  }else{
    waku[0] = Waku( Nhist, hist, xlabel[fl_axis], Form("ee",   xlabel[fl_axis]), "ee"        );
    waku[1] = Waku( Nhist, hist, xlabel[fl_axis], Form("mm",   xlabel[fl_axis]), "#mu#mu"    );
    waku[2] = Waku( Nhist, hist, xlabel[fl_axis], Form("eemm", xlabel[fl_axis]), "ee+#mu#mu" );
  }

  for(Int_t i=0; i<Nhist; i++ ){
    c1->cd(i/2+1);
    if( fl_axis==3 ) gPad->SetLogy(); // for CL(B-vtx)
    if( i%2==0 ) waku[i/2]->Draw();
    hist[i]->Draw("same");
  }

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[0], "gmc", "L" );
  legend1->AddEntry( hist[1], "rd",  "L" );
  legend1->Draw();

  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_factor%d_s0%s.eps",  fname[fl_axis], fl_factor, stream) );
    c1->Print( Form("pic/%s_factor%d_s0%s.root", fname[fl_axis], fl_factor, stream) );
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


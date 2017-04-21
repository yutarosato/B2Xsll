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
#include <TFile.h>
#include <TArrow.h>
#include <TMath.h>


Int_t main( Int_t argc, Char_t** argv ){
  using namespace sig_gmc_rd_cut2;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==6 || argc==7) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_Mbc (char*)stream (char*)setname (double)nstream (double)nset (int)fl_body [(int)fl_appRun]"
					<< std::endl, abort();
  Char_t*  stream       = argv[1];
  Char_t*  setname      = argv[2];
  Double_t used_nstream = atof(argv[3]);
  Double_t used_nset    = atof(argv[4]);
  Int_t    fl_body      = atoi(argv[5]);
  Int_t    fl_appRun    = 1;
  if( argc==7 ) fl_appRun = atoi( argv[6] );

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Ncategory     = 6; //[gmc(qq,non-peak,peak),sigmc(false except(t,t), false with (t,t), true)]
  const Int_t Nplot         = 4; //[bkg(qq,non-peak,peak),sigmc]
  const Int_t Nchain        = Ncategory*2; // x[ee,mm      ]
  const Int_t Nhist         = Nplot    *3; // x[ee,mm,ee+mm]
  const Int_t fl_sb[Nchain] = {0,0,0,1,1,1,
			       0,0,0,1,1,1}; // 0(bkg), 1(sig), 2(rd)
  const Int_t fl_mode_ll[Nchain] = {1,1,1,1,1,1,
				    0,0,0,0,0,0}; // 1(e), 0(mu)
  const Int_t nfile[Nchain] = {0};
  const Double_t scale_event_sig = sigmc_amount*(used_nset/(Double_t)nset); // sigmc : N -> N/alpha
  const Double_t scale_event_bkg = used_nstream;                            //   gmc : N -> N/alpha

  const Int_t add[Nhist][Nchain] ={
    {1,0,0,0,0,0}, // ee(qq)
    {1,1,0,1,0,0}, // ee(non-peak)
    {1,1,1,1,0,0}, // ee(peak)
    {1,1,1,1,1,1}, // ee(sig)
    {0,0,0,0,0,0,1,0,0,0,0,0}, // mm
    {0,0,0,0,0,0,1,1,0,1,0,0}, // mm
    {0,0,0,0,0,0,1,1,1,1,0,0}, // mm
    {0,0,0,0,0,0,1,1,1,1,1,1}, // mm
    {1,0,0,0,0,0,1,0,0,0,0,0}, // ee+mm
    {1,1,0,1,0,0,1,1,0,1,0,0}, // ee+mm
    {1,1,1,1,0,0,1,1,1,1,0,0}, // ee+mm
    {1,1,1,1,1,1,1,1,1,1,1,1}, // ee+mm
  };

  const Bool_t   flag_scale    = true;
  const Bool_t   flag_k4pi     = true; // 1(veto  K4pi    modes)
  const Bool_t   flag_unflavor = true; // 1(veto unflavor modes)

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if     ( fl_sb[i]==0 ) sTmp << indir[fl_sb[i]] << "gMC_*_s0[" << stream << "]";     // bkg
    else if( fl_sb[i]==1 ) sTmp << indir[fl_sb[i]] << "sigMC_*_set[" << setname << "]"; // sig
    else if( fl_sb[i]==2 ) sTmp << indir[fl_sb[i]] << "RD_";                            // rd
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Char_t* LRNB_cut = new Char_t[2048];
  strcpy( LRNB_cut, "1" ); // lrnb cut should be already applied.
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    sTmp << LRNB_cut;
    if     ( i%Ncategory==0 ){ // gmc(qq)
      sTmp << " && genbfl==0";
    }else if( i%Ncategory==1 ){ // gmc(bb-non-peak)
      sTmp << " && genbfl!=0"
	   << " && !( (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999)) )" // charmonium events(right-pi0)
	   << " && !(!(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3) "  // double miss-id
           << " && !(!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)" // single miss-id (cc)
	   << " && !(!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)"; // single miss-id with neutrino
    }else if( i%Ncategory==2 ){ // gmc(bb-peak)
      sTmp << " && genbfl!=0";
      sTmp << " && ("
	   << "    ( (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999)) )" // charmonium events(right-pi0)
	   << " || (!(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3) "  // double miss-id
	   << " || (!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)" // single miss-id (cc)
	   << " || (!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)"; // single miss-id wit
      sTmp << " )";
    }else if( i%Ncategory==3 ){ // false except (t,t)
      sTmp << LRNB_cut << " && self!=1 && ";
      sTmp << makeCut_q2fl( 1,1,1 );
    }else if( i%Ncategory==4 ){ // false with (t,t)
      sTmp << LRNB_cut << " && self!=1 && ";
      sTmp << makeCut_q2fl( 1,1 );
    }else if( i%Ncategory==5 ){ // true
      sTmp << LRNB_cut << " && self==1 ";
    }

    sTmp << " && " << makeCut_body( fl_body ).c_str();
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    strcpy( add_cut[i], (Char_t*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
    
  using namespace Mbc_bkg;
  //using namespace Mbc_bkg_wide;

  const Bool_t flag_save = true; // outfile.eps and outfile.root
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1", 3, 1 );

  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll[j] )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 );
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
    if( flag_scale ){
      /*
      tmphist[j]->Sumw2();
      if     ( fl_sb[j]==1 ) tmphist[j]->Scale( 1/scale_event_sig );
      else if( fl_sb[j]==0 ) tmphist[j]->Scale( 1/scale_event_bkg );
      */
      if( fl_sb[j]==1 ) tmphist[j]->Scale( scale_event_bkg/scale_event_sig );
      tmphist[j]->Sumw2();
      tmphist[j]->Scale( 1/scale_event_bkg );
    }
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    Deco( hist[i], 3, col_fil[i%Nplot], col_fil[i%Nplot] );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }

  // +++++++ display ++++++++++++++++++++++++++++++++++
  Double_t entry_sig[Nhist] = {0}; // # of events in signal box region
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
    for( Int_t j=hist[i]->FindBin(5.27+0.000000001); j<=xbin; j++ ) entry_sig[i] += hist[i]->GetBinContent(j);
    std::cout << "Added-Files( ";
    for( Int_t j=0; j<Nchain; j++ ) if( add[i][j] ) std::cout << j <<  ",";
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

  Double_t entry_sig_each[Nchain] = {0}; // # of events in signal box region (tmphist)
  for( Int_t i=0; i<Nchain; i++ ){
    std::cout << Form("<tmphist %d> ",i);
    for( Int_t j=tmphist[i]->FindBin(5.27+0.000000001); j<=xbin; j++ ) entry_sig_each[i] += tmphist[i]->GetBinContent(j);
    std::cout << entry_sig_each[i] << " events" << std::endl;
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  TH2D** waku = new TH2D*[Nhist];
  waku[0] = Waku( Nhist, hist, Form("%s (ee)",       xlabel) );
  waku[1] = Waku( Nhist, hist, Form("%s (#mu#mu)",   xlabel) );
  waku[2] = Waku( Nhist, hist, Form("%s (ee+#mu#mu)",xlabel) );
  
  for(Int_t i=Nhist-1; i>=0; i-- ){
    c1->cd(i/Nplot+1);
    if     ( 0*Nplot <= i && i < 1*Nplot ) hist[i]->SetTitle( Form("ee"       ) );
    else if( 1*Nplot <= i && i < 2*Nplot ) hist[i]->SetTitle( Form("#mu#mu"   ) );
    else if( 2*Nplot <= i && i < 3*Nplot ) hist[i]->SetTitle( Form("ee+#mu#mu") );

    hist[i]->SetXTitle(xlabel);
    hist[i]->SetYTitle(waku[i/Nplot]->GetYaxis()->GetTitle());
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    //waku[i]->Draw();

    if( i%Nplot==Nplot-1 ) hist[i]->Draw("hist");
    else                           hist[i]->Draw("hist same");
  }

  c1->Update();
  c1->Print( Form("pic/%s_%dbody_s0%s_set%s.eps", axis, fl_body, stream, setname) );
  
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] add_cut;
  delete   c1;

  return 0;
}


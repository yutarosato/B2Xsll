#include <iostream>
#include <math.h>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Set_fake_lepton/Nominal_cut_selection.h"
#include "../Set_fake_lepton/Branch.h"
#include "../Set_fake_lepton/makeCut.h"

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
#include <TF1.h>
#include <TArrow.h>
#include <TLeaf.h>

Int_t main( Int_t argc, Char_t** argv ){
  using namespace gmc_fl;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==2 || argc==3) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_Mbc_bb (int)fl_peak [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  const Int_t fl_peak   = atoi(argv[1]);
  Char_t*  stream       = "0";
  Double_t used_nstream = 1;
  Int_t    fl_appRun    = 1;
  if( argc==3 ) fl_appRun = atoi( argv[2] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t  Nchain        = 2; // [e,mu]
  const Int_t  Nhist         = 2; // [e,mu]
  const Int_t  nfile[Nchain] = {0};
  const Bool_t fl_weight     = true;
  const Bool_t flag_k4pi     = !true; // 1(veto  K4pi    modes)
  const Bool_t flag_unflavor = !true; // 1(veto unflavor modes)
  const Int_t  fl_q2         = 0;
  const Int_t  fl_mode_ll[Nchain] = {1, // 1(e)
				     0, // 0(mu)
  };
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    sTmp << indir << "/gMC_*s0[" << stream << "]";
    //sTmp << indir << "/gMC_*e031*s0[" << stream << "]"; // tmppp
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1},    // ee
    {0,1},  // mm
  };

  Char_t* LRNB_cut = new Char_t[2048];
  strcpy( LRNB_cut, "1" ); // lrnb cut should be already applied.
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    sTmp << "genbfl!=0 && "; // select only bb(mixed, charged) events

    ///*
    if     ( fl_peak== 0 ) sTmp << "  (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999))"; // charmonium events(right-pi0)
    else if( fl_peak== 1 ) sTmp << " !(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 ";  // double miss-id
    else if( fl_peak== 2 ) sTmp << " !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg "; // single miss-id (cc)
    else if( fl_peak== 3 ) sTmp << " !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg "; // single miss-id with neutrino
    //*/
    /*
    if     ( fl_peak== 0 ) sTmp << "( (lpgt==3 && lmgt==3 ) && lpself==1 && lmself==1 && rest_sw==0 && gm_bg1<10 ) "; // charmonium events
    else if( fl_peak== 1 ) sTmp << " !(lpgt==3 && lmgt==3 ) && lpself==0 && lmself==0                                                         && rest_sw==0 && gm_bg1<10 && lporg==lmorg "; // double miss-id
    else if( fl_peak== 2 ) sTmp << " !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<10              && ( abs(korg)==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg "; // single miss-id (cc)
    else if( fl_peak== 3 ) sTmp << " !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<10 && rm_xs>999 && ( abs(korg)==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg "; // single miss-id with neutrino
    */
    sTmp << " && " << LRNB_cut;
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    strcpy( add_cut[i], (Char_t*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  using namespace q2_theta_nonuniform;
  const Bool_t flag_save = true; // outfile.eps and outfile.root

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain    = new MChain*[Nchain];
  TH2D**    tmphist  = new TH2D*  [Nchain];
  TH2D**    tmphistE = new TH2D*  [Nchain];
  TH2D**    hist     = new TH2D*  [Nhist];
  TCanvas*  c1       = Canvas( "c1","c1", 2, 1 );
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
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 );
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
    tmphist[j]  = new TH2D( Form("tmphist%d", j), Form("%s",chain[j%Nchain]->GetChange()), xbin,xbins[fl_mode_ll[j]], ybin,ymin,ymax );
    tmphistE[j] = new TH2D( Form("tmphistE%d",j), Form("%s",chain[j%Nchain]->GetChange()), xbin,xbins[fl_mode_ll[j]], ybin,ymin,ymax );

    if( fl_weight ){
      chain[j]->GetTree()->Project( Form("tmphist%d", j), axis, Form("(%s)*(%s)", add_cut[j], "weight"         ) );
      chain[j]->GetTree()->Project( Form("tmphistE%d",j), axis, Form("(%s)*(%s)", add_cut[j], "weightE*weightE") );
      for ( Int_t kx=0; kx<xbin; kx++ ){
	for ( Int_t ky=0; ky<ybin; ky++ ){
	  tmphist[j]->SetBinError(kx,ky, sqrt(tmphistE[j]->GetBinContent(kx,ky)) );
	}
      }

    }else chain[j]->GetTree()->Project( Form("tmphist%d", j), axis, Form("%s",add_cut[j]) );

  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH2D( Form("hist%d", i), Form("lep%d",fl_mode_ll[i]), xbin,xbins[fl_mode_ll[i]], ybin,ymin,ymax );
    ((TGaxis*)hist[i]->GetXaxis())->SetMaxDigits(3);
    ((TGaxis*)hist[i]->GetYaxis())->SetMaxDigits(3);
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    hist[i]->SetXTitle( xlabel );
    hist[i]->SetYTitle( ylabel );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j], (Double_t) add[i][j] );
    }
  }

  // +++++++ display ++++++++++++++++++++++++++++++++++
   for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
    for( Int_t j=0; j<Nchain; j++ ) if( add[i][j] ) std::cout << j << ",";
    std::cout << ")["
	      <<  entry_all
	      <<  " events ( canvas : "
	      << entry_canvas
	      << " / under : "
	      << entry_under
	      << " / over  : "
	      << entry_over
	      << "]" << std::endl;
  }

  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  for( Int_t i=0; i<Nchain; i++ ){
    c1->cd(i+1);
    hist[i]->Draw("COLZ");
  }
  
  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/peak_q2_theta_peak%d.eps",  fl_peak)  );
    c1->Print( Form("pic/peak_q2_theta_peak%d.root", fl_peak) );
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

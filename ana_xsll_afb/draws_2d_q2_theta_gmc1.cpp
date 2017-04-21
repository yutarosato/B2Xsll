#include <iostream>

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
#include <TH2D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TArrow.h>


Int_t main( Int_t argc, Char_t** argv ){
  using namespace gmc;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_2d (char*)stream (double)used_nstream (int)fl_mode_ll [(int)fl_appRun]" << std::endl
					<< "[  stream  ] 0,1,2,3,4,5,01,0-5" << std::endl, abort();
  
  const Char_t*  stream       = argv[1];
  const Double_t used_nstream = atof(argv[2]); 
  const Int_t    fl_mode_ll   = atoi(argv[3]); // 1(e), 0(mu)
  Int_t          fl_appRun    = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t    Nchain          = 2; // signal-region and 1 Mbc sideband
  const Int_t    Ncategory       = 6; // qq, 4peaks, no-peak
  const Int_t    Ntmp            = Nchain*Ncategory;
  const Int_t    Nplot           = 3; // qq peak, no-peak
  const Int_t    Nhist           = Nchain*Nplot;
  const Int_t    nfile[Nchain]   = {0};
  const Double_t scale_event_bkg = used_nstream; // gmc : N -> N/alpha
  const Bool_t   flag_save       = true;
  const Bool_t   flag_k4pi       = true; // 1(veto  K4pi    modes)
  const Bool_t   flag_unflavor   = true; // 1(veto unflavor modes)
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    sTmp << indir << "/gMC_*_s0[" << stream << "]";
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Ntmp] = {
    { 1, 0, 0, 0, 0, 0}, // qq
    { 0, 1, 1, 1, 1, 0}, // bb (   peak)
    { 0,-1,-1,-1,-1, 1}, // bb (nonpeak)
    { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, // qq
    { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0}, // bb (   peak)
    { 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1, 1}, // bb (nonpeak)
  };
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t fl_q2 = 0;
  Char_t* LRNB_cut = new Char_t[2048];
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.91, 0.36, 0.94, 0.92, 0.87, 0.57, 0.92, 0.89).c_str() ); // 2d NB_lep(bcs=bb)
  LRNB_cut = "1";
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;

  Char_t** add_cut = new Char_t*[Ncategory];
  for( Int_t j=0; j<Ncategory; j++ ){
    add_cut[j] = new Char_t[4096];

    if     ( j%Ncategory== 0 ) sTmp << "genbfl==0 "; // qq
    if     ( j%Ncategory== 1 ) sTmp << "genbfl!=0 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999))"; // charmonium events
    else if( j%Ncategory== 2 ) sTmp << "genbfl!=0 &&  !(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 ";  // double miss-id
    else if( j%Ncategory== 3 ) sTmp << "genbfl!=0 && !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg "; // swtich miss-id
    else if( j%Ncategory== 4 ) sTmp << "genbfl!=0 && !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg "; // single miss-id with neutrino
    else if( j%Ncategory== 5 ) sTmp << "genbfl!=0 "; // total-bb

    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    sTmp << " && " << LRNB_cut;
    strcpy( add_cut[j], (Char_t*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  using namespace q2_theta_nonuniform;
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH2D**    tmphist = new TH2D*  [Ntmp];
  TH2D**    hist2   = new TH2D*  [Nhist];
  TH1D**    histx   = new TH1D*  [Nhist];
  TH1D**    histy   = new TH1D*  [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",Nchain+2, Nplot );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis) << std::endl;
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    if     ( j==1 ) chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.26 ); // Mbc sideband
  }
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j ) << std::endl;
    chain[j]->GetCut()->Display(0);
    chain[j]->MakeTree();
  }

  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  std::cout << " ************************ make tmphist *************************************" << std::endl;
  for( Int_t j=0; j<Ntmp; j++ ){
    //tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[j%Nchain]->GetChange()), xbin,xmin,xmax, ybin,ymin,ymax );
    tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[j%Nchain]->GetChange()), xbin,xbins[fl_mode_ll], ybin,ymin,ymax );
    chain[j/Ncategory]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j%Ncategory] ); 
    tmphist[j]->Sumw2();
    tmphist[j]->Scale( 1/scale_event_bkg );
    // +++++++ display ++++++++++++++++++++++++++++++++++
    std::cout << std::setw(10) << Form("<infile%d -> tmphist%d>",j/Ncategory,j);
    std::cout << std::setw(12) << tmphist[j]->GetEntries() << " events : " 
	      << add_cut[j%Ncategory]  << std::endl;
  }

  // +++++++ make hist2  ++++++++++++++++++++++++++++++++++
  std::cout << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    //hist2[i] = new TH2D( Form("hist%d",i),Form("hist%d",i), xbin,xmin,xmax, ybin,ymin,ymax );
    hist2[i] = new TH2D( Form("hist%d",i),Form("hist%d",i), xbin,xbins[fl_mode_ll], ybin,ymin,ymax );
    if     ( i== 0 ) hist2[i]->SetTitle("qq,           signal-box" );
    else if( i== 1 ) hist2[i]->SetTitle("bb(peak),     signal-box" );
    else if( i== 2 ) hist2[i]->SetTitle("bb(non-peak), signal-box" );
    else if( i== 3 ) hist2[i]->SetTitle("qq,           sideband"   );
    else if( i== 4 ) hist2[i]->SetTitle("bb(peak),     sideband"   );
    else if( i== 5 ) hist2[i]->SetTitle("bb(non-peak), sideband"   );
    ((TGaxis*)hist2[i]->GetXaxis())->SetMaxDigits(3);
    ((TGaxis*)hist2[i]->GetYaxis())->SetMaxDigits(3);
    hist2[i]->GetXaxis()->CenterTitle();
    hist2[i]->GetYaxis()->CenterTitle();
    hist2[i]->SetXTitle( xlabel );
    hist2[i]->SetYTitle( ylabel );
    for( Int_t j=0; j<Ntmp; j++ ){
      if( add[i][j] ) hist2[i]->Add( tmphist[j], (Double_t) add[i][j] );
    }
    // +++++++ display ++++++++++++++++++++++++++++++++++
    std::cout << std::setw(10) << Form("<hist%d>",i);
    //std::cout << std::setw(12) << hist2[i]->GetEntries() << " events, Added-Files( ";
    std::cout << std::setw(12) << hist2[i]->Integral() << " events, Added-Files( ";
    for( Int_t j=0; j<Ntmp; j++ ) if( add[i][j] ) std::cout << add[i][j]*j << ",";
    std::cout << " )" << std::endl; 
  }

  // +++++++ make hist[xy] ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nhist; i++ ){
    histx[i] = hist2[i]->ProjectionX();
    histy[i] = hist2[i]->ProjectionY();
    histx[i]->Scale( 1/histx[i]->Integral() );
    histy[i]->Scale( 1/histy[i]->Integral() );
    Deco( histx[i], 0, i/Nplot+1, i/Nplot+1 );
    Deco( histy[i], 0, i/Nplot+1, i/Nplot+1 );
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  Double_t allymax_x = GetYMax( Nhist, histx);
  Double_t allymax_y = GetYMax( Nhist, histy);
  for( Int_t i=0; i<Nhist; i++ ){
    histx[i]->SetMaximum( 1.1*allymax_x );
    histy[i]->SetMaximum( 1.1*allymax_y );
  }
  c1->Draw();
  for(Int_t i=0; i<Nplot; i++ ){
    for(Int_t j=0; j<Nchain; j++ ){
      c1->cd((Nchain+2)*i+1+j); hist2[j*Nplot+i]->Draw( "COLZ" );
      c1->cd((Nchain+2)*i+1+Nchain+0);
      if( j==0 ) histx[j*Nplot+i]->Draw();
      else       histx[j*Nplot+i]->Draw("same");
      c1->cd((Nchain+2)*i+1+Nchain+1);
      if( j==0 ) histy[j*Nplot+i]->Draw();
      else       histy[j*Nplot+i]->Draw("same");
    }
  }

  // +++++++ tlegend ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( histx[0*Nplot], "signal-box", "L" );
  legend1->AddEntry( histx[1*Nplot], "sideband",   "L" );

  legend1->Draw();
  legend1->Draw();

  // +++++++ save +++++++++++++++++++++++++++++++++++++++
  c1->Update();

  if( flag_save ) c1->Print( Form("pic/2d_q2_theta_gmc1_lep%d_s0%s.eps", fl_mode_ll, stream) );
  std::cout << "finish" << std::flush;  
  if( fl_appRun ) app.Run();
  
  delete[] chain;
  delete[] tmphist;
  delete[] hist2;
  delete[] histx;
  delete[] histy;
  delete   c1;
  delete   legend1;
  
  return 0;
}

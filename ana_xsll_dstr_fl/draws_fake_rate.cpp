#include <iostream>
#include <math.h>
#include <string.h>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Set_dstr_fl/Nominal_cut_selection.h"
#include "../Set_dstr_fl/Branch.h"

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
#include <TPaletteAxis.h>

Int_t main( Int_t argc, Char_t** argv ){
  using namespace gmc;
  using namespace p_cos;
  TApplication app( "app", &argc, argv );
  Style();
  const Bool_t fl_appRun = true;
  const Bool_t flag_save = true; // outfile.eps and outfile.root
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Char_t* infile  = "test.dat";
  const Int_t   Nhist   = 4; // [mu-K, mu-pi, e-K, e-pi]
  const Int_t   fl_adef = 1; // definition of momentum-axis (0) 0.4 GeV, (1) 0.6 GeV
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << "axis definition : " << fl_adef << std::endl;
  TH2D** hist_fake = new TH2D*  [Nhist];
  TH2D** hist_err  = new TH2D*  [Nhist];
  TH2D** hist_chi2 = new TH2D*  [Nhist+2];
  hist_fake[0] = new TH2D( "hist0_fake", "fake K(mu-id)",  xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_fake[1] = new TH2D( "hist1_fake", "fake K(e-id)",   xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_fake[2] = new TH2D( "hist2_fake", "fake pi(mu-id)", xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_fake[3] = new TH2D( "hist3_fake", "fake pi(e-id)",  xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_err [0] = new TH2D( "hist0_err ", "err  K(mu-id)",  xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_err [1] = new TH2D( "hist1_err ", "err  K(e-id)",   xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_err [2] = new TH2D( "hist2_err ", "err  pi(mu-id)", xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_err [3] = new TH2D( "hist3_err ", "err  pi(e-id)",  xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_chi2[0] = new TH2D( "hist0_chi2", "chi2 K(mu-id)",  xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_chi2[1] = new TH2D( "hist1_chi2", "chi2 K(e-id)",   xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_chi2[2] = new TH2D( "hist2_chi2", "chi2 pi(mu-id)", xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_chi2[3] = new TH2D( "hist3_chi2", "chi2 pi(e-id)",  xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_chi2[4] = new TH2D( "hist4_chi2", "chi2 K(no-pid)", xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
  hist_chi2[5] = new TH2D( "hist5_chi2", "chi2 pi(no-pid)",xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );

  for( Int_t i=0; i<Nhist; i++ ){
    hist_fake[i]->GetXaxis()->CenterTitle();
    hist_fake[i]->GetYaxis()->CenterTitle();
    hist_fake[i]->SetXTitle( xlabel );
    hist_fake[i]->SetYTitle( ylabel );
    hist_fake[i]->SetLabelSize( 0.03, "Z" );
    hist_fake[i]->GetZaxis()->SetTitleOffset(-1.3);
    hist_err [i]->GetXaxis()->CenterTitle();
    hist_err [i]->GetYaxis()->CenterTitle();
    hist_err [i]->SetXTitle( xlabel );
    hist_err [i]->SetYTitle( ylabel );
    hist_err [i]->SetLabelSize( 0.03, "Z" );
    hist_err [i]->GetZaxis()->SetTitleOffset(-1.3);
  }
  for( Int_t i=0; i<Nhist+2; i++ ){
    hist_chi2[i]->GetXaxis()->CenterTitle();
    hist_chi2[i]->GetYaxis()->CenterTitle();
    hist_chi2[i]->SetXTitle( xlabel );
    hist_chi2[i]->SetYTitle( ylabel );
    hist_chi2[i]->SetLabelSize( 0.03, "Z" );
    hist_chi2[i]->GetZaxis()->SetTitleOffset(-1.3);
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t   buf[1024]; 
  ifstream fin;
  fin.open(infile);
  while(!fin.eof()){
    fin.getline(buf,1024);
    if(buf[0]=='#') continue; // comment
    if(buf[0]=='*') break;    // finish
    std::istringstream sTmp(buf);
    Int_t    fl_pid, fl_par, fl_mom, fl_cos;
    Double_t fake_rate, fake_rateE, chi2;
    sTmp >> fl_pid >> fl_par >> fl_mom >> fl_cos >> chi2 >> fake_rate >> fake_rateE;

    // select signal lepton region
    if( fl_adef==1 ){ // wide(0.6 GeV)
      if     ( fl_cos >=1 && fl_cos <= 3 && fl_mom>=4 ) continue;
      else if( fl_cos >=4 && fl_cos <= 6 && fl_mom>=5 ) continue;
      else if( fl_cos >=7 && fl_cos <= 8 && fl_mom>=6 ) continue;
      else if( fl_cos >=9 && fl_cos <=10 && fl_mom>=7 ) continue;

    }else if( fl_adef== 0 ){ // narrow(0.4 GeV)
      if     ( fl_cos >= 1 && fl_cos <= 2 && fl_mom>= 5 ) continue;
      else if( fl_cos >= 3 && fl_cos <= 4 && fl_mom>= 6 ) continue;
      else if( fl_cos >= 5 && fl_cos <= 6 && fl_mom>= 7 ) continue;
      else if( fl_cos >= 7 && fl_cos <= 8 && fl_mom>= 8 ) continue;
      else if( fl_cos == 9                && fl_mom>= 9 ) continue;
      else if( fl_cos >=10                && fl_mom>=10 ) continue;
    }
    
    if( fl_pid !=2 ){
      hist_fake[fl_pid+2*fl_par]->SetBinContent( fl_cos, fl_mom, fake_rate  );
      hist_fake[fl_pid+2*fl_par]->SetBinError  ( fl_cos, fl_mom, fake_rateE );
      hist_err [fl_pid+2*fl_par]->SetBinContent( fl_cos, fl_mom, fake_rateE );
      hist_chi2[fl_pid+2*fl_par]->SetBinContent( fl_cos, fl_mom, chi2       );
    }else{
      hist_chi2[   4    +fl_par]->SetBinContent( fl_cos, fl_mom, chi2       );
    }
  }
    fin.close();

  // +++++++ draw ++++++++++++++++++++++++++++++++++

  TCanvas*  c1 = Canvas( "c1","c1", 2, 2 );
   c1->Draw();
  
  for(Int_t i=0; i<Nhist; i++ ){
    c1->cd(i+1);
    hist_fake[i]->DrawCopy( "COLZ" );
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  TCanvas*  c2 = Canvas( "c2","c2", 2, 2 );
  c2->Draw();
  
  for(Int_t i=0; i<Nhist; i++ ){
    c2->cd(i+1);
    hist_err [i]->DrawCopy( "COLZ" );
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  TCanvas*  c3 = Canvas( "c3","c3", 2, 3 );
  c3->Draw();
  
  for(Int_t i=0; i<Nhist+2; i++ ){
    c3->cd(i+1);
    hist_chi2[i]->DrawCopy( "COLZ" );
  }

  c1->Update();
  c2->Update();
  c3->Update();

  if( flag_save ){
    c1->Print( "pic/fake_rate_c1.eps"  );
    c2->Print( "pic/fake_rate_c2.eps"  );
    c3->Print( "pic/fake_rate_c3.eps"  );
    c1->Print( "pic/fake_rate_c1.root" );
    c2->Print( "pic/fake_rate_c2.root" );
    c3->Print( "pic/fake_rate_c3.root" );
    TFile outfile( "pic/fake_rate.root", "RECREATE" );
    for( Int_t i=0; i<Nhist;   i++ ) hist_fake[i]->Write();
    for( Int_t i=0; i<Nhist;   i++ ) hist_err [i]->Write();
    for( Int_t i=0; i<Nhist+2; i++ ) hist_chi2[i]->Write();
    outfile.Close();
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] hist_fake;
  delete[] hist_err ;
  delete[] hist_chi2;
  delete   c1;
  delete   c2;
  delete   c3;

  return 0;
}

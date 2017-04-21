#include <iostream>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Util/Manip.h"
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
  //using namespace sigmc_afterbgsup;
  using namespace sigmc_afterbgsup_c10sym;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_2d (int)fl_mode_ll (int)fl_xsid (char*)setname [(int)fl_appRun]" << std::endl
					<< "[fl_mode_ll] : 1(e), 0(mu)"                << std::endl
					<< "[ fl_xsid  ] : 0(all), 1(K), 2(K*), 3(Xs)" << std::endl, abort();
  
  Int_t   fl_mode_ll = atoi( argv[1] ); // 1(e), 0(mu)
  Int_t   fl_xsid    = atoi( argv[2] ); // 0(all), 1(K), 2(K*), 3(Xs)
  Char_t* setname    = argv[3]; 
  Int_t   fl_appRun  = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t  Nchain        = 2; // fixed [h511, h12]
  const Int_t  Nhist         = 3; // fixed [h511, h12, h511/h12]
  const Int_t  nfile[Nchain] = {0};
  
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) infile[i] = new Char_t[1024]; 
  
  std::stringstream sTmp;
  sTmp << indir << "sigMC_*_m9999m_*_set[" << setname << "]";
  strcpy( infile[0], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  sTmp << indir_gen << "sigMC_*_m9999m_*_set[" << setname << "]";
  strcpy( infile[1], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];

  sTmp << "self==1 && recbfl!=2";
  if     ( fl_xsid==1 ){ // K
    sTmp << "&& gm_fl_xs<0 && gm_m_xs<0.50"
	 << " && ( "
	 << Form( " (rm_xs==  1) || " ) // K+
	 << Form( " (rm_xs== 10) "    ) // Ks
	 << " ) ";
  }else if( fl_xsid==2 ){ // K*
    sTmp << "&& gm_fl_xs<0 && gm_m_xs>0.50"
	 << " && ( "
	 << Form( " (rm_xs==101 && abs(xs_m-%f)<0.050) || ", PDGmass::kstr0 ) // K+pi-
	 << Form( " (rm_xs==110 && abs(xs_m-%f)<0.050) || ", PDGmass::kstrp ) // Kspi-
	 << Form( "(rm_xs==1001 && abs(xs_m-%f)<0.050) || ", PDGmass::kstrp ) // K+pi0
	 << Form( "(rm_xs==1010 && abs(xs_m-%f)<0.050) ",    PDGmass::kstr0 ) // Kspi0
	 << " ) ";
  }else if( fl_xsid==3 ) sTmp << "&& gm_fl_xs>0 && gm_m_xs>1.10"; // Xs
  
  strcpy( add_cut[0], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();

  if     ( fl_xsid==0 ) add_cut[1] = "1"; // h12 (all)
  else if( fl_xsid==1 ) add_cut[1] = "1&&gm_fl_xs<0&&Xs_m<0.50"; // h12 (K),  gm_fl_xs : 1(xs), -1(k,k*), gm_xs, gm_bg, Xs_m
  else if( fl_xsid==2 ) add_cut[1] = "1&&gm_fl_xs<0&&Xs_m>0.50"; // h12 (K*),
  else if( fl_xsid==3 ) add_cut[1] = "1&&gm_fl_xs>0&&Xs_m>1.10"; // h12 (Xs),  
  
  const Bool_t  flag_save     = true;
  
  const Char_t*  tname[2] = { "h511", "h12" };
  const Char_t*  axis[2]  = { "coslp:cc_m*cc_m", "coslp:llg_m*llg_m"};
  const Int_t    ybin     =    24;
  const Double_t ymin     = -1.00;
  const Double_t ymax     =  1.00;
  const Char_t*  ylabel   = "cos#theta";
  const Double_t offset   =   0.0;
  const Double_t xmin     =   0.0;
  const Double_t xmax     =  24.0;
  const Char_t*  xlabel   = "q^{2} [GeV^{2}]";
  const Int_t    xbin     = 9;
  const Double_t xbins[2][xbin+1] = { // [mm,ee]
  {0.0, 2.0, 4.3,
  (PDGmass::jpsi -0.25-0.0*0.15)*(PDGmass::jpsi -0.25-0.0*0.15),
  (PDGmass::jpsi +0.10+0.0*0.05)*(PDGmass::jpsi +0.10+0.0*0.05),
  (PDGmass::psi2s-0.15-0.0*0.10)*(PDGmass::psi2s-0.15-0.0*0.10),
  (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ),
  16.0, 19.0, 24.0},
  {0.0, 2.0, 4.3,
  (PDGmass::jpsi -0.25-1.0*0.15)*(PDGmass::jpsi -0.25-1.0*0.15),
  (PDGmass::jpsi +0.10+1.0*0.05)*(PDGmass::jpsi +0.10+1.0*0.05),
  (PDGmass::psi2s-0.15-1.0*0.10)*(PDGmass::psi2s-0.15-1.0*0.10),
  (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ),
  16.0, 19.0, 24.0},
  };
  
  const Int_t nslicex = 4;
  const Int_t slicex[2*nslicex] = {1,2,
  3,3,
  5,5,
  7,9};
  const Int_t nslicey = 4;
  const Int_t slicey[2*nslicey] = {0*ybin/nslicey+1, 1*ybin/nslicey,
  1*ybin/nslicey+1, 2*ybin/nslicey,
  2*ybin/nslicey+1, 3*ybin/nslicey,
  3*ybin/nslicey+1, 4*ybin/nslicey,
  };
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain     = new MChain*[Nchain];
  TH2D**    tmphist   = new TH2D*  [Nchain];
  TH2D**    hist      = new TH2D*  [Nhist];
  TH1D**    histx_rec = new TH1D*  [nslicey];
  TH1D**    histx_gen = new TH1D*  [nslicey];
  TH1D**    histx_eff = new TH1D*  [nslicey];
  TH1D**    histy_rec = new TH1D*  [nslicex];
  TH1D**    histy_gen = new TH1D*  [nslicex];
  TH1D**    histy_eff = new TH1D*  [nslicex];
  TCanvas*  c1        = Canvas( "c1","c1",4, 3 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++

  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname[j],axis[j]) << std::endl;
    std::cout << Form( "<infile %d > ", j );
    if( j== 1 ) chain[j] = new MChain( infile[1], tname[j], branch_table(), nfile[j], tail );
    else        chain[j] = new MChain( infile[0], tname[j], branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname[j] );
  }
  
  // ++++++++++++++++++++++++
  // cut change
  //
  //chain[0]->GetCut()->Set(    441, 0 );
  //chain[0]->GetCut()->Set(    443, 0 );
  //chain[0]->GetCut()->Set( 100443, 0 );
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
  
  //chain[0]->GetTree()->Scan("cc_m:cc_morg:rm_l:rm_xs:coslp:lppt:lmpt:mpp:mmp:epp:emp:lpc:lmc", "self==1&&cc_m*cc_m<1.0 && (coslp>0.6||coslp<-0.6)" ); // tmppppp

  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin,xbins[fl_mode_ll], ybin,ymin,ymax );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis[j], add_cut[j] );
  }
  
  // +++++++ make hist[h511,h12]  ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nchain; i++ ){
    hist[i] = new TH2D( Form("hist_%s",tname[i]),Form("hist_%s",tname[i]), xbin,xbins[fl_mode_ll], ybin,ymin,ymax );
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    hist[i]->SetXTitle( xlabel );
    hist[i]->SetYTitle( ylabel );
    hist[i]->Add( tmphist[i] );
    // +++++++ display ++++++++++++++++++++++++++++++++++
    std::cout << std::setw(10) << Form("<hist%d>",i);
    std::cout << std::setw(9)  << hist[i]->GetEntries() << " events"
	      << std::endl; 
  }
  
  // +++++++ make hist[h511/h12]  ++++++++++++++++++
  hist[2] = divide( hist[0], hist[1], "h511_h12", "h511_h12" );
  hist[2]->SetAxisRange( 0.0,0.20,"Z" );

  // +++++++ projection hist ++++++++++++++++++

  std::cout << Form("[projection -> X axis(%s)]", xlabel) << std::endl;
  for( Int_t my=0; my<nslicey; my++ ){
    histx_rec[my] = new TH1D(*(hist[0]->ProjectionX( Form("_px%d",my+1), slicey[2*my], slicey[2*my+1], "e" )) );
    histx_gen[my] = new TH1D(*(hist[1]->ProjectionX( Form("_px%d",my+1), slicey[2*my], slicey[2*my+1], "e" )) );
    histx_rec[my]->SetLineColor  (my+1);
    histx_gen[my]->SetLineColor  (my+1);
    histx_rec[my]->SetMarkerColor(my+1);
    histx_gen[my]->SetMarkerColor(my+1);
    std::cout << Form( "%d-%d bin : ", slicey[2*my], slicey[2*my+1] )
	      << "(rec) " << histx_rec[my]->GetEntries() << ", "
	      << "(gen) " << histx_gen[my]->GetEntries() <<  std::endl;

    histx_eff[my] = divide( histx_rec[my], histx_gen[my] );
  }

  std::cout << Form("[projection -> Y axis(%s)]", ylabel) << std::endl;  
  for( Int_t mx=0; mx<nslicex; mx++ ){
    histy_rec[mx] = new TH1D(*(hist[0]->ProjectionY( Form("_py%d",mx+1), slicex[2*mx], slicex[2*mx+1], "e")) );
    histy_gen[mx] = new TH1D(*(hist[1]->ProjectionY( Form("_py%d",mx+1), slicex[2*mx], slicex[2*mx+1], "e")) );
    histy_rec[mx]->SetLineColor  (mx+1);
    histy_gen[mx]->SetLineColor  (mx+1);
    histy_rec[mx]->SetMarkerColor(mx+1);
    histy_gen[mx]->SetMarkerColor(mx+1);
    std::cout << Form( "%d-%d bin : ", slicex[2*mx], slicex[2*mx+1] )
	      << "(rec) " << histy_rec[mx]->GetEntries() << ", "
	      << "(gen) " << histy_gen[mx]->GetEntries() <<  std::endl;

    histy_eff[mx] = divide( histy_rec[mx], histy_gen[mx] );
  }

  // +++++++ make hist(afb) ++++++++++++++++++
  TH1D** hist_afb = new TH1D*[2];
  hist_afb[0] = new TH1D( "afb_rec", "afb_rec", xbin, xbins[fl_mode_ll]);
  hist_afb[1] = new TH1D( "afb_gen", "afb_gen", xbin, xbins[fl_mode_ll]);

  for( Int_t k=0; k<2; k++ ){
    hist_afb[k]->SetXTitle( xlabel   );
    hist_afb[k]->SetYTitle( "A_{FB}" );
    for( Int_t m=0; m<xbin; m++ ){
      Double_t theta_p = hist[k]->Integral( m+1,m+1,ybin/2+1, ybin   );
      Double_t theta_m = hist[k]->Integral( m+1,m+1,       1, ybin/2 );
      if( theta_p+theta_m==0 ) continue;
      Double_t afb     = (theta_p - theta_m ) / (theta_p + theta_m );
      Double_t afbE    = 2/(theta_p+theta_m)/(theta_p+theta_m)*sqrt(theta_p*theta_m*(theta_p+theta_m));
      hist_afb[k]->SetBinContent( m+1, afb  );
      hist_afb[k]->SetBinError  ( m+1, afbE );
    }
    hist_afb[k]->SetLineColor  (1+k);
    hist_afb[k]->SetMarkerColor(1+k);
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);
  hist[0]->Draw( "COLZ" );
  c1->cd(2);
  hist[1]->Draw( "COLZ" );
  c1->cd(3);
  hist[2]->Draw( "COLZ" );

  c1->cd(4);
  TH2D* waku = new TH2D("AFB","A_{FB}",2,0,25,2,-0.4,0.4);
  waku->SetXTitle( xlabel   );
  waku->SetYTitle( "A_{FB}" );
  waku->GetXaxis()->CenterTitle();
  waku->GetYaxis()->CenterTitle();
  waku->Draw();
  hist_afb[0]->Draw("PE0same");
  hist_afb[1]->Draw("PE0same");
  TLegend* legend3 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend3->SetHeader( Form("lep%d, xsid%d", fl_mode_ll, fl_xsid) );
  legend3->AddEntry( histx_eff[0],"Rec.", "P" );
  legend3->AddEntry( histx_eff[1],"Gen.", "P" );
  legend3->Draw();

  c1->cd(5);
  TH2D* wakurecx = Waku( nslicey, histx_rec, xlabel, Form("%s (rec)", xlabel), Form("%s (rec)", xlabel) );
  wakurecx->Draw();
  for( Int_t my=0; my<nslicey; my++ ) histx_rec[my]->Draw("Psame");
  TLegend* legend6 = new TLegend( 0.75,0.75,0.99,0.99 );
  for( Int_t my=0; my<nslicey; my++ ) legend6->AddEntry( histx_rec[my], Form("%d-%d bin", slicey[2*my], slicey[2*my+1]), "P" );
  legend6->Draw();

  c1->cd(6);
  for( Int_t my=0; my<nslicey; my++ ) histx_gen[my]->Scale(1/histx_gen[my]->Integral());
  TH2D* wakugenx = Waku( nslicey, histx_gen, xlabel, Form("%s (gen)", xlabel), Form("%s (gen)", xlabel) );
  wakugenx->Draw();
  for( Int_t my=0; my<nslicey; my++ ) histx_gen[my]->Draw("Psame");
  TLegend* legend7 = new TLegend( 0.75,0.75,0.99,0.99 );
  for( Int_t my=0; my<nslicey; my++ ) legend7->AddEntry( histx_gen[my], Form("%d-%d bin", slicey[2*my], slicey[2*my+1]), "P" );
  legend7->Draw();

  c1->cd(7);
  TH2D* wakurecy = Waku( nslicex, histy_rec, ylabel, Form("%s (rec)", ylabel), Form("%s (rec)", ylabel) );
  wakurecy->Draw();
  for( Int_t mx=0; mx<nslicex; mx++ ) histy_rec[mx]->Draw("Psame");
  TLegend* legend4 = new TLegend( 0.75,0.75,0.99,0.99 );
  for( Int_t mx=0; mx<nslicex; mx++ ) legend4->AddEntry( histy_rec[mx], Form("%d-%d bin", slicex[2*mx], slicex[2*mx+1]), "P" );
  legend4->Draw();

  c1->cd(8);
  for( Int_t mx=0; mx<nslicex; mx++ ) histy_gen[mx]->Scale(1/histy_gen[mx]->Integral());
  TH2D* wakugeny = Waku( nslicex, histy_gen, ylabel, Form("%s (gen)", ylabel), Form("%s (gen)", ylabel) );
  wakugeny->Draw();
  for( Int_t mx=0; mx<nslicex; mx++ ) histy_gen[mx]->Draw("Psame");
  TLegend* legend5 = new TLegend( 0.75,0.75,0.99,0.99 );
  for( Int_t mx=0; mx<nslicex; mx++ ) legend5->AddEntry( histy_gen[mx], Form("%d-%d bin", slicex[2*mx], slicex[2*mx+1]), "P" );
  legend5->Draw();
  
  c1->cd(9);
  TH2D* wakux = Waku( nslicey, histx_eff, xlabel, Form("eff. of %s", xlabel), Form("eff. of %s", xlabel) );
  wakux->Draw();
  for( Int_t my=0; my<nslicey; my++ ) histx_eff[my]->Draw("Psame");
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  for( Int_t my=0; my<nslicey; my++ ) legend1->AddEntry( histx_eff[my], Form("%d-%d bin", slicey[2*my], slicey[2*my+1]), "P" );
  legend1->Draw();

  c1->cd(10);
  TH2D* wakuy = Waku( nslicex, histy_eff, ylabel, Form("eff. of %s", ylabel), Form("eff. of %s", ylabel) );
  wakuy->Draw();
  for( Int_t mx=0; mx<nslicex; mx++ ) histy_eff[mx]->Draw("Psame");
  TLegend* legend2 = new TLegend( 0.75,0.75,0.99,0.99 );
  for( Int_t mx=0; mx<nslicex; mx++ ) legend2->AddEntry( histy_eff[mx], Form("%d-%d bin", slicex[2*mx], slicex[2*mx+1]), "P" );
  legend2->Draw();

  c1->Update();
  if( flag_save ) c1->Print( Form("pic/2d_q2_theta_lep%d_xsid%d_set%s.eps", fl_mode_ll, fl_xsid, setname) );
  std::cout << "finish" << std::flush;  
  if( fl_appRun ) app.Run();
  
  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete   c1;
    
  return 0;
}


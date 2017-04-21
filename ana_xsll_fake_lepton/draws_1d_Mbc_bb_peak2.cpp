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
#include <TLeaf.h>

Int_t main( Int_t argc, Char_t** argv ){
  using namespace gmc_fl;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==1 || argc==2) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_Mbc_bb [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Char_t*  stream       = "0";
  Double_t used_nstream = 1;
  Int_t    fl_appRun    = 1;
  if( argc==2 ) fl_appRun = atoi( argv[1] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t  Ncategory     =  5;
  const Int_t  Nplot         =  4;
  const Int_t  Nchain        = 2*Ncategory; // [e,mu      ]x[Ncategory]
  const Int_t  Nhist         = 3*Nplot;     // [e,mu,total]x[  Nplot  ]
  const Int_t  nfile[Nchain] = {0};
  const Bool_t fl_weight     = true;
  const Bool_t flag_k4pi     = !true; // 1(veto  K4pi    modes)
  const Bool_t flag_unflavor = !true; // 1(veto unflavor modes)
  const Int_t  fl_q2         = 0;
  const Int_t  fl_mode_ll[Nchain] = {1,1,1,1,1, // 1(e)
				     0,0,0,0,0, // 0(mu)
  };

  //const Double_t scale_event_bkg = used_nstream; // gmc : N -> N/alpha
  //const Bool_t flag_scale = true;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<2; i++ ){
    infile[i] = new Char_t[1024];
    sTmp << indir << "/gMC_*s0[" << stream << "]";
    //sTmp << indir << "/gMC_*e0[567]*s0[" << stream << "]"; // tmppp
    //sTmp << indir << "/gMC_*e0[012]*s0[" << stream << "]"; // tmppp
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    { 1},                  // ee
    { 0, 1},               // ee
    { 0, 0, 1,},           // ee
    { 0, 0, 0, 1,},        // ee
    { 0, 0, 0, 0, 0, 1},              // mm
    { 0, 0, 0, 0, 0, 0, 1},           // mm
    { 0, 0, 0, 0, 0, 0, 0, 1},        // mm
    { 0, 0, 0, 0, 0, 0, 0, 0, 1},     // mm
    { 1, 0, 0, 0, 0, 1},                   // ee+mm
    { 0, 1, 0, 0, 0, 0, 1},                // ee+mm
    { 0, 0, 1, 0, 0, 0, 0, 1},             // ee+mm
    { 0, 0, 0, 1, 0, 0, 0, 0, 1},          // ee+mm
  };

  Char_t* LRNB_cut = new Char_t[2048];
  strcpy( LRNB_cut, "1" ); // lrnb cut should be already applied.
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    sTmp << "genbfl!=0 && "; // select only bb(mixed, charged) events
    ///*
    if     ( i%Ncategory== 0 ) sTmp << "  (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999))";
    else if( i%Ncategory== 1 ) sTmp << " !(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 ";  // double miss-id
    else if( i%Ncategory== 2 ) sTmp << " !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg "; // single miss-id (cc)
    else if( i%Ncategory== 3 ) sTmp << " !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg "; // single miss-id with neutrino
    else if( i%Ncategory== 4 ) sTmp << " 1 "; // total
    //*/
    /*
    if     ( i%Ncategory== 0 ) sTmp << "( (lpgt==3 && lmgt==3 ) && lpself==1 && lmself==1 && rest_sw==0 && gm_bg1<10 ) "; // charmonium events
    else if( i%Ncategory== 1 ) sTmp << " !(lpgt==3 && lmgt==3 ) && lpself==0 && lmself==0                                                         && rest_sw==0 && gm_bg1<10 && lporg==lmorg "; // double miss-id
    else if( i%Ncategory== 2 ) sTmp << " !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<10              && ( abs(korg)==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg "; // single miss-id (cc)
    else if( i%Ncategory== 3 ) sTmp << " !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<10 && rm_xs>999 && ( abs(korg)==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg "; // single miss-id with neutrino
    else if( i%Ncategory== 4 ) sTmp << " 1 "; // total
    */
    sTmp << " && " << LRNB_cut;
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    strcpy( add_cut[i], (Char_t*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  using namespace Mbc_comb;
  //using namespace Mbc;

  const Bool_t flag_save = true; // outfile.eps and outfile.root
  const Bool_t flag_scan = !true;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain    = new MChain*[2];
  TH1D**    tmphist  = new TH1D*  [Nchain];
  TH1D**    tmphistE = new TH1D*  [Nchain];
  TH1D**    hist     = new TH1D*  [Nhist];
  TF1**     func     = new TF1*   [Nhist];
  TCanvas*  c1       = Canvas( "c1","c1", 3, 1 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<2; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], j )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<2; j++ ){
    chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 );
    //chain[j]->GetCut()->Set( "de",  1, -0.50, 0.0, 0.50 );
  }
  //
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t j=0; j<2; j++ ){
    std::cout << Form( "<infile %d > ", j ) << std::endl;
    chain[j]->GetCut()->Display(0);
    chain[j]->MakeTree();
  }
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<tmphist %d > add_cut : ", j ) << add_cut[j] << std::endl;
    tmphist [j] = new TH1D( Form("tmphist%d", j), Form("%s",chain[fl_mode_ll[j]]->GetChange()), xbin,offset+xmin,offset+xmax );
    tmphistE[j] = new TH1D( Form("tmphistE%d",j), Form("%s",chain[fl_mode_ll[j]]->GetChange()), xbin,offset+xmin,offset+xmax );
    if( fl_weight ){
      chain[fl_mode_ll[j]]->GetTree()->Project( Form("tmphist%d", j), axis, Form("(%s)*(%s)", add_cut[j], "weight"         ) );
      chain[fl_mode_ll[j]]->GetTree()->Project( Form("tmphistE%d",j), axis, Form("(%s)*(%s)", add_cut[j], "weightE*weightE") );
      for ( Int_t k=0; k<xbin; k++ ) tmphist[j]->SetBinError(k, sqrt(tmphistE[j]->GetBinContent(k)) );
    }else chain[fl_mode_ll[j]]->GetTree()->Project( Form("tmphist%d", j), axis, Form("%s",add_cut[j]) );

    if( flag_scan ){
      std::cout << "+++++++++++++++++++++++++++++++++++++++++++++< SCAN tmphist" << j << " : "
		<< chain[fl_mode_ll[j]]->GetTree()->GetEntries( Form("5.27<%s && 5.29>%s && %s", axis, axis, add_cut[j]) )
		<< " entries >++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		<< std::endl;
      chain[fl_mode_ll[j]]->GetTree()->SetScanField(0);
      chain[fl_mode_ll[j]]->GetTree()->Scan( "rm_xs:pi0self:lpgt:lmgt:lporg:lmorg:lpself:lmself:lpselfid:lmselfid:lpmoid:lmmoid:gb1_semi:gb2_semi:gb1nd:gb2nd:gb1d1_se:gb1d2_se:gb2d1_se:gb2d2_se:gm_bg1:gm_bg2:gm_b1:gm_b2:gm_l1:gm_l2:gm_nu1:gm_nu2:rest:rest_sw:rest2:rest2_sw:dntrk:ntrk:kmoid:pi1moid:pi2moid:pi3moid:pi4moid:korg:pi1org:pi2org:pi3org:pi4org:kselfid:pi1selfi:pi2selfi:pi3selfi:pi4selfi", Form("5.27<%s && 5.29>%s && %s", axis, axis, add_cut[j] ) );
      std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    //Deco( hist[i], 3, col_fil[i%Nplot], col_fil[i%Nplot] );
    Deco( hist[i], 1, i%Nplot+2, i%Nplot+2 );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j], (Double_t) add[i][j] );
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
    entry_sig[i] = hist[i]->Integral( hist[i]->FindBin(5.27+0.0000001), xbin );
    TString sTmp = "Added-Files( ";
    for( Int_t j=0; j<Nchain; j++ ) if( add[i][j] ) sTmp += j, sTmp += ",";
    sTmp += ")[";
    sTmp += entry_all;
    sTmp += " events ( canvas : ";
    sTmp += entry_canvas;
    sTmp += " / under : ";
    sTmp += entry_under;
    sTmp += " / over  : ";
    sTmp += entry_over;
    sTmp += " / sig  : ";
    sTmp += entry_sig[i];
    sTmp += "]";
    
    hist[i]->SetTitle( sTmp.Data() );
    std::cout << sTmp.Data() << std::endl;
  }
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
	    << std::setw(40) << std::right << "B.G. source"
	    << " || "
	    << std::setw( 7) << std::right << "ee" 
	    << std::setw( 7) << std::right << "mm" 
	    << std::setw( 7) << std::right << "ee+mm" << std::endl
	    << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  for( Int_t i=0; i<Nchain/2; i++ ){
    if     ( i== 0 ) std::cout << std::setw(40) << std::right << "charmonium events ";
    else if( i== 1 ) std::cout << std::setw(40) << std::right << "double miss-pid ";
    else if( i== 2 ) std::cout << std::setw(40) << std::right << "single miss-pid(cc) ";
    else if( i== 3 ) std::cout << std::setw(40) << std::right << "single miss-pid(nu) ";
    else if( i== 4 ) std::cout << std::setw(40) << std::right << "total ";
    std::cout << " || "
	      << std::setw(7) << std::right << tmphist[          i]->Integral( tmphist[i]->FindBin(5.27+0.00000001),tmphist[i]->FindBin(5.29-0.000000001) ) // ee
	      << std::setw(7) << std::right << tmphist[Ncategory+i]->Integral( tmphist[i]->FindBin(5.27+0.00000001),tmphist[i]->FindBin(5.29-0.000000001) ) // mm
	      << std::setw(7) << std::right
	      << tmphist[i]->Integral( tmphist[i]->FindBin(5.27+0.00000001),tmphist[i]->FindBin(5.29-0.000000001) ) + tmphist[Ncategory+i]->Integral( tmphist[i]->FindBin(5.27+0.00000001),tmphist[i]->FindBin(5.29-0.000000001) ) // ee+mm
	      << std::endl;
  }
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();

  TH2D* waku = Waku( Nhist, hist, xlabel );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s (stream=%s)", axis, stream) );

  //=================================================================
  TArrow* ar1 = new TArrow( 5.27, 1.02*((TAxis*)waku->GetYaxis())->GetXmax(),
			    5.29, 1.02*((TAxis*)waku->GetYaxis())->GetXmax(),
			    0.01,"<|>" );
  ar1->SetLineColor(2);
  ar1->SetFillColor(2);
  ar1->SetLineWidth(2);
  //=================================================================

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend** legend1 = new TLegend*[3];
  for( Int_t i=0; i<3; i++ ){
    legend1[i] = new TLegend( 0.75,0.75,0.99,0.99 );
    if     ( i== 0 ) legend1[i]->SetHeader("ee"       );
    else if( i== 1 ) legend1[i]->SetHeader("#mu#mu"   );
    else if( i== 2 ) legend1[i]->SetHeader("ee+#mu#mu");
    legend1[i]->AddEntry( hist[Nplot*i+3],"single miss-id from pi0", "P" );
    legend1[i]->AddEntry( hist[Nplot*i+2],"single miss-id from cc",  "P" );
    legend1[i]->AddEntry( hist[Nplot*i+1],"double miss-id",          "P" );
    legend1[i]->AddEntry( hist[Nplot*i+0],"charmonium",              "P" );

    c1->cd(i+1);
    waku      ->Draw();
    ar1       ->Draw();
    legend1[i]->Draw();
    for(Int_t j=Nplot-1; j>=0; j-- ) hist[Nplot*i+j]->Draw("same");
  }
  
  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_peak2_s0%s.eps",  axis, stream) );
    c1->Print( Form("pic/%s_peak2_s0%s.root", axis, stream) );
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

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
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TF1.h>
#include <TArrow.h>
#include <TPaveStats.h>
Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==8 || argc==9) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_2d (char*)stream (char*)setname (double)used_nstream (double)used_nset (char*)tag (char*)bcs (char*)brname [(int)fl_appRun]" << std::endl
					<< "[  stream  ] 0,1,2,3,4,5,01,0-5" << std::endl
					<< "[  setname ] A,B,..,U, AB, A-U"  << std::endl
					<< std::endl, abort();
  Char_t*  stream       = argv[1];
  Char_t*  setname      = argv[2];
  Double_t used_nstream = atof(argv[3]); 
  Double_t used_nset    = atof(argv[4]);
  Char_t*  tag          = argv[5];
  Char_t*  bcs          = argv[6];
  Char_t*  brname       = argv[7];
  Int_t    fl_appRun    = 1;
  if( argc==9 ) fl_appRun = atoi( argv[8] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Int_t fl_lrnb = 0; // 0(LR), 1(NB)
  if(      strncmp(brname, "lr", 2) == 0 ) fl_lrnb=0; // LR
  else if( strncmp(brname, "nb", 2) == 0 ) fl_lrnb=1; // NB
  else std::cerr << "[ABORT] Wrong branch-name : " << brname << std::endl, abort();
  if( fl_lrnb==0 ) std::cerr << "[ABORT] " << argv[0] << " is not supported for LR " << std::endl, abort();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t    Nchain             = 28; // fixed !! [7 categories]x[2 lepton-type]x[2 Xs-region], [ true, pre-true, false, uds, charm, mixed, charged ]-> [low-xs(ee)] : [high-xs(ee)] : [low-xs(mm)] : [high-xs(mm)]
  const Int_t    Nhist              = 30; // fixed !! [5 categories]x[2 lepton-type]x[3 Xs-region], [ true, pre-true, false, qq, bb                     ]-> [low-xs(ee)] : [high-xs(ee)] : [low-xs(mm)] : [high-xs(mm)]
  const Int_t    nfile[Nchain]      = {0};
  const Int_t    fl_mode_ll[Nchain] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,   // 1(ee)
					0,0,0,0,0,0,0,0,0,0,0,0,0,0 }; // 0(mm)
  const Char_t*  indir_sig[2]       = { "", Form("NB_lep_calib/hbk/hbk_%s_%s_522_bcs_merge/",tag,bcs) };
  const Char_t*  indir_bkg[2]       = { "", Form("NB_lep_calib/hbk/hbk_%s_%s_522_bcs/",      tag,bcs) };
  const Char_t*  tail               = "*.root";
  const Double_t scale_event_sig    = sigmc_amount*(used_nset/(Double_t)nset); // sigmc : N -> N/alpha
  const Double_t scale_event_bkg    = used_nstream;                            //   gmc : N -> N/alpha 
  const Bool_t   fl_message         = !true;
  const Bool_t   flag_k4pi          = !true; // 1(veto  K4pi    modes)
  const Bool_t   flag_unflavor      = !true; // 1(veto unflavor modes)
  const Bool_t   flag_ccpi0         = true; // 1(veto ccpi0 peak    )
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if(      i%7 < 2.5 ) sTmp << indir_sig[fl_lrnb] << "/sigMC_*_set[" << setname << "]"; // sig
    else if( i%7 > 2.5 ) sTmp << indir_bkg[fl_lrnb] << "/gMC_" << bkgtype[i%7-3] << "_*_s0["    << stream  << "]"; // bkg
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // sig               (low-xs,ee)
    {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // false with(t,t)   (low-xs,ee)
    {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // false except(t,t) (low-xs,ee)
    {0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // qq                (low-xs,ee)
    {0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // bb                (low-xs,ee)
    {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // sig               (high-xs,ee)
    {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // false with(t,t)   (high-xs,ee)
    {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // false except(t,t) (high-xs,ee)
    {0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // qq                (high-xs,ee)
    {0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // bb                (high-xs,ee)
    {1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // sig               (total-xs,ee)
    {0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // false with(t,t)   (total-xs,ee)
    {0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // false except(t,t) (total-xs,ee)
    {0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // qq                (total-xs,ee)
    {0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // bb                (total-xs,ee)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0}, // sig               (low-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}, // false with(t,t)   (low-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0}, // false except(t,t) (low-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0}, // qq                (low-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0}, // bb                (low-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0}, // sig               (high-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}, // false with(t,t)   (high-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0}, // false except(t,t) (high-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0}, // qq                (high-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}, // bb                (high-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0}, // sig               (total-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0}, // false with(t,t)   (total-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0}, // false except(t,t) (total-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0}, // qq                (total-xs,mm)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1}, // bb                (total-xs,mm)
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( int i=0; i<Nchain; i++) add_cut[i] = new Char_t[1024];

  for( int i=0; i<Nchain; i++ ) {
    if(      i%7==0 ) sTmp << "self==1";                              // true
    else if( i%7==1 ) sTmp << "self!=1 && " << makeCut_q2fl( 1,1   ); // false  with  q2=true, fl=true
    else if( i%7==2 ) sTmp << "self!=1 && " << makeCut_q2fl( 1,1,1 ); // false except q2=true, fl=true
    else              sTmp << "1"; // gmc
    
    if(      (i/7)%2==0 ) sTmp << " && xs_m < 1.1"; // low-Xs
    else if( (i/7)%2==1 ) sTmp << " && xs_m > 1.1"; // high-Xs
    
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    if( flag_ccpi0    ) sTmp << " && " << cut_ccpi0;
    strcpy( add_cut[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  
  const Char_t*  tname     = "h511";
  const Char_t*  axis      = brname;
  const Double_t yoffset   =  0.0;
  const Int_t    ybin      =   50;
  //const Int_t    ybin      =  100;
  const Double_t ymin  [2] = { 0.0,       -1.0    }; // {LR, NB}
  const Double_t ymax  [2] = { 1.0,        1.0    }; // {LR, NB}
  const Char_t*  ylabel[2] = { "LR(qq)", "NB(qq)" }; // {LR, NB}
  const Double_t xoffset   =  0.0;
  const Int_t    xbin      =   50;
  //const Int_t    xbin      =  100;
  const Double_t xmin  [2] = { 0.0,       -1.0     }; // {LR, NB}
  const Double_t xmax  [2] = { 1.0,        1.0     }; // {LR, NB}
  const Char_t*  xlabel[2] = { "LR(bb)",  "NB(bb)" }; // {LR, NB}

  
  const Bool_t flag_save  = true; // outfile.eps and outfile.root
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH2D**    tmphist = new TH2D*  [Nchain];
  TH2D**    hist    = new TH2D*  [Nhist];
  TH1D**    hist_X  = new TH1D*  [Nhist];
  TH1D**    hist_Y  = new TH1D*  [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",Nhist ); // LR or NB
  TCanvas*  c2      = Canvas( "c2","c2",3, 4  ); // LR or NB (projection)
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s, %s ) *************************************",tname,axis,bcs) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll[j] )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
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
    tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin,xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin,yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), Form("%s_qq:%s_bb", Form(axis,fl_mode_ll[j]), Form(axis,fl_mode_ll[j])), add_cut[j] );
    //////////////////////////////// correction events with LR=1 or NB=1
    if( chain[j]->GetTree()->GetEntries( Form("(%s)&&(%s_qq>1||%s_bb>1)", add_cut[j],Form(axis,fl_mode_ll[j]),Form(axis,fl_mode_ll[j])) ) ){
      std::cerr << "[ABORT] Strange Evnets having LR>1 or NB>1  Exist !" << std::endl;
      std::cerr << Form( "(%s)&&(%s_qq>1||%s_bb>1)", add_cut[j],Form(axis,fl_mode_ll[j]),Form(axis,fl_mode_ll[j]) ) << std::endl;
      chain[j]->GetTree()->SetScanField(0);
      chain[j]->GetTree()->Scan( "*", Form("(%s)&&(%s>1)", add_cut[j],Form(axis,fl_mode_ll[j]) ) ); 
      abort();
    }
    Double_t tmp_entry = tmphist[j]->GetEntries();
    for( Int_t nx=0; nx<xbin; nx++ ){
      Double_t tmp_bin  = tmphist[j]->GetBinContent(nx+1, ybin   );
      Double_t tmp_over = tmphist[j]->GetBinContent(nx+1, ybin+1 );
      tmphist[j]->SetBinContent( nx+1, ybin,   tmp_bin+tmp_over );
      tmphist[j]->SetBinContent( nx+1, ybin+1, 0                );
    }
    for( Int_t ny=0; ny<ybin+1; ny++ ){
      Double_t tmp_bin  = tmphist[j]->GetBinContent(xbin,   ny+1 );
      Double_t tmp_over = tmphist[j]->GetBinContent(xbin+1, ny+1 );
      tmphist[j]->SetBinContent( xbin,   ny+1, tmp_bin+tmp_over );
      tmphist[j]->SetBinContent( xbin+1, ny+1, 0                );
    }
    tmphist[j]->SetEntries( tmp_entry );
    ////////////////////////////////
    tmphist[j]->Sumw2();
    if(      j%7 < 2.5 ) tmphist[j]->Scale( 1/scale_event_sig );
    else if( j%7 > 2.5 ) tmphist[j]->Scale( 1/scale_event_bkg );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH2D( Form("hist%d",i), Form("hist%d",i), xbin,xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin,yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    hist[i]->SetXTitle( xlabel[fl_lrnb] );
    hist[i]->SetYTitle( ylabel[fl_lrnb] );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }

  // +++++++ display ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = 0; 
    Double_t entry_over   = 0;
    for( Int_t n=0; n<xbin; n++ ){
      entry_under += hist[i]->GetBinContent(0,      n+1);
      entry_over  += hist[i]->GetBinContent(xbin+1, n+1);
    }
    for( Int_t n=0; n<ybin; n++ ){
      entry_under += hist[i]->GetBinContent(n+1, 0      );
      entry_over  += hist[i]->GetBinContent(n+1, ybin+1 );
    }
    entry_under += hist[i]->GetBinContent( 0,      0      );
    entry_over  += hist[i]->GetBinContent( xbin+1, ybin+1 );
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
    sTmp += "]";
      
    //hist[i]->SetTitle( sTmp.Data() );
    std::cout << sTmp.Data() << std::endl;
    hist[i]->Scale( 1/(hist[i]->Integral()) );
  }
  hist[ 0]->SetTitle( "true,     low-Xs,   ee" );
  hist[ 1]->SetTitle( "pre-true, low-Xs,   ee" );
  hist[ 2]->SetTitle( "false,    low-Xs,   ee" );
  hist[ 3]->SetTitle( "qq,       low-Xs,   ee" );
  hist[ 4]->SetTitle( "bb,       low-Xs,   ee" );
  hist[ 5]->SetTitle( "true,     high-Xs,  ee" );
  hist[ 6]->SetTitle( "pre-true, high-Xs,  ee" );
  hist[ 7]->SetTitle( "false,    high-Xs,  ee" );
  hist[ 8]->SetTitle( "qq,       high-Xs,  ee" );
  hist[ 9]->SetTitle( "bb,       high-Xs,  ee" );
  hist[10]->SetTitle( "true,     total-Xs, ee" );
  hist[11]->SetTitle( "pre-true, total-Xs, ee" );
  hist[12]->SetTitle( "false,    total-Xs, ee" );
  hist[13]->SetTitle( "qq,       total-Xs, ee" );
  hist[14]->SetTitle( "bb,       total-Xs, ee" );

  hist[15]->SetTitle( "true,     low-Xs,   mm" );
  hist[16]->SetTitle( "pre-true, low-Xs,   mm" );
  hist[17]->SetTitle( "false,    low-Xs,   mm" );
  hist[18]->SetTitle( "qq,       low-Xs,   mm" );
  hist[19]->SetTitle( "bb,       low-Xs,   mm" );
  hist[20]->SetTitle( "true,     high-Xs,  mm" );
  hist[21]->SetTitle( "pre-true, high-Xs,  mm" );
  hist[22]->SetTitle( "false,    high-Xs,  mm" );
  hist[23]->SetTitle( "qq,       high-Xs,  mm" );
  hist[24]->SetTitle( "bb,       high-Xs,  mm" );
  hist[25]->SetTitle( "true,     total-Xs, mm" );
  hist[26]->SetTitle( "pre-true, total-Xs, mm" );
  hist[27]->SetTitle( "false,    total-Xs, mm" );
  hist[28]->SetTitle( "qq,       total-Xs, mm" );
  hist[29]->SetTitle( "bb,       total-Xs, mm" );

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  for( Int_t i=0; i<Nhist; i++ ){
    c1->cd(i+1);
    hist[i]->Draw( "COLZ" );
    //hist[i]->Draw( "SURF2Z" );
  }


  c2->Draw();
  for( Int_t i=0; i<Nhist; i++ ){
    hist_X[i] = hist[i]->ProjectionX();
    hist_Y[i] = hist[i]->ProjectionY();
    //Deco( hist_X[i], 2, i%5+1, i%5+1 );
    //Deco( hist_Y[i], 2, i%5+1, i%5+1 );
    if( i%5==0 ){
      Deco( hist_X[i], 2, 2, 2 );
      Deco( hist_Y[i], 2, 2, 2 );
    }else if( i%5==1 ){
      Deco( hist_X[i], 2, 5, 5 );
      Deco( hist_Y[i], 2, 5, 5 );
    }else if( i%5==2 ){
      Deco( hist_X[i], 2, 6, 6 );
      Deco( hist_Y[i], 2, 6, 6 );
    }else if( i%5==3 ){
      Deco( hist_X[i], 2, 3, 3 );
      Deco( hist_Y[i], 2, 3, 3 );
    }else if( i%5==4 ){
      Deco( hist_X[i], 2, 4, 4 );
      Deco( hist_Y[i], 2, 4, 4 );
    }

  }

  TH2D** waku_X = new TH2D*[6]; // [low-xs, high-xs, total-xs] x [ee, mm]
  TH2D** waku_Y = new TH2D*[6]; // [low-xs, high-xs, total-xs] x [ee, mm]
  waku_X[0] =  Waku( Nhist, hist_X, Form( "%s  (Low-Xs,ee)", xlabel[fl_lrnb]) );
  waku_Y[0] =  Waku( Nhist, hist_Y, Form( "%s  (Low-Xs,ee)", ylabel[fl_lrnb]) );
  waku_X[1] =  Waku( Nhist, hist_X, Form( "%s (High-Xs,ee)", xlabel[fl_lrnb]) );
  waku_Y[1] =  Waku( Nhist, hist_Y, Form( "%s (High-Xs,ee)", ylabel[fl_lrnb]) );
  waku_X[2] =  Waku( Nhist, hist_X, Form( "%s (Total-Xs,ee)",xlabel[fl_lrnb]) );
  waku_Y[2] =  Waku( Nhist, hist_Y, Form( "%s (Total-Xs,ee)",ylabel[fl_lrnb]) );

  waku_X[3] =  Waku( Nhist, hist_X, Form( "%s  (Low-Xs,mm)", xlabel[fl_lrnb]) );
  waku_Y[3] =  Waku( Nhist, hist_Y, Form( "%s  (Low-Xs,mm)", ylabel[fl_lrnb]) );
  waku_X[4] =  Waku( Nhist, hist_X, Form( "%s (High-Xs,mm)", xlabel[fl_lrnb]) );
  waku_Y[4] =  Waku( Nhist, hist_Y, Form( "%s (High-Xs,mm)", ylabel[fl_lrnb]) );
  waku_X[5] =  Waku( Nhist, hist_X, Form( "%s (Total-Xs,mm)",xlabel[fl_lrnb]) );
  waku_Y[5] =  Waku( Nhist, hist_Y, Form( "%s (Total-Xs,mm)",ylabel[fl_lrnb]) );

  waku_X[0]->SetTitle( Form("  Low-Xs,ee (stream=%s,set=%s)", stream,setname) );
  waku_X[1]->SetTitle( Form(" High-Xs,ee (stream=%s,set=%s)", stream,setname) );
  waku_X[2]->SetTitle( Form("Total-Xs,ee (stream=%s,set=%s)", stream,setname) );
  waku_X[3]->SetTitle( Form("  Low-Xs,mm (stream=%s,set=%s)", stream,setname) );
  waku_X[4]->SetTitle( Form(" High-Xs,mm (stream=%s,set=%s)", stream,setname) );
  waku_X[5]->SetTitle( Form("Total-Xs,mm (stream=%s,set=%s)", stream,setname) );
  waku_Y[0]->SetTitle( Form("  Low-Xs,ee (stream=%s,set=%s)", stream,setname) );
  waku_Y[1]->SetTitle( Form(" High-Xs,ee (stream=%s,set=%s)", stream,setname) );
  waku_Y[2]->SetTitle( Form("Total-Xs,ee (stream=%s,set=%s)", stream,setname) );
  waku_Y[3]->SetTitle( Form("  Low-Xs,mm (stream=%s,set=%s)", stream,setname) );
  waku_Y[4]->SetTitle( Form(" High-Xs,mm (stream=%s,set=%s)", stream,setname) );
  waku_Y[5]->SetTitle( Form("Total-Xs,mm (stream=%s,set=%s)", stream,setname) );

  for( Int_t i=0; i<6; i++ ){
    c2->cd(i+1); waku_X[i]->Draw();
    c2->cd(i+7); waku_Y[i]->Draw();
  }
  
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist_X[0],"Sig.(true)",     "PL" );
  legend1->AddEntry( hist_X[1],"Sig.(pre-true)", "PL" );
  legend1->AddEntry( hist_X[2],"Sig.(false)",    "PL" );
  legend1->AddEntry( hist_X[3],"Bkg.(qq)",       "PL" );
  legend1->AddEntry( hist_X[4],"Bkg.(bb)",       "PL" );
  legend1->Draw();

  

  for( Int_t i=0; i<Nhist; i++ ){
    c2->cd(i/5+1);
    if( i%5==0 || i%5==3 || i%5==4 ){
      if( i%5==0 || i%5==4 ) hist_X[i]->SetMarkerStyle(24), hist_X[i]->SetMarkerSize(0.8);
      hist_X[i]->Draw("same");
    }
    c2->cd(i/5+7);
    if( i%5==0 || i%5==3 || i%5==4 ){
      if( i%5==0 || i%5==3 ) hist_Y[i]->SetMarkerStyle(24), hist_Y[i]->SetMarkerSize(0.8);
      hist_Y[i]->Draw("same");
    }
  }
  
  c1->Update();
  c2->Update();
  if( flag_save ){
    c1->Print(     Form("pic/2d_%s_lep_%s_%s_s0%s_set%s_proj_c1.eps", lrnb[fl_lrnb],tag,bcs,stream,setname)             );
    c2->Print(     Form("pic/2d_%s_lep_%s_%s_s0%s_set%s_proj_c2.eps", lrnb[fl_lrnb],tag,bcs,stream,setname)             );
    TFile outfile( Form("pic/2d_%s_lep_%s_%s_s0%s_set%s_proj.root",   lrnb[fl_lrnb],tag,bcs,stream,setname), "RECREATE" );
    c1->Write();
    c2->Write();
    for( Int_t i=0; i<Nhist; i++ ) hist[i]->Write(), hist_X[i]->Write(), hist_Y[i]->Write();
    outfile.Close();
  }

  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] hist_X;
  delete[] hist_Y;
  delete[] add_cut;
  delete legend1;
  delete c1;

  return 0;
}

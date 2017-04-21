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
  const Int_t    Nchain             = 16; // fixed !! [4 category]x[2 lepton-type]x[2 Xs-region], [ sigMC, sigMC(t,t), sigMC(except-t,t), gMC ]-> [low-xs(ee)] : [high-xs(ee)] : [low-xs(mm)] : [high-xs(mm)]
  const Int_t    Nhist              =  8; // fixed !! [ sig/bkg  ]x[2 lepton-type]x[2 Xs-region], [ sig, bkg ] -> [ee,low-xs] : [ee,high-xs] : [mm,low-xs] : [mm,high-xs]
  const Int_t    nfile[Nchain]      = {0};
  const Int_t    fl_mode_ll[Nchain] = { 1,1,1,1,1,1,1,1,   // 1(ee)
					0,0,0,0,0,0,0,0 }; // 0(mm)
  const Char_t*  indir_sig[2]       = { "", Form("NB_lep_calib/hbk/hbk_%s_%s_522_bcs_merge/",tag,bcs) };
  const Char_t*  indir_bkg[2]       = { "", Form("NB_lep_calib/hbk/hbk_%s_%s_522_bcs/",      tag,bcs) };
  const Char_t*  tail               = "*.root";
  const Double_t scale_event_sig    = sigmc_amount*(used_nset/(Double_t)nset); // sigmc : N -> N/alpha
  const Double_t scale_event_bkg    = used_nstream;                            //   gmc : N -> N/alpha 
  const Bool_t   fl_message         = !true;
  const Bool_t   flag_k4pi          = !true; // 1(veto  K4pi    modes)
  const Bool_t   flag_unflavor      = !true; // 1(veto unflavor modes)
  const Bool_t   flag_ccpi0         = !true; // 1(veto ccpi0 peak    )
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if(      i%4 < 2.5 ) sTmp << indir_sig[fl_lrnb] << "/sigMC_*_set[" << setname << "]"; // sig
    else if( i%4 > 2.5 ) sTmp << indir_bkg[fl_lrnb] << "/gMC_*_s0["    << stream  << "]"; // bkg
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1,1},     // sig(ee, low-xs)
    {0,0,1,1}, // bkg(ee, low-xs)
    {0,0,0,0,1,1},     // sig(ee, high-xs)
    {0,0,0,0,0,0,1,1}, // bkg(ee, high-xs)
    {0,0,0,0,0,0,0,0,1,1},     // sig(mm, low-xs)
    {0,0,0,0,0,0,0,0,0,0,1,1}, // bkg(mm, low-xs)
    {0,0,0,0,0,0,0,0,0,0,0,0,1,1},     // sig(mm, high-xs)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}, // bkg(mm, high-xs)
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( int i=0; i<Nchain; i++) add_cut[i] = new Char_t[1024];

  for( int i=0; i<Nchain; i++ ) {
    if(      i%4==0 ) sTmp << "self==1";                              // true
    else if( i%4==1 ) sTmp << "self!=1 && " << makeCut_q2fl( 1,1   ); // false  with  q2=true, fl=true
    else if( i%4==2 ) sTmp << "self!=1 && " << makeCut_q2fl( 1,1,1 ); // false except q2=true, fl=true
    else if( i%4==3 ) sTmp << "1";
    
    if(      (i/4)%2==0 ) sTmp << " && xs_m < 1.1"; // low-Xs
    else if( (i/4)%2==1 ) sTmp << " && xs_m > 1.1"; // high-Xs

    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    if( flag_ccpi0    ) sTmp << " && " << cut_ccpi0;
    
    //sTmp << " && !(lpgt==3 && lmgt==3 && (dzll3dcalib>0.0190 || kfbclcalib<1.0e-18) )";
    //sTmp << " && dzll3d<0.0150"; // for test
    
    strcpy( add_cut[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  
  const Char_t*  tname     = "h511";
  const Char_t*  axis      = brname;
  const Double_t yoffset   =  0.0;
  const Int_t    ybin  [2] = {  50,        100    }; // {LR, NB}
  //const Int_t    ybin  [2] = { 100,        200    }; // {LR, NB}
  const Double_t ymin  [2] = { 0.0,       -1.0    }; // {LR, NB}
  const Double_t ymax  [2] = { 1.0,        1.0    }; // {LR, NB}
  const Char_t*  ylabel[2] = { "LR(qq)", "NB(qq)" }; // {LR, NB}
  const Double_t xoffset   =  0.0;
  const Int_t    xbin  [2] = {  50,        100     }; // {LR, NB}
  //const Int_t    xbin  [2] = { 100,        200     }; // {LR, NB}
  const Double_t xmin  [2] = { 0.0,       -1.0     }; // {LR, NB}
  const Double_t xmax  [2] = { 1.0,        1.0     }; // {LR, NB}
  const Char_t*  xlabel[2] = { "LR(bb)",  "NB(bb)" }; // {LR, NB}

  const Bool_t flag_save  = true; // outfile.eps and outfile.root
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH2D**    tmphist = new TH2D*  [Nchain];
  TH2D**    hist    = new TH2D*  [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",8 ); // LR or NB
  TCanvas*  c2      = Canvas( "c2","c2",4 ); // sig-eff v.s. bkg-rej-eff
  TCanvas*  c3      = Canvas( "c3","c3",5 ); // significance
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
    tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
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
    for( Int_t nx=0; nx<xbin[fl_lrnb]; nx++ ){
      Double_t tmp_bin  = tmphist[j]->GetBinContent(nx+1, ybin[fl_lrnb]   );
      Double_t tmp_over = tmphist[j]->GetBinContent(nx+1, ybin[fl_lrnb]+1 );
      tmphist[j]->SetBinContent( nx+1, ybin[fl_lrnb],   tmp_bin+tmp_over );
      tmphist[j]->SetBinContent( nx+1, ybin[fl_lrnb]+1, 0                );
    }
    for( Int_t ny=0; ny<ybin[fl_lrnb]+1; ny++ ){
      Double_t tmp_bin  = tmphist[j]->GetBinContent(xbin[fl_lrnb],   ny+1 );
      Double_t tmp_over = tmphist[j]->GetBinContent(xbin[fl_lrnb]+1, ny+1 );
      tmphist[j]->SetBinContent( xbin[fl_lrnb],   ny+1, tmp_bin+tmp_over );
      tmphist[j]->SetBinContent( xbin[fl_lrnb]+1, ny+1, 0                );
    }
    tmphist[j]->SetEntries( tmp_entry );
    ////////////////////////////////
    tmphist[j]->Sumw2();
    if(      j%4 < 2.5 ) tmphist[j]->Scale( 1/scale_event_sig );
    else if( j%4 > 2.5 ) tmphist[j]->Scale( 1/scale_event_bkg );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    if(      i==0 ) hist[i] = new TH2D( Form("hist_%s_lep%d_%s_%s_lowxs_sig", lrnb[fl_lrnb],1,tag,bcs), Form("hist_%s_lep%d_%s_%s_lowxs_sig",  lrnb[fl_lrnb],1,tag,bcs), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
    else if( i==1 ) hist[i] = new TH2D( Form("hist_%s_lep%d_%s_%s_lowxs_bkg", lrnb[fl_lrnb],1,tag,bcs), Form("hist_%s_lep%d_%s_%s_lowxs_bkg",  lrnb[fl_lrnb],1,tag,bcs), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
    else if( i==2 ) hist[i] = new TH2D( Form("hist_%s_lep%d_%s_%s_highxs_sig",lrnb[fl_lrnb],1,tag,bcs), Form("hist_%s_lep%d_%s_%s_highxs_sig", lrnb[fl_lrnb],1,tag,bcs), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
    else if( i==3 ) hist[i] = new TH2D( Form("hist_%s_lep%d_%s_%s_highxs_bkg",lrnb[fl_lrnb],1,tag,bcs), Form("hist_%s_lep%d_%s_%s_highxs_bkg", lrnb[fl_lrnb],1,tag,bcs), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
    else if( i==4 ) hist[i] = new TH2D( Form("hist_%s_lep%d_%s_%s_lowxs_sig", lrnb[fl_lrnb],0,tag,bcs), Form("hist_%s_lep%d_%s_%s_lowxs_sig",  lrnb[fl_lrnb],0,tag,bcs), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
    else if( i==5 ) hist[i] = new TH2D( Form("hist_%s_lep%d_%s_%s_lowxs_bkg", lrnb[fl_lrnb],0,tag,bcs), Form("hist_%s_lep%d_%s_%s_lowxs_bkg",  lrnb[fl_lrnb],0,tag,bcs), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
    else if( i==6 ) hist[i] = new TH2D( Form("hist_%s_lep%d_%s_%s_highxs_sig",lrnb[fl_lrnb],0,tag,bcs), Form("hist_%s_lep%d_%s_%s_highxs_sig", lrnb[fl_lrnb],0,tag,bcs), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
    else if( i==7 ) hist[i] = new TH2D( Form("hist_%s_lep%d_%s_%s_highxs_bkg",lrnb[fl_lrnb],0,tag,bcs), Form("hist_%s_lep%d_%s_%s_highxs_bkg", lrnb[fl_lrnb],0,tag,bcs), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );

    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    hist[i]->SetXTitle( xlabel[fl_lrnb] );
    hist[i]->SetYTitle( ylabel[fl_lrnb] );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }
  // +++++++ calculation +++++++++++++++++++++
  TGraph** graph_eff = new TGraph*[Nhist/2]; // low-xs, high-xs
  TH2D**   hist_sign = new TH2D*  [Nhist/2]; // low-xs, high-xs
  hist_sign[0] = new TH2D( Form("sign_lowxs_lep%d",  1), Form("sign_lowxs_lep%d",  1), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
  hist_sign[1] = new TH2D( Form("sign_highxs_lep%d", 1), Form("sign_highxs_lep%d", 1), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
  hist_sign[2] = new TH2D( Form("sign_lowxs_lep%d",  0), Form("sign_lowxs_lep%d",  0), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
  hist_sign[3] = new TH2D( Form("sign_highxs_lep%d", 0), Form("sign_highxs_lep%d", 0), xbin[fl_lrnb],xoffset+xmin[fl_lrnb],xoffset+xmax[fl_lrnb], ybin[fl_lrnb],yoffset+ymin[fl_lrnb],yoffset+ymax[fl_lrnb] );
  for( Int_t i=0; i<Nhist/2; i++ ){
    graph_eff[i]  = new TGraph();
    Deco( graph_eff[i],  0, i+1, i+1 );
    hist_sign[i]->GetXaxis()->CenterTitle();
    hist_sign[i]->GetYaxis()->CenterTitle();
    hist_sign[i]->SetXTitle( xlabel[fl_lrnb] );
    hist_sign[i]->SetYTitle( ylabel[fl_lrnb] );
    hist_sign[i]->SetAxisRange( 0,10,"Z" );
  }
  graph_eff[0] ->SetName( "eff_lowxs_lep1"  );
  graph_eff[1] ->SetName( "eff_highxs_lep1"  );
  graph_eff[2] ->SetName( "eff_lowxs_lep0" );
  graph_eff[3] ->SetName( "eff_highxs_lep0"  );

  Double_t entry_sig[Nhist/2] = { hist[0]->Integral(), hist[2]->Integral(), hist[4]->Integral(), hist[6]->Integral() }; // { low-xs(ee), high-xs(ee), low-xs(mm), high-xs(mm) }
  Double_t entry_bkg[Nhist/2] = { hist[1]->Integral(), hist[3]->Integral(), hist[5]->Integral(), hist[7]->Integral() }; // { low-xs(ee), high-xs(ee), low-xs(mm), high-xs(mm) }

  std::cout << "=================================" << std::endl
	    << "[Scale Factor]"                    << std::endl
	    << "sig : " << scale_event_sig
	    << "("      << used_nset << " set)"    << std::endl
	    << "bkg : " << scale_event_bkg         << std::endl
	    << "=================================" << std::endl;
  if( fl_message ){
    std::cout << "===================================================" << std::endl
	      << "nsig(  low-xs, ee ) = " << entry_sig[0] << std::endl
	      << "nsig( high-xs, ee ) = " << entry_sig[1] << std::endl
	      << "nsig(  low-xs, mm ) = " << entry_sig[2] << std::endl
	      << "nsig( high-xs, mm ) = " << entry_sig[3] << std::endl
	      << "nbkg(  low-xs, ee ) = " << entry_bkg[0] << std::endl
	      << "nbkg( high-xs, ee ) = " << entry_bkg[1] << std::endl
	      << "nbkg(  low-xs, mm ) = " << entry_bkg[2] << std::endl
	      << "nbkg( high-xs, mm ) = " << entry_bkg[3] << std::endl
	      << "===================================================" << std::endl;
  }


  Double_t max_sign [Nhist/2] = {0}; // maximum significance             [low-xs, high-xs]
  Double_t max_sig  [Nhist/2] = {0}; // # of sig at maximum significance [low-xs, high-xs]
  Double_t max_bkg  [Nhist/2] = {0}; // # of bkg at maximum significance [low-xs, high-xs]
  Double_t max_lr_qq[Nhist/2] = {0}; // LR value at maximum significance [low-xs, high-xs]
  Double_t max_lr_bb[Nhist/2] = {0}; // LR value at maximum significance [low-xs, high-xs]
  Int_t    tmp_ind=0;
  for( Int_t kx=0; kx<xbin[fl_lrnb]; kx++ ){
    for( Int_t ky=0; ky<ybin[fl_lrnb]; ky++ ){
      for( Int_t i=0; i<Nhist/2; i++ ){
	Double_t nsig_err   = 0;
	Double_t nbkg_err   = 0;
	Double_t nsig       = hist[2*i  ]->IntegralAndError( kx+1,xbin[fl_lrnb], ky+1,ybin[fl_lrnb], nsig_err );
	Double_t nbkg       = hist[2*i+1]->IntegralAndError( kx+1,xbin[fl_lrnb], ky+1,ybin[fl_lrnb], nbkg_err );
	Double_t sign       = nsig==0 ? 0 : nsig / sqrt(nsig+nbkg);
	Double_t round_sig  = nsig==0 ? 0 : sign/nsig - sign*sign*sign/nsig/nsig/2.0;
	Double_t round_bkg  = nsig==0 ? 0 : -sign*sign*sign/nsig/nsig/2.0;
	Double_t sign_error = sqrt( round_sig*round_sig*nsig_err*nsig_err + round_bkg*round_bkg*nbkg_err*nbkg_err );
	graph_eff[i]->SetPoint     ( tmp_ind++, nsig/entry_sig[i], 1.0-nbkg/entry_bkg[i]    );
	hist_sign[i]->SetBinContent( kx+1, ky+1, sign );
	if( sign > max_sign[i] ){
	  max_sign [i] = sign;
	  max_sig  [i] = nsig;
	  max_bkg  [i] = nbkg;
	  max_lr_qq[i] = hist[2*i]->GetBinLowEdge(ky+1);
	  max_lr_bb[i] = hist[2*i]->GetBinLowEdge(kx+1);
	}
	if( fl_message ){
	  Double_t tmp_x, tmp_y;
	  std::cout << "Bin = (" << kx+1 << ", " << ky+1 << "), i = " << i << std::endl;
	  std::cout << " S = " << std::setw(7) << std::right << nsig << ","
		    << " B = " << std::setw(7) << std::right << nbkg << std::endl;
	  graph_eff[i]->GetPoint     ( kx+ky, tmp_x, tmp_y ); 
	  std::cout << "[ eff(sig) = " << std::setw(8) << std::right << tmp_x
		    << ", eff(bkg) = " << std::setw(8) << std::right << tmp_y << " ] ";
	  std::cout << "[ LH(qq) = "   << std::setw(4) << std::right << hist_sign[i]->GetBinLowEdge(ky+1)
		    << ", LH(bb) = "   << std::setw(4) << std::right << hist_sign[i]->GetBinLowEdge(kx+1)
		    << " -> sign = "   << std::setw(8) << std::right << hist_sign[i]->GetBinContent( kx+1, ky+1 ) << "] " << std::endl;
	}
      }
    }
  }
  // +++++++ display ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = 0; 
    Double_t entry_over   = 0;
    for( Int_t n=0; n<xbin[fl_lrnb]; n++ ){
      entry_under += hist[i]->GetBinContent(0,      n+1);
      entry_over  += hist[i]->GetBinContent(xbin[fl_lrnb]+1, n+1);
    }

    for( Int_t n=0; n<ybin[fl_lrnb]; n++ ){
      entry_under += hist[i]->GetBinContent(n+1, 0      );
      entry_over  += hist[i]->GetBinContent(n+1, ybin[fl_lrnb]+1 );
    }
    entry_under += hist[i]->GetBinContent( 0,      0      );
    entry_over  += hist[i]->GetBinContent( xbin[fl_lrnb]+1, ybin[fl_lrnb]+1 );
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
    sTmp += ") ]";
      
    hist[i]->SetTitle( sTmp.Data() );
    std::cout << sTmp.Data() << std::endl;
    hist[i]->Scale( 1/(hist[i]->Integral()) );
  }

  std::cout << "[ Maximum Significane ]"  << std::endl;
  for( Int_t i=0; i<Nhist/2; i++ ){
    if(      i==0 ) std::cout << " (  low-xs, ee ) ";
    else if( i==1 ) std::cout << " ( high-xs, ee ) ";
    else if( i==2 ) std::cout << " (  low-xs, mm ) ";
    else if( i==3 ) std::cout << " ( high-xs, mm ) ";
    std::cout << std::setw(10) << std::right << max_sign[i] << " ( S = "
	      << std::setw( 7) << std::right << max_sig[i]                 << ", B = "
	      << std::setw( 7) << std::right << max_bkg[i]                 << ", LR(qq) = "
	      << std::setw( 5) << std::right << max_lr_qq[i] << ", LR(bb) = "
	      << std::setw( 5) << std::right << max_lr_bb[i] << " )"
	      << std::endl;
  }

  std::cout << " (   total, ee ) "
	    << std::setw(10) << std::right << (max_sig[0] + max_sig[1])/sqrt(max_sig[0]+max_sig[1]+max_bkg[0]+max_bkg[1]) << " ( S = "
	    << std::setw( 7) << std::right <<  max_sig[0] + max_sig[1] << ", B = "
	    << std::setw( 7) << std::right <<  max_bkg[0] + max_bkg[1] << " )"
	    << std::endl;
  std::cout << " (   total, mm ) "
	    << std::setw(10) << std::right << (max_sig[2] + max_sig[3])/sqrt(max_sig[2]+max_sig[3]+max_bkg[2]+max_bkg[3]) << " ( S = "
	    << std::setw( 7) << std::right <<  max_sig[2] + max_sig[3] << ", B = "
	    << std::setw( 7) << std::right <<  max_bkg[2] + max_bkg[3] << " )"
	    << std::endl;
  std::cout << " (   total,  low-xs ) "
	    << std::setw(10) << std::right << (max_sig[0] + max_sig[2])/sqrt(max_sig[0]+max_sig[2]+max_bkg[0]+max_bkg[2]) << " ( S = "
	    << std::setw( 7) << std::right <<  max_sig[0] + max_sig[2] << ", B = "
	    << std::setw( 7) << std::right <<  max_bkg[0] + max_bkg[2] << " )"
	    << std::endl;
  std::cout << " (   total, high-xs ) "
	    << std::setw(10) << std::right << (max_sig[1] + max_sig[3])/sqrt(max_sig[1]+max_sig[3]+max_bkg[1]+max_bkg[3]) << " ( S = "
	    << std::setw( 7) << std::right <<  max_sig[1] + max_sig[3] << ", B = "
	    << std::setw( 7) << std::right <<  max_bkg[1] + max_bkg[3] << " )"
	    << std::endl;
  std::cout << " (   total          ) "
	    << std::setw(10) << std::right << (max_sig[0] + max_sig[1] + max_sig[2] + max_sig[3])/sqrt(max_sig[0]+max_sig[1]+max_sig[2]+max_sig[3]+max_bkg[0]+max_bkg[1]+max_bkg[2]+max_bkg[3]) << " ( S = "
	    << std::setw( 7) << std::right <<  max_sig[0] + max_sig[1] + max_sig[2] + max_sig[3] << ", B = "
	    << std::setw( 7) << std::right <<  max_bkg[0] + max_bkg[1] + max_bkg[2] + max_bkg[3] << " )"
	    << std::endl;
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  for( Int_t i=0; i<Nhist; i++ ){
    c1->cd(i+1);
    hist[i]->Draw( "COLZ" );
  }
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c2->Draw();

  TH2D* waku_eff  = new TH2D( "Performance of Likelihood Ratio Cut", "Performance of Likelihood Ratio Cut", 2,0,1,2,0,1  );
  waku_eff->GetXaxis()->CenterTitle();
  waku_eff->GetYaxis()->CenterTitle();
  waku_eff->SetXTitle( "Sig. Eff."           );
  waku_eff->SetYTitle( "Bkg. Rejection Eff." );
  for( Int_t i=0; i<Nhist/2; i++ ){
    c2->cd(i+1);
    waku_eff->Draw();
    graph_eff[i]->Draw("Psame");
  }
  // +++++++ tlegend2 ++++++++++++++++++++++++++++++++++
  TLegend* legend2 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend2->AddEntry( graph_eff[0],"sig(ee,Low-X_{s})",  "L" );
  legend2->AddEntry( graph_eff[1],"sig(ee,High-X_{s})", "L" );
  legend2->AddEntry( graph_eff[2],"sig(ee,Low-X_{s})",  "L" );
  legend2->AddEntry( graph_eff[3],"sig(ee,High-X_{s})", "L" );
  legend2->Draw();
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c3->Draw();
  for( Int_t i=0; i<Nhist/2; i++ ){
    c3->cd(i+1);
    hist_sign[i]->Draw("COLZ");
  }
  // +++++++ tlegend ++++++++++++++++++++++++++++++++++
  c3->cd(6);
  TPaveText* box    = new TPaveText( 0.0,0.4,1.0,1.0 );
  box->SetTextSize(0.015);
  TLegend*   legend = new TLegend  ( 0.0,0.0,1.0,0.4 );
  legend->SetTextSize(0.015);
  for( Int_t j=0; j<Nchain; j++ ){
    box->AddText( Form("tmphist%d : Mode(%s : %s : %d files)", j,chain[j]->GetModeName(),gSystem->BaseName(chain[j]->GetFileName()),chain[j]->GetNfile()) );
    box->AddText( Form("   %s", add_cut[j]) );
  }

  for( Int_t i=0; i<Nhist; i++ ){
    legend->AddEntry( hist[i],"","PL" );
  }
  box->Draw();
  legend->Draw();

  c1->Update();
  c2->Update();
  c3->Update();
  if( flag_save ){
    c1->Print(     Form("pic/2d_%s_lep_%s_%s_s0%s_set%s_c1.eps", lrnb[fl_lrnb],tag,bcs,stream,setname)             );
    c2->Print(     Form("pic/2d_%s_lep_%s_%s_s0%s_set%s_c2.eps", lrnb[fl_lrnb],tag,bcs,stream,setname)             );
    c3->Print(     Form("pic/2d_%s_lep_%s_%s_s0%s_set%s_c3.eps", lrnb[fl_lrnb],tag,bcs,stream,setname)             );
    TFile outfile( Form("pic/2d_%s_lep_%s_%s_s0%s_set%s.root",   lrnb[fl_lrnb],tag,bcs,stream,setname), "RECREATE" );
    c1->Write();
    c2->Write(); 
    c3->Write();
    for( Int_t i=0; i<Nhist;   i++ ) hist[i]->Write();
    for( Int_t i=0; i<Nhist/2; i++ ) graph_eff[i]->Write(), hist_sign[i]->Write();
    outfile.Close();
  }

  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] add_cut;
  delete legend;
  delete legend2;
  delete c1;
  delete c2;
  delete c3;

  return 0;
}

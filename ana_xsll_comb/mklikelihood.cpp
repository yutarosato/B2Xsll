#include "mklikelihood.h"

Int_t main( Int_t argc, Char_t** argv ){
  using namespace gmc;
  Style();
  if( argc!=4 ) std::cerr << "wrong input" << std::endl
					     << " Usage   : ./mklikelihood (char*)filename (char*)outdir (char*)tag" << std::endl
					     << " Example : ./mklikelihood sigMC_uds_..root        bkg/         1ks" << std::endl
					     << " Example : ./mklikelihood       plot               1           1  (plot the PDFs)"
					     << std::endl, abort();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Char_t* filename  = argv[1];
  const Char_t* treename  = "h511";
  const Char_t* outdir    = argv[2];
  const Char_t* basename  = gSystem->BaseName( filename );
  const Char_t* tag       = argv[3];
  const Char_t* pdf_dir    = "pdf_20120823/"; // PDF using setK-U and stream3-5 for ksfw0 (hbk5) kfbcl-cut remove deroe
  //const Char_t* pdf_dir    = "pdf_20120606/"; // PDF using setK-U and stream3-5 for ksfw0 (hbk3) kfbcl-cut (almost same with 20120505)  remove deroe
  //const Char_t* pdf_dir    = "pdf_20120529/"; // PDF using setK-U and stream3-5 for ksfw0 (hbk3) kfbcl-cut, remove-kid to lepton
  //const Char_t* pdf_dir    = "pdf_20120524/"; // PDF using setK-U and stream3-5 for ksfw1 (hbk3) kfbcl-cut, remove deroe
  //const Char_t* pdf_dir    = "pdf_20120518/"; // PDF using setK-U and stream3-5 for ksfw1 (hbk3) kfbcl-cut
  //const Char_t* pdf_dir    = "pdf_20120508/"; // PDF using setK-U and stream3-5 for ksfw0 (hbk3) kfbcl-cut
  //const Char_t* pdf_dir    = "pdf_20120416/"; // PDF using setK-U and stream3-5
  //const Char_t* pdf_dir    = "pdf_20120326/"; // PDF using setK-U and stream3-5
  //const Char_t* pdf_dir    = "pdf_20120403/"; // PDF ssing setK-U and stream0-2
  const Char_t* stream     = "3-5";
  //const Char_t* stream     = "0-2";
  const Char_t* setname    = "K-U";
  
  const Bool_t fl_message = !true;

  const Int_t fl_de      = 2; // 0 -> not use, 1 -> use hist-PDF, 2 -> use fit-PDF
  const Int_t fl_kfbchi2 = 2;
  const Int_t fl_dzll3d  = 2;
  const Int_t fl_bccm    = 2;
  const Int_t fl_fmiss   = 1; // fit-PDF is not supported yet!!!!
  const Int_t fl_ksfw    = 2;
  
  Int_t fl_mc=0;
  if( strncmp(basename, "sigMC", 5)==0 ) fl_mc=1; // sigmc
  if( strncmp(basename,   "gMC", 3)==0 ) fl_mc=2; // gmc
  if( strcmp (basename,  "plot"   )==0 ) fl_mc=3; // PDF plot
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // read PDF
  TFile file_de_lep1     ( Form("%s/pdf_de_lep1_s0%s_set%s.root",  pdf_dir, stream, setname  )  );
  TFile file_de_lep0     ( Form("%s/pdf_de_lep0_s0%s_set%s.root",  pdf_dir, stream, setname  )  );
  TFile file_dzll3d      ( Form("%s/pdf_dzll3d_s0%s_set%s.root",   pdf_dir, stream, setname  )  );
  TFile file_bccm        ( Form("%s/pdf_bccm_s0%s_set%s.root",     pdf_dir, stream, setname  )  );
  TFile file_sig_kfbchi2 ( Form("%s/pdf_kfbchi2_set%s_sig.root",   pdf_dir, setname          )  );
  TFile file_bkg_kfbchi2 ( Form("%s/pdf_kfbchi2_s0%s_bkg.root",    pdf_dir, stream           )  );
  TFile file_sig_fmiss   ( Form("%s/ksfw_fmiss%s_tot_deroe_sig.root",    pdf_dir, tag              )  ); // deroe
  TFile file_sig_fmiss_qq( Form("%s/ksfw_fmiss%s_qq_deroe_sig.root",     pdf_dir, tag              )  ); // deroe
  TFile file_sig_fmiss_bb( Form("%s/ksfw_fmiss%s_bb_deroe_sig.root",     pdf_dir, tag              )  ); // deroe
  TFile file_bkg_fmiss   ( Form("%s/ksfw_fmiss%s_tot_deroe_bkg.root",    pdf_dir, tag              )  ); // deroe
  TFile file_bkg_fmiss_qq( Form("%s/ksfw_fmiss%s_qq_deroe_bkg.root",     pdf_dir, tag              )  ); // deroe
  TFile file_bkg_fmiss_bb( Form("%s/ksfw_fmiss%s_bb_deroe_bkg.root",     pdf_dir, tag              )  ); // deroe
  //TFile file_sig_fmiss   ( Form("%s/ksfw_fmiss%s_tot_sig.root",    pdf_dir, tag              )  );
  //TFile file_sig_fmiss_qq( Form("%s/ksfw_fmiss%s_qq_sig.root",     pdf_dir, tag              )  );
  //TFile file_sig_fmiss_bb( Form("%s/ksfw_fmiss%s_bb_sig.root",     pdf_dir, tag              )  );
  //TFile file_bkg_fmiss   ( Form("%s/ksfw_fmiss%s_tot_bkg.root",    pdf_dir, tag              )  );
  //TFile file_bkg_fmiss_qq( Form("%s/ksfw_fmiss%s_qq_bkg.root",     pdf_dir, tag              )  );
  //TFile file_bkg_fmiss_bb( Form("%s/ksfw_fmiss%s_bb_bkg.root",     pdf_dir, tag              )  );

  if( file_de_lep1.IsZombie()      ) std::cerr << "[ABORT] can not find input-file : " << file_de_lep1.GetName()      << std::endl, abort();
  if( file_de_lep0.IsZombie()      ) std::cerr << "[ABORT] can not find input-file : " << file_de_lep0.GetName()      << std::endl, abort();
  if( file_dzll3d.IsZombie()       ) std::cerr << "[ABORT] can not find input-file : " << file_dzll3d.GetName()       << std::endl, abort();
  if( file_bccm.IsZombie()         ) std::cerr << "[ABORT] can not find input-file : " << file_bccm.GetName()         << std::endl, abort();
  if( file_sig_kfbchi2.IsZombie()  ) std::cerr << "[ABORT] can not find input-file : " << file_sig_kfbchi2.GetName()  << std::endl, abort();
  if( file_bkg_kfbchi2.IsZombie()  ) std::cerr << "[ABORT] can not find input-file : " << file_bkg_kfbchi2.GetName()  << std::endl, abort();
  if( file_sig_fmiss.IsZombie()    ) std::cerr << "[ABORT] can not find input-file : " << file_sig_fmiss.GetName()    << std::endl, abort();
  if( file_sig_fmiss_qq.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_sig_fmiss_qq.GetName() << std::endl, abort();
  if( file_sig_fmiss_bb.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_sig_fmiss_bb.GetName() << std::endl, abort();
  if( file_bkg_fmiss.IsZombie()    ) std::cerr << "[ABORT] can not find input-file : " << file_bkg_fmiss.GetName()    << std::endl, abort();
  if( file_bkg_fmiss_qq.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_bkg_fmiss_qq.GetName() << std::endl, abort();
  if( file_bkg_fmiss_bb.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_bkg_fmiss_bb.GetName() << std::endl, abort();
  
  TH1D* hist_sig_de_lep1         = (TH1D*)file_de_lep1.Get     ( "pdf_hist_de_lep1_sig"         );
  TH1D* hist_sig_de_lep1_pi0     = (TH1D*)file_de_lep1.Get     ( "pdf_hist_de_lep1_pi0_sig"     );
  TH1D* hist_sig_de_lep0         = (TH1D*)file_de_lep0.Get     ( "pdf_hist_de_lep0_sig"         );
  TH1D* hist_sig_de_lep0_pi0     = (TH1D*)file_de_lep0.Get     ( "pdf_hist_de_lep0_pi0_sig"     );
  TH1D* hist_bkg_de_lep1         = (TH1D*)file_de_lep1.Get     ( "pdf_hist_de_lep1_tot_bkg"     );
  TH1D* hist_bkg_de_lep1_qq      = (TH1D*)file_de_lep1.Get     ( "pdf_hist_de_lep1_qq_bkg"      );
  TH1D* hist_bkg_de_lep1_bb      = (TH1D*)file_de_lep1.Get     ( "pdf_hist_de_lep1_bb_bkg"      );
  TH1D* hist_bkg_de_lep0         = (TH1D*)file_de_lep0.Get     ( "pdf_hist_de_lep0_tot_bkg"     );
  TH1D* hist_bkg_de_lep0_qq      = (TH1D*)file_de_lep0.Get     ( "pdf_hist_de_lep0_qq_bkg"      );
  TH1D* hist_bkg_de_lep0_bb      = (TH1D*)file_de_lep0.Get     ( "pdf_hist_de_lep0_bb_bkg"      );
  TH1D* hist_sig_dzll3d          = (TH1D*)file_dzll3d.Get      ( "pdf_hist_dzll3d_sig"          );
  TH1D* hist_bkg_dzll3d          = (TH1D*)file_dzll3d.Get      ( "pdf_hist_dzll3d_tot_bkg"      );
  TH1D* hist_bkg_dzll3d_qq       = (TH1D*)file_dzll3d.Get      ( "pdf_hist_dzll3d_qq_bkg"       );
  TH1D* hist_bkg_dzll3d_bb       = (TH1D*)file_dzll3d.Get      ( "pdf_hist_dzll3d_bb_bkg"       );
  TH1D* hist_sig_bccm            = (TH1D*)file_bccm.Get        ( "pdf_hist_bccm_sig"            );
  TH1D* hist_bkg_bccm            = (TH1D*)file_bccm.Get        ( "pdf_hist_bccm_tot_bkg"        );
  TH1D* hist_bkg_bccm_qq         = (TH1D*)file_bccm.Get        ( "pdf_hist_bccm_qq_bkg"         );
  TH1D* hist_bkg_bccm_bb         = (TH1D*)file_bccm.Get        ( "pdf_hist_bccm_bb_bkg"         );
  TH1D* hist_sig_kfbchi2_0tracks = (TH1D*)file_sig_kfbchi2.Get ( "pdf_hist_kfbchi2_0tracks_sig" );
  TH1D* hist_sig_kfbchi2_1tracks = (TH1D*)file_sig_kfbchi2.Get ( "pdf_hist_kfbchi2_1tracks_sig" );
  TH1D* hist_sig_kfbchi2_2tracks = (TH1D*)file_sig_kfbchi2.Get ( "pdf_hist_kfbchi2_2tracks_sig" );
  TH1D* hist_sig_kfbchi2_3tracks = (TH1D*)file_sig_kfbchi2.Get ( "pdf_hist_kfbchi2_3tracks_sig" );
  TH1D* hist_sig_kfbchi2_4tracks = (TH1D*)file_sig_kfbchi2.Get ( "pdf_hist_kfbchi2_4tracks_sig" );
  TH1D* hist_bkg_kfbchi2         = (TH1D*)file_bkg_kfbchi2.Get ( "pdf_hist_kfbchi2_tot_bkg"     );
  TH1D* hist_bkg_kfbchi2_qq      = (TH1D*)file_bkg_kfbchi2.Get ( "pdf_hist_kfbchi2_qq_bkg"      );
  TH1D* hist_bkg_kfbchi2_bb      = (TH1D*)file_bkg_kfbchi2.Get ( "pdf_hist_kfbchi2_bb_bkg"      );
  TH1D* hist_sig_fmiss           = (TH1D*)file_sig_fmiss.Get   ( "ksfw0_ksfw0"                  );
  TH1D* hist_sig_fmiss_qq        = (TH1D*)file_sig_fmiss_qq.Get( "ksfw0_ksfw0"                  );
  TH1D* hist_sig_fmiss_bb        = (TH1D*)file_sig_fmiss_bb.Get( "ksfw0_ksfw0"                  );
  TH1D* hist_bkg_fmiss           = (TH1D*)file_bkg_fmiss.Get   ( "ksfw0_ksfw0"                  );
  TH1D* hist_bkg_fmiss_qq        = (TH1D*)file_bkg_fmiss_qq.Get( "ksfw0_ksfw0"                  );
  TH1D* hist_bkg_fmiss_bb        = (TH1D*)file_bkg_fmiss_bb.Get( "ksfw0_ksfw0"                  );
  
  if( hist_sig_de_lep1         == NULL ) std::cerr << "[ABORT] can not find histgram  1 : " << std::endl, abort();
  if( hist_sig_de_lep1_pi0     == NULL ) std::cerr << "[ABORT] can not find histgram  2 : " << std::endl, abort();
  if( hist_sig_de_lep0         == NULL ) std::cerr << "[ABORT] can not find histgram  3 : " << std::endl, abort();
  if( hist_sig_de_lep0_pi0     == NULL ) std::cerr << "[ABORT] can not find histgram  4 : " << std::endl, abort();
  if( hist_bkg_de_lep1         == NULL ) std::cerr << "[ABORT] can not find histgram  5 : " << std::endl, abort();
  if( hist_bkg_de_lep1_qq      == NULL ) std::cerr << "[ABORT] can not find histgram  6 : " << std::endl, abort();
  if( hist_bkg_de_lep1_bb      == NULL ) std::cerr << "[ABORT] can not find histgram  7 : " << std::endl, abort();
  if( hist_bkg_de_lep0         == NULL ) std::cerr << "[ABORT] can not find histgram  8 : " << std::endl, abort();
  if( hist_bkg_de_lep0_qq      == NULL ) std::cerr << "[ABORT] can not find histgram  9 : " << std::endl, abort();
  if( hist_bkg_de_lep0_bb      == NULL ) std::cerr << "[ABORT] can not find histgram 10 : " << std::endl, abort();
  if( hist_sig_dzll3d          == NULL ) std::cerr << "[ABORT] can not find histgram 11 : " << std::endl, abort();
  if( hist_bkg_dzll3d          == NULL ) std::cerr << "[ABORT] can not find histgram 12 : " << std::endl, abort();
  if( hist_bkg_dzll3d_qq       == NULL ) std::cerr << "[ABORT] can not find histgram 13 : " << std::endl, abort();
  if( hist_bkg_dzll3d_bb       == NULL ) std::cerr << "[ABORT] can not find histgram 14 : " << std::endl, abort();
  if( hist_sig_bccm            == NULL ) std::cerr << "[ABORT] can not find histgram 15 : " << std::endl, abort();
  if( hist_bkg_bccm            == NULL ) std::cerr << "[ABORT] can not find histgram 16 : " << std::endl, abort();
  if( hist_bkg_bccm_qq         == NULL ) std::cerr << "[ABORT] can not find histgram 17 : " << std::endl, abort();
  if( hist_bkg_bccm_bb         == NULL ) std::cerr << "[ABORT] can not find histgram 18 : " << std::endl, abort();
  if( hist_sig_kfbchi2_0tracks == NULL ) std::cerr << "[ABORT] can not find histgram 19 : " << std::endl, abort();
  if( hist_sig_kfbchi2_1tracks == NULL ) std::cerr << "[ABORT] can not find histgram 20 : " << std::endl, abort();
  if( hist_sig_kfbchi2_2tracks == NULL ) std::cerr << "[ABORT] can not find histgram 21 : " << std::endl, abort();
  if( hist_sig_kfbchi2_3tracks == NULL ) std::cerr << "[ABORT] can not find histgram 22 : " << std::endl, abort();
  if( hist_sig_kfbchi2_4tracks == NULL ) std::cerr << "[ABORT] can not find histgram 23 : " << std::endl, abort();
  if( hist_bkg_kfbchi2         == NULL ) std::cerr << "[ABORT] can not find histgram 24 : " << std::endl, abort();
  if( hist_bkg_kfbchi2_qq      == NULL ) std::cerr << "[ABORT] can not find histgram 25 : " << std::endl, abort();
  if( hist_bkg_kfbchi2_bb      == NULL ) std::cerr << "[ABORT] can not find histgram 26 : " << std::endl, abort();
  if( hist_sig_fmiss           == NULL ) std::cerr << "[ABORT] can not find histgram 27 : " << std::endl, abort();
  if( hist_sig_fmiss_qq        == NULL ) std::cerr << "[ABORT] can not find histgram 28 : " << std::endl, abort();
  if( hist_sig_fmiss_bb        == NULL ) std::cerr << "[ABORT] can not find histgram 29 : " << std::endl, abort();
  if( hist_bkg_fmiss           == NULL ) std::cerr << "[ABORT] can not find histgram 30 : " << std::endl, abort();
  if( hist_bkg_fmiss_qq        == NULL ) std::cerr << "[ABORT] can not find histgram 31 : " << std::endl, abort();
  if( hist_bkg_fmiss_bb        == NULL ) std::cerr << "[ABORT] can not find histgram 32 : " << std::endl, abort();
  
  TF1* func_sig_de_lep1         = (TF1*)file_de_lep1.Get     ( "pdf_func_de_lep1_sig"         );
  TF1* func_sig_de_lep1_pi0     = (TF1*)file_de_lep1.Get     ( "pdf_func_de_lep1_pi0_sig"     );
  TF1* func_sig_de_lep0         = (TF1*)file_de_lep0.Get     ( "pdf_func_de_lep0_sig"         );
  TF1* func_sig_de_lep0_pi0     = (TF1*)file_de_lep0.Get     ( "pdf_func_de_lep0_pi0_sig"     );
  TF1* func_bkg_de_lep1         = (TF1*)file_de_lep1.Get     ( "pdf_func_de_lep1_tot_bkg"     );
  TF1* func_bkg_de_lep1_qq      = (TF1*)file_de_lep1.Get     ( "pdf_func_de_lep1_qq_bkg"      );
  TF1* func_bkg_de_lep1_bb      = (TF1*)file_de_lep1.Get     ( "pdf_func_de_lep1_bb_bkg"      );
  TF1* func_bkg_de_lep0         = (TF1*)file_de_lep0.Get     ( "pdf_func_de_lep0_tot_bkg"     );
  TF1* func_bkg_de_lep0_qq      = (TF1*)file_de_lep0.Get     ( "pdf_func_de_lep0_qq_bkg"      );
  TF1* func_bkg_de_lep0_bb      = (TF1*)file_de_lep0.Get     ( "pdf_func_de_lep0_bb_bkg"      );
  TF1* func_sig_dzll3d          = (TF1*)file_dzll3d.Get      ( "pdf_func_dzll3d_sig"          );
  TF1* func_bkg_dzll3d          = (TF1*)file_dzll3d.Get      ( "pdf_func_dzll3d_tot_bkg"      );
  TF1* func_bkg_dzll3d_qq       = (TF1*)file_dzll3d.Get      ( "pdf_func_dzll3d_qq_bkg"       );
  TF1* func_bkg_dzll3d_bb       = (TF1*)file_dzll3d.Get      ( "pdf_func_dzll3d_bb_bkg"       );
  TF1* func_sig_bccm            = (TF1*)file_bccm.Get        ( "pdf_func_bccm_sig"            );
  TF1* func_bkg_bccm            = (TF1*)file_bccm.Get        ( "pdf_func_bccm_tot_bkg"        );
  TF1* func_bkg_bccm_qq         = (TF1*)file_bccm.Get        ( "pdf_func_bccm_qq_bkg"         );
  TF1* func_bkg_bccm_bb         = (TF1*)file_bccm.Get        ( "pdf_func_bccm_bb_bkg"         );
  TF1* func_sig_kfbchi2_0tracks = (TF1*)file_sig_kfbchi2.Get ( "pdf_func_kfbchi2_0tracks_sig" );
  TF1* func_sig_kfbchi2_1tracks = (TF1*)file_sig_kfbchi2.Get ( "pdf_func_kfbchi2_1tracks_sig" );
  TF1* func_sig_kfbchi2_2tracks = (TF1*)file_sig_kfbchi2.Get ( "pdf_func_kfbchi2_2tracks_sig" );
  TF1* func_sig_kfbchi2_3tracks = (TF1*)file_sig_kfbchi2.Get ( "pdf_func_kfbchi2_3tracks_sig" );
  TF1* func_sig_kfbchi2_4tracks = (TF1*)file_sig_kfbchi2.Get ( "pdf_func_kfbchi2_4tracks_sig" );
  TF1* func_bkg_kfbchi2         = (TF1*)file_bkg_kfbchi2.Get ( "pdf_func_kfbchi2_tot_bkg"     );
  TF1* func_bkg_kfbchi2_qq      = (TF1*)file_bkg_kfbchi2.Get ( "pdf_func_kfbchi2_qq_bkg"      );
  TF1* func_bkg_kfbchi2_bb      = (TF1*)file_bkg_kfbchi2.Get ( "pdf_func_kfbchi2_bb_bkg"      );
  //TF1* func_sig_fmiss         = (TF1*)file_sig_fmiss.Get   ( "pdf_func_fmiss_sig"           );
  //TF1* func_sig_fmiss_qq      = (TF1*)file_sig_fmiss_qq.Get( "pdf_func_fmiss_qq_sig"        );
  //TF1* func_sig_fmiss_bb      = (TF1*)file_sig_fmiss_bb.Get( "pdf_func_fmiss_bb_sig"        );
  //TF1* func_bkg_fmiss         = (TF1*)file_bkg_fmiss.Get   ( "pdf_func_fmiss_bkg"           );
  //TF1* func_bkg_fmiss_qq      = (TF1*)file_bkg_fmiss_qq.Get( "pdf_func_fmiss_qq_bkg"        );
  //TF1* func_bkg_fmiss_bb      = (TF1*)file_bkg_fmiss_bb.Get( "pdf_func_fmiss_bb_bkg"        );
  
  if( func_sig_de_lep1         == NULL ) std::cerr << "[ABORT] can not find func  1 : " << std::endl, abort();
  if( func_sig_de_lep1_pi0     == NULL ) std::cerr << "[ABORT] can not find func  2 : " << std::endl, abort();
  if( func_sig_de_lep0         == NULL ) std::cerr << "[ABORT] can not find func  3 : " << std::endl, abort();
  if( func_sig_de_lep0_pi0     == NULL ) std::cerr << "[ABORT] can not find func  4 : " << std::endl, abort();
  if( func_bkg_de_lep1         == NULL ) std::cerr << "[ABORT] can not find func  5 : " << std::endl, abort();
  if( func_bkg_de_lep1_qq      == NULL ) std::cerr << "[ABORT] can not find func  6 : " << std::endl, abort();
  if( func_bkg_de_lep1_bb      == NULL ) std::cerr << "[ABORT] can not find func  7 : " << std::endl, abort();
  if( func_bkg_de_lep0         == NULL ) std::cerr << "[ABORT] can not find func  8 : " << std::endl, abort();
  if( func_bkg_de_lep0_qq      == NULL ) std::cerr << "[ABORT] can not find func  9 : " << std::endl, abort();
  if( func_bkg_de_lep0_bb      == NULL ) std::cerr << "[ABORT] can not find func 10 : " << std::endl, abort();
  if( func_sig_dzll3d          == NULL ) std::cerr << "[ABORT] can not find func 11 : " << std::endl, abort();
  if( func_bkg_dzll3d          == NULL ) std::cerr << "[ABORT] can not find func 12 : " << std::endl, abort();
  if( func_bkg_dzll3d_qq       == NULL ) std::cerr << "[ABORT] can not find func 13 : " << std::endl, abort();
  if( func_bkg_dzll3d_bb       == NULL ) std::cerr << "[ABORT] can not find func 14 : " << std::endl, abort();
  if( func_sig_bccm            == NULL ) std::cerr << "[ABORT] can not find func 15 : " << std::endl, abort();
  if( func_bkg_bccm            == NULL ) std::cerr << "[ABORT] can not find func 16 : " << std::endl, abort();
  if( func_bkg_bccm_qq         == NULL ) std::cerr << "[ABORT] can not find func 17 : " << std::endl, abort();
  if( func_bkg_bccm_bb         == NULL ) std::cerr << "[ABORT] can not find func 18 : " << std::endl, abort();
  if( func_sig_kfbchi2_0tracks == NULL ) std::cerr << "[ABORT] can not find func 19 : " << std::endl, abort();
  if( func_sig_kfbchi2_1tracks == NULL ) std::cerr << "[ABORT] can not find func 20 : " << std::endl, abort();
  if( func_sig_kfbchi2_2tracks == NULL ) std::cerr << "[ABORT] can not find func 21 : " << std::endl, abort();
  if( func_sig_kfbchi2_3tracks == NULL ) std::cerr << "[ABORT] can not find func 22 : " << std::endl, abort();
  if( func_sig_kfbchi2_4tracks == NULL ) std::cerr << "[ABORT] can not find func 23 : " << std::endl, abort();
  if( func_bkg_kfbchi2         == NULL ) std::cerr << "[ABORT] can not find func 24 : " << std::endl, abort();
  if( func_bkg_kfbchi2_qq      == NULL ) std::cerr << "[ABORT] can not find func 25 : " << std::endl, abort();
  if( func_bkg_kfbchi2_bb      == NULL ) std::cerr << "[ABORT] can not find func 26 : " << std::endl, abort();
  //if( func_sig_fmiss         == NULL ) std::cerr << "[ABORT] can not find func 27 : " << std::endl, abort();
  //if( func_sig_fmiss_sig     == NULL ) std::cerr << "[ABORT] can not find func 28 : " << std::endl, abort();
  //if( func_sig_fmiss_bkg     == NULL ) std::cerr << "[ABORT] can not find func 29 : " << std::endl, abort();
  //if( func_bkg_fmiss         == NULL ) std::cerr << "[ABORT] can not find func 30 : " << std::endl, abort();
  //if( func_bkg_fmiss_qq      == NULL ) std::cerr << "[ABORT] can not find func 31 : " << std::endl, abort();
  //if( func_bkg_fmiss_bb      == NULL ) std::cerr << "[ABORT] can not find func 32 : " << std::endl, abort();

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  if( fl_mc==3 ){
    TCanvas*  c1      = Canvas( "c1","c1",8 );
    TCanvas*  c2      = Canvas( "c2","c2",8 );
    if( fl_de ){
      c1->cd(1);
      Deco( hist_sig_de_lep1,     2, 2, 2 ); hist_sig_de_lep1    ->Draw(      );
      Deco( hist_sig_de_lep1_pi0, 2, 3, 3 ); hist_sig_de_lep1_pi0->Draw("same");
      Deco( hist_bkg_de_lep1,     2, 7, 7 ); hist_bkg_de_lep1    ->Draw("same");
      Deco( hist_bkg_de_lep1_qq,  2, 8, 8 ); hist_bkg_de_lep1_qq ->Draw("same");
      Deco( hist_bkg_de_lep1_bb,  2, 9, 9 ); hist_bkg_de_lep1_bb ->Draw("same");
      if( fl_de==2 ){
	Deco( func_sig_de_lep1,     0, 2, 2 ); func_sig_de_lep1    ->Draw("same");
	Deco( func_sig_de_lep1_pi0, 0, 3, 3 ); func_sig_de_lep1_pi0->Draw("same");
	Deco( func_bkg_de_lep1,     0, 7, 7 ); func_bkg_de_lep1    ->Draw("same");
	Deco( func_bkg_de_lep1_qq,  0, 8, 8 ); func_bkg_de_lep1_qq ->Draw("same");
	Deco( func_bkg_de_lep1_bb,  0, 9, 9 ); func_bkg_de_lep1_bb ->Draw("same");
	
	c2->cd(1);
	Deco( func_sig_de_lep1,     1, 2, 2 ); func_sig_de_lep1    ->Draw(      );
	Deco( func_sig_de_lep1_pi0, 1, 3, 3 ); func_sig_de_lep1_pi0->Draw("same");
	Deco( func_bkg_de_lep1,     1, 7, 7 ); func_bkg_de_lep1    ->Draw("same");
	Deco( func_bkg_de_lep1_qq,  1, 8, 8 ); func_bkg_de_lep1_qq ->Draw("same");
	Deco( func_bkg_de_lep1_bb,  1, 9, 9 ); func_bkg_de_lep1_bb ->Draw("same");
      }
      c1->cd(2);
      Deco( hist_sig_de_lep0,     2, 2, 2 ); hist_sig_de_lep0    ->Draw(      );
      Deco( hist_sig_de_lep0_pi0, 2, 3, 3 ); hist_sig_de_lep0_pi0->Draw("same");
      Deco( hist_bkg_de_lep0,     2, 7, 7 ); hist_bkg_de_lep0    ->Draw("same");
      Deco( hist_bkg_de_lep0_qq,  2, 8, 8 ); hist_bkg_de_lep0_qq ->Draw("same");
      Deco( hist_bkg_de_lep0_bb,  2, 9, 9 ); hist_bkg_de_lep0_bb ->Draw("same");
      
      if( fl_de==2 ){
	Deco( func_sig_de_lep0,     0, 2, 2 ); func_sig_de_lep0    ->Draw("same");
	Deco( func_sig_de_lep0_pi0, 0, 3, 3 ); func_sig_de_lep0_pi0->Draw("same");
	Deco( func_bkg_de_lep0,     0, 7, 7 ); func_bkg_de_lep0    ->Draw("same");
	Deco( func_bkg_de_lep0_qq,  0, 8, 8 ); func_bkg_de_lep0_qq ->Draw("same");
	Deco( func_bkg_de_lep0_bb,  0, 9, 9 ); func_bkg_de_lep0_bb ->Draw("same");
	c2->cd(2);
	Deco( func_sig_de_lep0,     1, 2, 2 ); func_sig_de_lep0    ->Draw(      );
	Deco( func_sig_de_lep0_pi0, 1, 3, 3 ); func_sig_de_lep0_pi0->Draw("same");
	Deco( func_bkg_de_lep0,     1, 7, 7 ); func_bkg_de_lep0    ->Draw("same");
	Deco( func_bkg_de_lep0_qq,  1, 8, 8 ); func_bkg_de_lep0_qq ->Draw("same");
	Deco( func_bkg_de_lep0_bb,  1, 9, 9 ); func_bkg_de_lep0_bb ->Draw("same");
      }
    }
    if( fl_kfbchi2 ){
      c1->cd(3);
      Deco( hist_sig_kfbchi2_0tracks, 2, 2, 2 ); hist_sig_kfbchi2_0tracks->Draw(      );
      Deco( hist_sig_kfbchi2_1tracks, 2, 3, 3 ); hist_sig_kfbchi2_1tracks->Draw("same");
      Deco( hist_sig_kfbchi2_2tracks, 2, 4, 4 ); hist_sig_kfbchi2_2tracks->Draw("same");
      Deco( hist_sig_kfbchi2_3tracks, 2, 5, 5 ); hist_sig_kfbchi2_3tracks->Draw("same");
      Deco( hist_sig_kfbchi2_4tracks, 2, 6, 6 ); hist_sig_kfbchi2_4tracks->Draw("same");
      Deco( hist_bkg_kfbchi2,         2, 7, 7 ); hist_bkg_kfbchi2        ->Draw("same");
      Deco( hist_bkg_kfbchi2_qq,      2, 8, 8 ); hist_bkg_kfbchi2_qq     ->Draw("same");
      Deco( hist_bkg_kfbchi2_bb,      2, 9, 9 ); hist_bkg_kfbchi2_bb     ->Draw("same");
      if( fl_kfbchi2==2 ){
	Deco( func_sig_kfbchi2_0tracks, 0, 2, 2 ); func_sig_kfbchi2_0tracks->Draw("same");
	Deco( func_sig_kfbchi2_1tracks, 0, 3, 3 ); func_sig_kfbchi2_1tracks->Draw("same");
	Deco( func_sig_kfbchi2_2tracks, 0, 4, 4 ); func_sig_kfbchi2_2tracks->Draw("same");
	Deco( func_sig_kfbchi2_3tracks, 0, 5, 5 ); func_sig_kfbchi2_3tracks->Draw("same");
	Deco( func_sig_kfbchi2_4tracks, 0, 6, 6 ); func_sig_kfbchi2_4tracks->Draw("same");
	Deco( func_bkg_kfbchi2,         0, 7, 7 ); func_bkg_kfbchi2        ->Draw("same");
	Deco( func_bkg_kfbchi2_qq,      0, 8, 8 ); func_bkg_kfbchi2_qq     ->Draw("same");
	Deco( func_bkg_kfbchi2_bb,      0, 9, 9 ); func_bkg_kfbchi2_bb     ->Draw("same");
	c2->cd(3);
	Deco( func_sig_kfbchi2_0tracks, 1, 2, 2 ); func_sig_kfbchi2_0tracks->Draw(      );
	Deco( func_sig_kfbchi2_1tracks, 1, 3, 3 ); func_sig_kfbchi2_1tracks->Draw("same");
	Deco( func_sig_kfbchi2_2tracks, 1, 4, 4 ); func_sig_kfbchi2_2tracks->Draw("same");
	Deco( func_sig_kfbchi2_3tracks, 1, 5, 5 ); func_sig_kfbchi2_3tracks->Draw("same");
	Deco( func_sig_kfbchi2_4tracks, 1, 6, 6 ); func_sig_kfbchi2_4tracks->Draw("same");
	Deco( func_bkg_kfbchi2,         1, 7, 7 ); func_bkg_kfbchi2        ->Draw("same");
	Deco( func_bkg_kfbchi2_qq,      1, 8, 8 ); func_bkg_kfbchi2_qq     ->Draw("same");
	Deco( func_bkg_kfbchi2_bb,      1, 9, 9 ); func_bkg_kfbchi2_bb     ->Draw("same");
      }
    }
    if( fl_dzll3d ){
      c1->cd(4);
      Deco( hist_sig_dzll3d,    2, 2, 2 ); hist_sig_dzll3d   ->Draw(      );
      Deco( hist_bkg_dzll3d,    2, 7, 7 ); hist_bkg_dzll3d   ->Draw("same");
      Deco( hist_bkg_dzll3d_qq, 2, 8, 8 ); hist_bkg_dzll3d_qq->Draw("same");
      Deco( hist_bkg_dzll3d_bb, 2, 9, 9 ); hist_bkg_dzll3d_bb->Draw("same");
      // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
      TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
      legend1->AddEntry( hist_bkg_dzll3d,    "bkg(total)", "L" );
      legend1->AddEntry( hist_bkg_dzll3d_qq, "bkg( qq  )", "L" );
      legend1->AddEntry( hist_bkg_dzll3d_bb, "bkg( bb  )", "L" );
      legend1->Draw();
      if( fl_dzll3d==2 ){
	Deco( func_sig_dzll3d,    0, 2, 2 ); func_sig_dzll3d   ->Draw("same");
	Deco( func_bkg_dzll3d,    0, 7, 7 ); func_bkg_dzll3d   ->Draw("same");
	Deco( func_bkg_dzll3d_qq, 0, 8, 8 ); func_bkg_dzll3d_qq->Draw("same");
	Deco( func_bkg_dzll3d_bb, 0, 9, 9 ); func_bkg_dzll3d_bb->Draw("same");
	c2->cd(4);
	Deco( func_sig_dzll3d,    1, 2, 2 ); func_sig_dzll3d   ->Draw(      );
	Deco( func_bkg_dzll3d,    1, 7, 7 ); func_bkg_dzll3d   ->Draw("same");
	Deco( func_bkg_dzll3d_qq, 1, 8, 8 ); func_bkg_dzll3d_qq->Draw("same");
	Deco( func_bkg_dzll3d_bb, 1, 9, 9 ); func_bkg_dzll3d_bb->Draw("same");
      }
    }
    if( fl_bccm ){
      c1->cd(5);
      Deco( hist_sig_bccm,    2, 2, 2 ); hist_sig_bccm   ->Draw(      );
      Deco( hist_bkg_bccm,    2, 7, 7 ); hist_bkg_bccm   ->Draw("same");
      Deco( hist_bkg_bccm_qq, 2, 8, 8 ); hist_bkg_bccm_qq->Draw("same");
      Deco( hist_bkg_bccm_bb, 2, 9, 9 ); hist_bkg_bccm_bb->Draw("same");
      if( fl_bccm==2 ){
	Deco( func_sig_bccm,    0, 2, 2 ); func_sig_bccm   ->Draw("same");
	Deco( func_bkg_bccm,    0, 7, 7 ); func_bkg_bccm   ->Draw("same");
	Deco( func_bkg_bccm_qq, 0, 8, 8 ); func_bkg_bccm_qq->Draw("same");
	Deco( func_bkg_bccm_bb, 0, 9, 9 ); func_bkg_bccm_bb->Draw("same");
	c2->cd(5);
	Deco( func_sig_bccm,    1, 2, 2 ); func_sig_bccm   ->Draw(      );
	Deco( func_bkg_bccm,    1, 7, 7 ); func_bkg_bccm   ->Draw("same");
	Deco( func_bkg_bccm_qq, 1, 8, 8 ); func_bkg_bccm_qq->Draw("same");
	Deco( func_bkg_bccm_bb, 1, 9, 9 ); func_bkg_bccm_bb->Draw("same");
      }
    }
    if( fl_fmiss ){
      c1->cd(6);
      Deco( hist_sig_fmiss_qq, 2,  2,  2 ); hist_sig_fmiss_qq->DrawNormalized(      );
      Deco( hist_sig_fmiss_bb, 2,  3,  3 ); hist_sig_fmiss_bb->DrawNormalized("same");
      Deco( hist_sig_fmiss,    2,  1,  1 ); hist_sig_fmiss   ->DrawNormalized("same");
      Deco( hist_bkg_fmiss_qq, 2,  6,  6 ); hist_bkg_fmiss_qq->DrawNormalized("same");
      Deco( hist_bkg_fmiss_bb, 2,  7,  7 ); hist_bkg_fmiss_bb->DrawNormalized("same");
      Deco( hist_bkg_fmiss,    2, 11, 11 ); hist_bkg_fmiss   ->DrawNormalized("same");
      // +++++++ tlegend2 ++++++++++++++++++++++++++++++++++
      TLegend* legend2 = new TLegend( 0.75,0.75,0.99,0.99 );
      legend2->AddEntry( hist_sig_fmiss_qq, "sig( qq  )", "L" );
      legend2->AddEntry( hist_sig_fmiss_bb, "sig( bb  )", "L" );
      legend2->AddEntry( hist_sig_fmiss,    "sig(total)", "L" );
      legend2->AddEntry( hist_bkg_fmiss_qq, "bkg( qq  )", "L" );
      legend2->AddEntry( hist_bkg_fmiss_bb, "bkg( bb  )", "L" );
      legend2->AddEntry( hist_bkg_fmiss,    "bkg(total)", "L" );
      legend2->Draw();
      //if( fl_fmiss==2 ){
      //Deco( func_sig_fmiss,    0, 2, 2 ); func_sig_fmiss   ->Draw("same");
      //Deco( func_bkg_fmiss,    0, 7, 7 ); func_bkg_fmiss   ->Draw("same");
      //Deco( func_bkg_fmiss_qq, 0, 8, 8 ); func_bkg_fmiss_qq->Draw("same");
      //Deco( func_bkg_fmiss_bb, 0, 9, 9 ); func_bkg_fmiss_bb->Draw("same");
      c2->cd(6);
      //Deco( func_sig_fmiss,    1, 2, 2 ); func_sig_fmiss   ->Draw(      );
      //Deco( func_bkg_fmiss,    1, 7, 7 ); func_bkg_fmiss   ->Draw("same");
      //Deco( func_bkg_fmiss_qq, 1, 8, 8 ); func_bkg_fmiss_qq->Draw("same");
      //Deco( func_bkg_fmiss_bb, 1, 9, 9 ); func_bkg_fmiss_bb->Draw("same");
      //}
    }
    c1->Print( Form("pic/pic_mklikelihood_hist_%s.root", tag) );
    c2->Print( Form("pic/pic_mklikelihood_func_%s.root", tag) );
    c1->Print( Form("pic/pic_mklikelihood_hist_%s.eps",  tag) );
    c2->Print( Form("pic/pic_mklikelihood_func_%s.eps",  tag) );
    return 0;
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TFile* rootf = new TFile( Form("%s%s",outdir, basename), "RECREATE" );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  TChain* chain = new TChain( treename );
  Int_t   nfile = chain->Add( filename );
  TTree*  tree  = new TTree();
  tree = chain->CloneTree();
  /*
  MCut_array* cut = new MCut_array( branch_table() );
  make_cut_ee( cut, treename );
  cut->Set( "Mbc",  0 );
  cut->Set( "de",   0 );
  TCut cut_ee = cut->Output();
  make_cut_mm( cut, treename );
  cut->Set( "Mbc",  0 );
  cut->Set( "de",   0 );
  TCut cut_mm = cut->Output();
  tree = chain->CopyTree( cut_ee || cut_mm );
  */
  //+++++++++++++++++++++++++++++++++++++++++++++++++
  Float_t lr    = 0;
  Float_t lr_qq = 0;
  Float_t lr_bb = 0;
  Float_t ls    = 1;
  Float_t ls_qq = 1;
  Float_t ls_bb = 1;
  Float_t lb    = 1;
  Float_t lb_qq = 1;
  Float_t lb_bb = 1;
  tree->Branch( Form("lr%s_tot",tag), &lr,    Form("lr%s_tot/F",tag) );
  tree->Branch( Form("ls%s_tot",tag), &ls,    Form("ls%s_tot/F",tag) );
  tree->Branch( Form("lb%s_tot",tag), &lb,    Form("lb%s_tot/F",tag) );
  tree->Branch( Form("lr%s_qq", tag), &lr_qq, Form("lr%s_qq/F", tag) );
  tree->Branch( Form("ls%s_qq", tag), &ls_qq, Form("ls%s_qq/F", tag) );
  tree->Branch( Form("lb%s_qq", tag), &lb_qq, Form("lb%s_qq/F", tag) );
  tree->Branch( Form("lr%s_bb", tag), &lr_bb, Form("lr%s_bb/F", tag) );
  tree->Branch( Form("ls%s_bb", tag), &ls_bb, Form("ls%s_bb/F", tag) );
  tree->Branch( Form("lb%s_bb", tag), &lb_bb, Form("lb%s_bb/F", tag) );
  TTree* newtree = tree->CloneTree( 0, "newtree" );
  //+++++++++++++++++++++++++++++++++++++++++++++++++

  Float_t  de, kfbchi, kfbdgf, dzll3d, bccm;
  Float_t  event, rm_l, rm_xs;
  Float_t ksfwsig, ksfwsig_qq, ksfwsig_bb;
  Float_t ksfwbkg, ksfwbkg_qq, ksfwbkg_bb;
  Float_t fmiss,   fmiss_qq,   fmiss_bb;

  tree->SetBranchAddress( "event",      &event      );
  tree->SetBranchAddress( "rm_l",       &rm_l       );
  tree->SetBranchAddress( "rm_xs",      &rm_xs      );
  
  tree->SetBranchAddress( "de",         &de         );
  tree->SetBranchAddress( "kfbchi",     &kfbchi     );
  tree->SetBranchAddress( "kfbdgf",     &kfbdgf     );
  tree->SetBranchAddress( "dzll3d",     &dzll3d     );
  tree->SetBranchAddress( "bccm",       &bccm       );
  tree->SetBranchAddress( "ksfwsig_tot",&ksfwsig    );
  tree->SetBranchAddress( "ksfwsig_qq", &ksfwsig_qq );
  tree->SetBranchAddress( "ksfwsig_bb", &ksfwsig_bb );
  tree->SetBranchAddress( "ksfwbkg_tot",&ksfwbkg    );
  tree->SetBranchAddress( "ksfwbkg_qq", &ksfwbkg_qq );
  tree->SetBranchAddress( "ksfwbkg_bb", &ksfwbkg_bb );
  tree->SetBranchAddress( Form("fmiss%s_tot_deroe",tag),  &fmiss    ); // deroe
  tree->SetBranchAddress( Form("fmiss%s_qq_deroe", tag),  &fmiss_qq ); // deroe
  tree->SetBranchAddress( Form("fmiss%s_bb_deroe", tag),  &fmiss_bb ); // deroe
  //tree->SetBranchAddress( Form("fmiss%s_tot",tag),  &fmiss    );
  //tree->SetBranchAddress( Form("fmiss%s_qq", tag),  &fmiss_qq );
  //tree->SetBranchAddress( Form("fmiss%s_bb", tag),  &fmiss_bb );

  //+++++++++++++++++++++++++++++++++++++++++++++++++
  Int_t    nev = 0;
  while (tree->GetEntry(nev, 0)) {
    ls     = 1;
    lb     = 1;
    ls_qq  = 1;
    lb_qq  = 1;
    ls_bb  = 1;
    lb_bb  = 1;
    if( fl_message ) std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
			       << "[ "  << std::setw(7) << std::right << nev
			       << "(evtNo = " << std::setw(7) << std::right << event << ") ] " << std::endl
			       << std::setw(10) << std::right <<          "de = " << std::setw(10) << std::right << de         << ", "
			       << std::setw(10) << std::right <<      "kfbchi = " << std::setw(10) << std::right << kfbchi     << ", "
			       << std::setw(10) << std::right <<      "kfbdgf = " << std::setw(10) << std::right << kfbdgf     << ", "
			       << std::setw(10) << std::right <<      "dzll3d = " << std::setw(10) << std::right << dzll3d     << ", "
			       << std::setw(10) << std::right <<        "bccm = " << std::setw(10) << std::right << bccm       << std::endl
			       << std::setw(10) << std::right <<     "ksfwsig = " << std::setw(10) << std::right << ksfwsig    << ", "
      			       << std::setw(10) << std::right <<  "ksfwsig_qq = " << std::setw(10) << std::right << ksfwsig_qq << ", "
      			       << std::setw(10) << std::right <<  "ksfwsig_bb = " << std::setw(10) << std::right << ksfwsig_bb << ", "
			       << std::setw(10) << std::right << "ksfwbkg_tot = " << std::setw(10) << std::right << ksfwbkg    << ", "
      			       << std::setw(10) << std::right <<  "ksfwbkg_qq = " << std::setw(10) << std::right << ksfwbkg_qq << ", "
      			       << std::setw(10) << std::right <<  "ksfwbkg_bb = " << std::setw(10) << std::right << ksfwbkg_bb << std::endl
			       << std::setw(10) << std::right <<   "fmiss_tot = " << std::setw(10) << std::right << fmiss      << ", "
    			       << std::setw(10) << std::right <<    "fmiss_qq = " << std::setw(10) << std::right << fmiss_qq   << ", "
    			       << std::setw(10) << std::right <<    "fmiss_bb = " << std::setw(10) << std::right << fmiss_bb   << std::endl;
    
    // de(sig)
    if( fl_message ) std::cout << "[ de(sig) ]" << std::endl;
    if( rm_l==1 ){ // ee
      if( rm_xs > 999 ){
	ls *= Eval( fl_de, de, hist_sig_de_lep1_pi0, func_sig_de_lep1_pi0, fl_message ); // w/  pi0
      }else{
	ls *= Eval( fl_de, de, hist_sig_de_lep1,     func_sig_de_lep1,     fl_message ); // w/o pi0
      }
    }else if( rm_l==0 ){ // mu
      if( rm_xs > 999 ){
	ls *= Eval( fl_de, de, hist_sig_de_lep0_pi0, func_sig_de_lep0_pi0, fl_message ); // w/  pi0
      }else{
	ls *= Eval( fl_de, de, hist_sig_de_lep0,     func_sig_de_lep0,     fl_message ); // w/o pi0
      }
    }

    // de(bkg)
    if( fl_message ) std::cout << "[ de(bkg) ]" << std::endl;
    if( rm_l==1 ){ // ee
      lb    *= Eval( fl_de, de, hist_bkg_de_lep1,    func_bkg_de_lep1,    fl_message );
      lb_qq *= Eval( fl_de, de, hist_bkg_de_lep1_qq, func_bkg_de_lep1_qq, fl_message );
      lb_bb *= Eval( fl_de, de, hist_bkg_de_lep1_bb, func_bkg_de_lep1_bb, fl_message );
    }else if( rm_l==0 ){ // mu
      lb    *= Eval( fl_de, de, hist_bkg_de_lep0,    func_bkg_de_lep0,    fl_message );
      lb_qq *= Eval( fl_de, de, hist_bkg_de_lep0_qq, func_bkg_de_lep0_qq, fl_message );
      lb_bb *= Eval( fl_de, de, hist_bkg_de_lep0_bb, func_bkg_de_lep0_bb, fl_message );
    }
   
    // kfbchi2(sig)
    if( fl_message ) std::cout << "[ kfbchi2(sig) ]" << std::endl;
    if(                                    rm_xs==  10 || rm_xs==1010 ) ls *= Eval( fl_kfbchi2, kfbchi/kfbdgf, hist_sig_kfbchi2_0tracks, func_sig_kfbchi2_0tracks, fl_message );// 0tracks
    else if( rm_xs==   1 || rm_xs==1001 || rm_xs== 110 || rm_xs==1110 ) ls *= Eval( fl_kfbchi2, kfbchi/kfbdgf, hist_sig_kfbchi2_1tracks, func_sig_kfbchi2_1tracks, fl_message );// 1tracks
    else if( rm_xs== 101 || rm_xs==1101 || rm_xs== 210 || rm_xs==1210 ) ls *= Eval( fl_kfbchi2, kfbchi/kfbdgf, hist_sig_kfbchi2_2tracks, func_sig_kfbchi2_2tracks, fl_message ); // 2tracks
    else if( rm_xs== 201 || rm_xs==1201 || rm_xs== 310 || rm_xs==1310 ) ls *= Eval( fl_kfbchi2, kfbchi/kfbdgf, hist_sig_kfbchi2_3tracks, func_sig_kfbchi2_3tracks, fl_message ); // 3tracks 
    else if( rm_xs== 301 || rm_xs==1301 || rm_xs== 410 || rm_xs== 401 ) ls *= Eval( fl_kfbchi2, kfbchi/kfbdgf, hist_sig_kfbchi2_4tracks, func_sig_kfbchi2_4tracks, fl_message ); // [45]tracks
   
    // kfbchi2(bkg)
    if( fl_message ) std::cout << "[ kfbchi2(bkg) ]" << std::endl;
    lb    *= Eval( fl_kfbchi2, kfbchi/kfbdgf, hist_bkg_kfbchi2,    func_bkg_kfbchi2,    fl_message );
    lb_qq *= Eval( fl_kfbchi2, kfbchi/kfbdgf, hist_bkg_kfbchi2_qq, func_bkg_kfbchi2_qq, fl_message );
    lb_bb *= Eval( fl_kfbchi2, kfbchi/kfbdgf, hist_bkg_kfbchi2_bb, func_bkg_kfbchi2_bb, fl_message );
   
    // dzll3d(sig&bkg)
    if( fl_message ) std::cout << "[ dzll3d ]" << std::endl;
    ls    *= Eval( fl_dzll3d, dzll3d, hist_sig_dzll3d,    func_sig_dzll3d,    fl_message );
    lb    *= Eval( fl_dzll3d, dzll3d, hist_bkg_dzll3d,    func_bkg_dzll3d,    fl_message );
    lb_qq *= Eval( fl_dzll3d, dzll3d, hist_bkg_dzll3d_qq, func_bkg_dzll3d_qq, fl_message );
    lb_bb *= Eval( fl_dzll3d, dzll3d, hist_bkg_dzll3d_bb, func_bkg_dzll3d_bb, fl_message );
   
    // cos-theta_B(sig&bkg)
    if( fl_message ) std::cout << "[ cos-theta_B ]" << std::endl;
    ls    *= Eval( fl_bccm, bccm, hist_sig_bccm,    func_sig_bccm,    fl_message );
    lb    *= Eval( fl_bccm, bccm, hist_bkg_bccm,    func_bkg_bccm,    fl_message );
    lb_qq *= Eval( fl_bccm, bccm, hist_bkg_bccm_qq, func_bkg_bccm_qq, fl_message );
    lb_bb *= Eval( fl_bccm, bccm, hist_bkg_bccm_bb, func_bkg_bccm_bb, fl_message );
   
    ////////////////////////
    ls_qq= ls;
    ls_bb= ls;
    ////////////////////////
    
    // KSFW(sig&bkg)
    if( fl_ksfw ){
      if( fl_message ) std::cout << "[ KSFW ]"   << std::endl
				 << "  (total) " << std::setw(10) << std::right << ksfwsig    << ", " << std::setw(10) << std::right << ksfwbkg    << std::endl
				 << "  ( qq  ) " << std::setw(10) << std::right << ksfwsig_qq << ", " << std::setw(10) << std::right << ksfwbkg_qq << std::endl
				 << "  ( bb  ) " << std::setw(10) << std::right << ksfwsig_bb << ", " << std::setw(10) << std::right << ksfwbkg_bb << std::endl;
      ls    *= ksfwsig;
      ls_qq *= ksfwsig_qq;
      ls_bb *= ksfwsig_bb;
      lb    *= ksfwbkg;
      lb_qq *= ksfwbkg_qq;
      lb_bb *= ksfwbkg_bb;
    }
    
    // Fmiss(sig&bkg)
    if( fl_message ) std::cout << "[ Fmiss ]" << std::endl;
    ls    *= Eval( fl_fmiss, fmiss,    hist_sig_fmiss,    NULL, fl_message );
    ls_qq *= Eval( fl_fmiss, fmiss_qq, hist_sig_fmiss_qq, NULL, fl_message );
    ls_bb *= Eval( fl_fmiss, fmiss_bb, hist_sig_fmiss_bb, NULL, fl_message );
    lb    *= Eval( fl_fmiss, fmiss,    hist_bkg_fmiss,    NULL, fl_message );
    lb_qq *= Eval( fl_fmiss, fmiss_qq, hist_bkg_fmiss_qq, NULL, fl_message );
    lb_bb *= Eval( fl_fmiss, fmiss_bb, hist_bkg_fmiss_bb, NULL, fl_message );
    
    if( ls || lb ) lr = ls/( ls+lb );
    else lr = -1;
    if( ls_qq || lb_qq ) lr_qq = ls_qq/( ls_qq+lb_qq );
    else lr_qq = -1;
    if( ls_bb || lb_bb ) lr_bb = ls_bb/( ls_bb+lb_bb );
    else lr_bb = -1;
    if( fl_message ){
      std::cout << "lr    = " << ls    << "/(" << ls    << " + " << lb    << " ) = " << lr    << std::endl
		<< "lr_qq = " << ls_qq << "/(" << ls_qq << " + " << lb_qq << " ) = " << lr_qq << std::endl
		<< "lr_bb = " << ls_bb << "/(" << ls_bb << " + " << lb_qq << " ) = " << lr_bb << std::endl;
    }
    if( lr    > 1.0 ) std::cerr << "Strange Events having LR>1.0 Exist (lr = " << lr << " ) " << filename << std::endl;
    if( lr_qq > 1.0 ) std::cerr << "Strange Events having LR>1.0 Exist (lr = " << lr << " ) " << filename << std::endl;
    if( lr_bb > 1.0 ) std::cerr << "Strange Events having LR>1.0 Exist (lr = " << lr << " ) " << filename << std::endl;
    
    nev++;
    newtree->Fill();
  }
  
  std::cout << "mklikelihood : "
	    << std::setw(40) << std::right << basename
	    << " : "
	    << nev <<  " events processed"
	    << std::endl;
  
  newtree->Write();
  rootf->Close();
  return 0;
}


Double_t Eval( Int_t fl_pdf, Double_t fvar, TH1D* hist, TF1* func,
	       Bool_t fl_message, Bool_t fl_range, Double_t fxmin, Double_t fxmax ){
  if(      fl_pdf == 1 ) return Eval_hist( fvar, hist, fl_message, fl_range, fxmin, fxmax );
  else if( fl_pdf == 2 ) return Eval_func( fvar, func, fl_message, fl_range, fxmin, fxmax );
  else return 1;
}

Double_t Eval_hist( Double_t fvar, TH1D* hist, Bool_t fl_message, Bool_t fl_range, Double_t fxmin, Double_t fxmax ){
  Double_t area = hist->Integral( 1, hist->GetNbinsX() );
  if( fl_range ) area = hist->Integral( hist->FindBin(fxmin+0.0000001), hist->FindBin(fxmax-0.0000001) );

  Double_t bin  = hist->GetBinWidth( 1 );
  Double_t var  = hist->GetBinContent( hist->FindBin (fvar) );
  if( fl_message ) std::cout << "       "
			     << "var  = "  << std::setw(10) << std::right << var  
			     << " ("       << std::setw(10) << std::right << fvar << "), "
			     << "bin  = "  << std::setw( 4) << std::right << bin  << ", " 
			     << "entry = " << std::setw(10) << std::right << hist->GetEntries() << ", "
			     << "area = "  << std::setw(10) << std::right << area
			     << " -----> " << std::setw(10) << std::right << var/area*bin
			     << std::endl;
  return var/area*bin;
}

Double_t Eval_func( Double_t fvar, TF1*  func, Bool_t fl_message, Bool_t fl_range, Double_t fxmin, Double_t fxmax ){
  Double_t area = func->Integral( func->GetXmin(), func->GetXmax() );
  if( fl_range ) area = func->Integral( fxmin, fxmax );
  Double_t var  = func->Eval( fvar );
  if( fl_message ) std::cout << "       "
			     << "var  = " << std::setw(10) << std::right << var  
			     << " ("      << std::setw(10) << std::right << fvar << "), "
			     << "area = " << std::setw(10) << std::right << area 
			     << " ----> " << std::setw(10) << std::right << var/area
			     << std::endl;
  return var/area;  
}

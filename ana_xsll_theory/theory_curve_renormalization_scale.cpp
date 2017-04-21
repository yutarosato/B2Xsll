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
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TArrow.h>
#include <TMath.h>
#include <Math/DistFunc.h>
#include <TRandom.h>
#include <TComplex.h>


const Int_t fl_renorm = 0; // 1(sh>0.25 -> mu=2.5GeV), 0(sh<0.25 -> mu=5.0GeV)


const Double_t q2_bound[2][7] = { // [mm,ee]
  {0.0, 4.3,
   //{1.0, 6.0,
   (PDGmass::jpsi -0.25)*(PDGmass::jpsi -0.25), // ~ 8.10
   (PDGmass::jpsi +0.10)*(PDGmass::jpsi +0.10), // ~10.22
   (PDGmass::psi2s-0.15)*(PDGmass::psi2s-0.15), // ~12.50
   (PDGmass::psi2s+0.10)*(PDGmass::psi2s+0.10), // ~14.33
     25.0},
  {0.0,  4.3,
   //{1.0, 6.0,
   (PDGmass::jpsi -0.40)*(PDGmass::jpsi -0.40), // ~ 7.27
   (PDGmass::jpsi +0.15)*(PDGmass::jpsi +0.15), // ~10.54
   (PDGmass::psi2s-0.25)*(PDGmass::psi2s-0.25), // ~11.81
   (PDGmass::psi2s+0.10)*(PDGmass::psi2s+0.10), // ~14.33
   25.0},
};

TComplex GetC7Eff0(Double_t sh, Bool_t nnlo=true) 
{
  // This function returns the zeroth-order alpha_s part of C7
  
  if (!nnlo) return -0.313;
  
  Double_t A7;
  
  // use energy scale of 2.5 GeV as a computational trick (G.Hiller)
  // at least for shat > 0.25
  A7 = -0.353 + 0.023;
  
  TComplex c7eff;
  if (fl_renorm)
    { 
      c7eff = A7;
      return c7eff;
    }
  
  // change energy scale to 5.0 for full NNLO calculation below shat = 0.25
  A7 = -0.312 + 0.008;
  c7eff = A7;
  
  return c7eff;
}



TComplex GetC7Eff1(Double_t sh, Double_t mbeff, Bool_t nnlo=true) 
{
  // This function returns the first-order alpha_s part of C7
  
  if (!nnlo) return 0.0;
  Double_t logsh;
  logsh = log(sh);
  
  TComplex uniti(0.0,1.0);
  
  TComplex c7eff = 0.0;
  if (fl_renorm)
    { 
      return c7eff;
    }
  
  // change energy scale to 5.0 for full NNLO calculation below shat = 0.25
  Double_t muscale = 5.0;
  Double_t alphas = 0.215;
  //Double_t A7 = -0.312 + 0.008;
  Double_t A8 = -0.148;
  //Double_t A9 = 4.174 + (-0.035);
  //Double_t A10 = -4.592 + 0.379;
  Double_t C1 = -0.487;
  Double_t C2 = 1.024;
  //Double_t T9 = 0.374 + 0.252;
  //Double_t U9 = 0.033 + 0.015;
  //Double_t W9 = 0.032 + 0.012;
  Double_t Lmu = log(muscale/mbeff);
  
  TComplex F71;
  TComplex f71;
  TComplex k7100(-0.68192,-0.074998);
  TComplex k7101(0.0,0.0);
  TComplex k7110(-0.23935,-0.12289);
  TComplex k7111(0.0027424,0.019676);
  TComplex k7120(-0.0018555,-0.175);
  TComplex k7121(0.022864,0.011456);
  TComplex k7130(0.28248,-0.12783);
  TComplex k7131(0.029027,-0.0082265);
  f71 = k7100 + k7101*logsh + sh*(k7110 + k7111*logsh) +
    sh*sh*(k7120 + k7121*logsh) + 
    sh*sh*sh*(k7130 + k7131*logsh); 
  F71 = (-208.0/243.0)*Lmu + f71;
  
  TComplex F72;
  TComplex f72;
  TComplex k7200(4.0915,0.44999);
  TComplex k7201(0.0,0.0);
  TComplex k7210(1.4361,0.73732);
  TComplex k7211(-0.016454,-0.11806);
  TComplex k7220(0.011133,1.05);
  TComplex k7221(-0.13718,-0.068733);
  TComplex k7230(-1.6949,0.76698);
  TComplex k7231(-0.17416,0.049359);
  f72 = k7200 + k7201*logsh + sh*(k7210 + k7211*logsh) +
    sh*sh*(k7220 + k7221*logsh) + 
    sh*sh*sh*(k7230 + k7231*logsh); 
  F72 = (416.0/81.0)*Lmu + f72;
  
  TComplex F78;
  F78 = (-32.0/9.0)*Lmu + 8.0*TMath::Pi()*TMath::Pi()/27.0 + (-44.0/9.0) 
    + (-8.0*TMath::Pi()/9.0)*uniti +
    (4.0/3.0*TMath::Pi()*TMath::Pi() - 40.0/3.0)*sh +
    (32.0*TMath::Pi()*TMath::Pi()/9.0 - 316.0/9.0)*sh*sh +
    (200.0*TMath::Pi()*TMath::Pi()/27.0 - 658.0/9.0)*sh*sh*sh +
    (-8.0*logsh/9.0)*(sh + sh*sh + sh*sh*sh);
  
  c7eff = - alphas/(4.0*TMath::Pi())*(C1*F71 + C2*F72 + A8*F78);
  
  return c7eff;
}


TComplex GetC9Eff0(Double_t sh, Double_t mbeff,
		   Bool_t nnlo=true, Bool_t btod=false) 
{
  // This function returns the zeroth-order alpha_s part of C9
  
  if (!nnlo) return 4.344;
  Double_t logsh;
  logsh = log(sh);
  Double_t mch = 0.29;
  
  
  Double_t muscale;
  muscale = 2.5;
  Double_t alphas;
  alphas = 0.267;
  Double_t A8;
  A8 = -0.164;
  Double_t A9;
  A9 = 4.287 + (-0.218);
  Double_t A10;
  A10 = -4.592 + 0.379;
  Double_t C1;
  C1 = -0.697;
  Double_t C2;
  C2 = 1.046;
  Double_t T9;
  T9 = 0.114 + 0.280;
  Double_t U9;
  U9 = 0.045 + 0.023;
  Double_t W9;
  W9 = 0.044 + 0.016;
  
  Double_t Lmu;
  Lmu = log(muscale/mbeff);
  
  
  TComplex uniti(0.0,1.0);
  
  TComplex hc;
  Double_t xarg;
  xarg = 4.0*mch/sh;
  hc = -4.0/9.0*log(mch*mch) + 8.0/27.0 + 4.0*xarg/9.0;
  if (xarg < 1.0)
    { 
      hc = hc - 2.0/9.0*(2.0 + xarg)*sqrt(fabs(1.0 - xarg))*
	(log((sqrt(1.0 - xarg)+1.0)/(sqrt(1.0 - xarg) - 1.0)) - 
	 uniti*TMath::Pi());
    } 
  else
    {
      hc = hc - 2.0/9.0*(2.0 + xarg)*sqrt(fabs(1.0 - xarg))*
	2.0*atan(1.0/sqrt(xarg - 1.0));
    }
  
  TComplex h1;
  xarg = 4.0/sh;
  h1 = 8.0/27.0 + 4.0*xarg/9.0;
  if (xarg < 1.0)
    { 
      h1 = h1 - 2.0/9.0*(2.0 + xarg)*sqrt(fabs(1.0 - xarg))*
	(log((sqrt(1.0 - xarg)+1.0)/(sqrt(1.0 - xarg) - 1.0)) - 
	 uniti*TMath::Pi());
    } 
  else
    {
      h1 = h1 - 2.0/9.0*(2.0 + xarg)*sqrt(fabs(1.0 - xarg))*
	2.0*atan(1.0/sqrt(xarg - 1.0));
    }
  
  TComplex h0;
  h0 = 8.0/27.0 - 4.0*log(2.0)/9.0 + 4.0*uniti*TMath::Pi()/9.0;
  
  
  // X=V_{ud}^* V_ub / V_{td}^* V_tb * (4/3 C_1 +C_2) * (h(\hat m_c^2, hat s)-
  // h(\hat m_u^2, hat s))
  TComplex Vudstar(1.0 - 0.2279*0.2279/2.0, 0.0);
  TComplex Vub((0.118+0.273)/2.0, -1.0*(0.305+0.393)/2.0);
  TComplex Vtdstar(1.0 - (0.118+0.273)/2.0,(0.305+0.393)/2.0);
  TComplex Vtb(1.0,0.0);
  
  TComplex Xd;
  Xd = (Vudstar * Vub / Vtdstar * Vtb) * (4.0/3.0*C1 + C2) * (hc - h0);
  
  TComplex c9eff = 4.344;
  if (fl_renorm)
    { 
      c9eff =  A9 + T9*hc + U9*h1 + W9*h0;
      if (btod)
	{
	  c9eff += Xd; 
	}
      return c9eff;
    }
  
  // change energy scale to 5.0 for full NNLO calculation below shat = 0.25
  muscale = 5.0;
  alphas = 0.215;
  A9 = 4.174 + (-0.035);
  C1 = -0.487;
  C2 = 1.024;
  A8 = -0.148;
  T9 = 0.374 + 0.252;
  U9 = 0.033 + 0.015;
  W9 = 0.032 + 0.012;
  Lmu = log(muscale/mbeff);
  
  Xd = (Vudstar * Vub / Vtdstar * Vtb) * (4.0/3.0*C1 + C2) * (hc - h0);
  
  c9eff = A9 + T9*hc + U9*h1 + W9*h0;
  
  if (btod)
    {
      c9eff += Xd; 
    }
  
  return c9eff;
}


TComplex GetC9Eff1(Double_t sh, Double_t mbeff,
		   Bool_t nnlo=true, Bool_t btod=false) 
{
  // This function returns the first-order alpha_s part of C9
  
  if (!nnlo) return 0.0;
  Double_t logsh;
  logsh = log(sh);
  Double_t mch = 0.29;
  
  TComplex uniti(0.0,1.0);
  
  TComplex c9eff = 0.0;
  if (fl_renorm)
    { 
      return c9eff;
    }
  
  // change energy scale to 5.0 for full NNLO calculation below shat = 0.25
  Double_t muscale = 5.0;
  Double_t alphas = 0.215;
  Double_t C1 = -0.487;
  Double_t C2 = 1.024;
  Double_t A8 = -0.148;
  Double_t Lmu = log(muscale/mbeff);
  
  TComplex F91;
  TComplex f91;
  TComplex k9100(-11.973,0.16371);
  TComplex k9101(-0.081271,-0.059691);
  TComplex k9110(-28.432,-0.25044);
  TComplex k9111(-0.040243,0.016442);
  TComplex k9120(-57.114,-0.86486);
  TComplex k9121(-0.035191,0.027909);
  TComplex k9130(-128.8,-2.5243);
  TComplex k9131(-0.017587,0.050639);
  f91 = k9100 + k9101*logsh + sh*(k9110 + k9111*logsh) +
    sh*sh*(k9120 + k9121*logsh) + 
    sh*sh*sh*(k9130 + k9131*logsh); 
  F91 = (-1424.0/729.0 + 16.0*uniti*TMath::Pi()/243.0 
         + 64.0/27.0*log(mch))*Lmu - 16.0*Lmu*logsh/243.0 +
    (16.0/1215.0 - 32.0/135.0/mch/mch)*Lmu*sh +
    (4.0/2835.0 - 8.0/315.0/mch/mch/mch/mch)*Lmu*sh*sh +
    (16.0/76545.0 - 32.0/8505.0/mch/mch/mch/mch/mch/mch)*
    Lmu*sh*sh*sh -256.0*Lmu*Lmu/243.0 + f91;
  
  TComplex F92;
  TComplex f92;
  TComplex k9200(6.6338,-0.98225);
  TComplex k9201(0.48763,0.35815);
  TComplex k9210(3.3585,1.5026);
  TComplex k9211(0.24146,-0.098649);
  TComplex k9220(-1.1906,5.1892);
  TComplex k9221(0.21115,-0.16745);
  TComplex k9230(-17.12,15.146);
  TComplex k9231(0.10552,-0.30383);
  f92 = k9200 + k9201*logsh + sh*(k9210 + k9211*logsh) +
    sh*sh*(k9220 + k9221*logsh) + 
    sh*sh*sh*(k9230 + k9231*logsh); 
  F92 = (256.0/243.0 - 32.0*uniti*TMath::Pi()/81.0 
         - 128.0/9.0*log(mch))*Lmu + 32.0*Lmu*logsh/81.0 +
    (-32.0/405.0 + 64.0/45.0/mch/mch)*Lmu*sh +
    (-8.0/945.0 + 16.0/105.0/mch/mch/mch/mch)*Lmu*sh*sh +
    (-32.0/25515.0 + 64.0/2835.0/mch/mch/mch/mch/mch/mch)*
    Lmu*sh*sh*sh + 512.0*Lmu*Lmu/81.0 + f92;
  
  TComplex F98;
  F98 = 104.0/9.0 - 32.0*TMath::Pi()*TMath::Pi()/27.0 + 
    (1184.0/27.0 - 40.0*TMath::Pi()*TMath::Pi()/9.0)*sh +
    (14212.0/135.0 - 32.0*TMath::Pi()*TMath::Pi()/3.0)*sh*sh +
    (193444.0/945.0 - 560.0*TMath::Pi()*TMath::Pi()/27.0)*sh*sh*sh +
    16.0*logsh/9.0*(1.0 + sh + sh*sh + sh*sh*sh);
  
  c9eff = - alphas/(4.0*TMath::Pi())*(C1*F91 + C2*F92 + A8*F98);
  
  return c9eff;
}


TComplex GetC10Eff(Double_t sh, Bool_t nnlo=true) 
{
  
  if (!nnlo) return -4.669;
  TComplex A10;
  A10 = -4.592 + 0.379;
  
  TComplex c10eff;
  c10eff = A10;
  
  return c10eff;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TF1* fLi2 = new TF1("fLi2", "-log(1-x)/x", 0, 1 );
Double_t Li2( Double_t x ){
  Double_t f = fLi2->Integral(0,x);
  return f;
};

Double_t Getf710( Double_t sh, Double_t mbeff ){
  Double_t f = 0;
  f +=  6.0*sh*(3.0+ 9.0*sh-2.0*sh*sh)*Li2(sh); // 1st term
  f -= 12.0*sh*(1.0+13.0*sh-4.0*sh*sh)*Li2(sqrt(sh)); // 2nd term
  f += 3.0*(1.0-23.0*sh+23.0*sh*sh-sh*sh*sh)*log(1.0-sh); // 3rd term
  f += 6.0*sh*(13.0-16.0*sh+3.0*sh*sh)*log(1.0-sqrt(sh)); // 4th term
  f += sh*( 5.0*TMath::Pi()*TMath::Pi()*(1.0+sh) - 3.0*(5.0-20.0*sqrt(sh)+sh)*(1.0-sqrt(sh))*(1.0-sqrt(sh)) ); // 5th term
  f += 24.0*sh*(1.0-sh)*(1.0-sh)*log(4.8/mbeff); // 6th term
  f /= -18.0*sh*(1.0-sh)*(1.0-sh);

  //f = 0; // convert to omega

  return f;
};

Double_t GetFuncf710( Double_t* x, Double_t* p ){
  Double_t mbeff  = 4.8;
  Double_t sh = x[0]/(mbeff*mbeff);
  //Double_t z  = x[1];
  return Getf710( sh, mbeff );
}

Double_t Getf910( Double_t sh, Double_t mbeff ){
  Double_t f = 0;

  f +=  6.0*sh*(1.0+ 3.0*sh-sh*sh)*Li2(sh); // 1st term
  f -= 12.0*sh*sh*(5.0-2.0*sh)*Li2(sqrt(sh)); // 2nd term
  f += 3.0*(1.0-10.0*sh+11.0*sh*sh-2.0*sh*sh*sh)*log(1.0-sh); // 3rd term
  f += 6.0*sh*(5.0-7.0*sh+2.0*sh*sh)*log(1.0-sqrt(sh)); // 4th term
  f += sh*( 3.0*(4.0*sqrt(sh)-3.0)*(1.0-sqrt(sh))*(1.0-sqrt(sh)) + TMath::Pi()*TMath::Pi()*(2.0+sh) ); // 5th term
  f /= -9.0*sh*(1.0-sh)*(1.0-sh);
  
  //f = 0; // convert to omega
  return f;
};

Double_t GetFuncf910( Double_t* x, Double_t* p ){
  Double_t mbeff  = 4.8;
  Double_t sh = x[0]/(mbeff*mbeff);
  //Double_t z  = x[1];
  return Getf910( sh, mbeff );
}

Double_t Getf79( Double_t sh, Double_t z, Double_t mbeff ){
  Double_t f = 0;

  f += 3.0*sh*(1.0+sqrt(sh))*(1.0+sqrt(sh))
    *( 3.0*(5.0+z*z) - 3.0*sqrt(sh)*(11.0-z*z) + 16.0*sh ) * Li2(sh);
  f += 12.0*sh*sqrt(sh)*(3.0+sh)*(1.0-3.0*z*z)*Li2(sqrt(sh));
  f += 3.0*(1.0-sh)*(1.0-sh)*( 3.0-z*z+sh*(9.0+z*z) )*log(1.0-sh);
  f += 3.0*sh*sh*( 13.0-15.0*z*z-sh*(5.0+z*z) )*log(sh);
  f += 3.0*sh*(7.0+3.0*z*z+8.0*sh-sqrt(sh)*(17.0-3.0*z*z))
    *(1.0+sqrt(sh))*(1.0+sqrt(sh))*log(1.0-sh)*log(sh);
  f += 6.0*sh*sqrt(sh)*(3.0+sh)*(1.0-3.0*z*z)*log(1.0-sqrt(sh))*log(sh);
  f -= 6.0*sh*(1.0-sh)*(5.0*z*z-sh*(4.0-3.0*z*z));
  f += sh*TMath::Pi()*TMath::Pi()*(7.0+3.0*z*z+8.0*sh*sh-sh*(19.0-9.0*z*z));
  f += 48.0*sh*(1.0-sh)*(1.0-sh)*log(4.8/mbeff);
  f /= -36.0*sh*(1.0-sh)*(1.0-sh);

  // omega79
  /*
  f = -4.0/3.0*log(4.8/mbeff)
    -4.0/3.0*Li2(sh) 
    -2.0/9.0*TMath::Pi()*TMath::Pi()
    -2.0/3.0*log(sh)*log(1.0-sh)
    -1.0/9.0*(2.0+7.0*sh)*log(1.0 - sh)/sh
    -2.0/9.0*sh*(3.0 - 2.0*sh)*log(sh)/pow((1.0 - sh),2) 
    +1.0/18.0*(5.0 - 9.0*sh)/(1.0 - sh);
  */
  return f;
};

Double_t GetFuncf79( Double_t* x, Double_t* p ){
  Double_t mbeff  = 4.8;
  Double_t sh = x[0]/(mbeff*mbeff);
  Double_t z  = x[1];
  return Getf79( sh, z, mbeff );
}

Double_t Getf77( Double_t sh, Double_t z, Double_t mbeff ){
  Double_t f = 0;

  f += 12.0*sqrt(sh)*(3.0+6.0*sh-sh*sh)*(1.0-3.0*z*z)*Li2(sqrt(sh));
  f += 3.0*(1.0+sqrt(sh))*(1.0+sqrt(sh))
    *( 8.0*(1.0+z*z)
       -sqrt(sh)*( 19.0-14.0*sqrt(sh)+15.0*sh-8.0*sh*sqrt(sh)
		   + (7.0-sqrt(sh)*(2.0-sqrt(sh))*(3.0+8.0*sqrt(sh)))*z*z)
       )*Li2(sh);
  f += 6.0*(1.0-sh)*(1.0-sh)*(5.0+z*z+sh*(1.0-z*z))*log(1.0-sh);
  f += 6.0*sh*( 5.0-7.0*z*z+sh*(1.0-11.0*z*z)-2.0*sh*sh*(1.0-z*z) )*log(sh);
  f += 3.0*(1.0+sqrt(sh))*(1.0+sqrt(sh))
    *( 4.0*(1.0+z*z)-sqrt(sh)*( 11.0-z*z-sqrt(sh)*( 6.0-7.0*sqrt(sh)+4.0*sh+(2.0-sqrt(sh))*(3.0+4.0*sqrt(sh))*z*z )))
    *log(1.0-sh)*log(sh);
  f += 6.0*sqrt(sh)*(3.0+6.0*sh-sh*sh)*(1.0-3.0*z*z)*log(1.0-sqrt(sh))*log(sh);
  f += 2.0*( 2.0*TMath::Pi()*TMath::Pi()*(1.0+z*z-3.0*sh*(1.0-z*z)-sh*sh*(1.0-3.0*z*z)+sh*sh*sh*(1.0-z*z))
	   +(1.0-sh)*(sh*(19.0-68.0*z*z)+4.0*(1.0+z*z)-sh*sh*(11.0-16.0*z*z))
	     );
  f += 48.0*(1.0-sh)*(1.0-sh)*(1.0+z*z+sh*(1.0-z*z))*log(4.8/mbeff);

  f /= -18.0*(1.0-sh)*(1.0-sh)*(1.0+z*z+sh*(1.0-z*z));

  // omega77
  /*
  f = -8.0/3.0*log(4.8/mbeff)
    -4.0/3.0*Li2(sh) 
    -2.0/9.0*TMath::Pi()*TMath::Pi()
    -2.0/3.0*log(sh)*log(1.0-sh)
    -log(1.0-sh)*(8.0+sh)/(2.0+sh)/3.0 
    -2.0/3.0*sh*(2.0 - 2.0*sh - sh*sh)*log(sh)/pow((1.0 - sh),2)/(2.0 + sh)
    -(16.0 - 11.0*sh - 17.0*sh*sh)/18.0/(2.0 + sh)/(1.0 - sh);
  */
  return f;
};

Double_t GetFuncf77( Double_t* x, Double_t* p ){
  Double_t mbeff  = 4.8;
  Double_t sh = x[0]/(mbeff*mbeff);
  Double_t z  = x[1];
  return Getf77( sh, z, mbeff );
}

Double_t Getf99( Double_t sh, Double_t z, Double_t mbeff ){
  Double_t f = 0;

  f += 12.0*sqrt(sh)*(5.0+12.0*sh-sh*sh)*(1.0-3.0*z*z)*Li2(sqrt(sh));
  f -= 3.0*(1.0+sqrt(sh))*(1.0+sqrt(sh))
    *( 8.0-11.0*sqrt(sh)+20.0*sh-17.0*sh*sqrt(sh)+8.0*sh*sh-(1.0+sqrt(sh))*(8.0-sqrt(sh)*(9.0-sqrt(sh)*(21.0-8.0*sqrt(sh))))*z*z )
    *Li2(sh);
  f += 6.0*sh*( 3.0-13.0*z*z+sh*(9.0-23.0*z*z)+2.0*sh*sh*(1.0+z*z) )*log(sh);
  f -= 12.0*(1.0-sh)*(1.0-sh)*( 2.0-z*z+sh*(1.0+z*z) )*log(1.0-sh);
  f -= 3.0*(1.0+sqrt(sh))*(1.0+sqrt(sh))
    *( 4.0-4.0*z*z-sqrt(sh)*(3.0+7.0*z*z) + 12.0*sh*(1.0-z*z) -sh*sqrt(sh)*(9.0+5.0*z*z)+4.0*sh*sh*(1.0+z*z) )
    *log(1.0-sh)*log(sh);
  f += 6.0*sqrt(sh)*(5.0+12.0*sh-sh*sh)*(1.0-3.0*z*z)*log(sh)*log(1.0-sqrt(sh));
  f += 3.0*(1.0-sh)*(5.0-5.0*z*z+sh*(28.0-66.0*z*z)-sh*sh*(5.0-3.0*z*z));
  f -= 2.0*TMath::Pi()*TMath::Pi()
    *(2.0-2.0*z*z+5.0*sh*(1.0-3.0*z*z)-sh*sh*(1.0+9.0*z*z)+2.0*sh*sh*sh*(1.0+z*z));
  f /= 18.0*(1.0-sh)*(1.0-sh)*(1.0+sh-z*z*(1.0-sh));

  // omega99
  /*
  f = -2.0/9.0*TMath::Pi()*TMath::Pi() - 4.0/3.0*Li2(sh)
    - 2.0/3.0*log(sh)*log(1.0-sh)
    - (5.0+4.0*sh)/(3.0*(1.0+2.0*sh)) * log(1.0-sh)
    - 2.0*sh*(1.0+sh)*(1.0-2.0*sh)
    /(3.0*pow(1.0-sh,2)*(1.0+2.0*sh)) * log(sh)
    + (5.0+9.0*sh-6.0*sh*sh)/(6.0*(1.0-sh)*(1.0+2.0*sh));
  */
  return f;
};

Double_t GetFuncf99( Double_t* x, Double_t* p ){
  Double_t mbeff  = 4.8;
  Double_t sh = x[0]/(mbeff*mbeff);
  Double_t z  = x[1];
  return Getf99( sh, z, mbeff );
}

Double_t d2gdsdz( const Double_t* x, const Double_t* p){ // x[0] = s(=q^2), x[1] = z(=cos_theta)
  // [PRD 66, 094013 (2002)]
  Bool_t btod = false;
  Bool_t nnlo = true;

  Double_t mb     = 4.8;
  Double_t mbeff  = 4.8;
  Double_t sh     = x[0]/(mbeff*mbeff);
  Double_t z      = x[1];
  //Double_t alphas = 0.119/( 1 + 0.119*log(pow(4.8,2)/pow(91.1867,2))*23.0/12.0/TMath::Pi() ); // ?????
  Double_t        alphas = 0.215; // muscale = 5.0 GeV
  if( fl_renorm ) alphas = 0.267; // muscale = 2.5 GeV
  TComplex c7eff0 = GetC7Eff0( sh,        nnlo       );
  TComplex c7eff1 = GetC7Eff1( sh, mbeff, nnlo       );
  TComplex c9eff0 = GetC9Eff0( sh, mbeff, nnlo, btod );
  TComplex c9eff1 = GetC9Eff1( sh, mbeff, nnlo, btod );
  TComplex c10eff = GetC10Eff( sh,        nnlo       );
 
  TComplex c7eff = c7eff0 + c7eff1;
  TComplex c9eff = c9eff0 + c9eff1;
  Double_t c7c7   = TComplex::Abs(  c7eff ) * TComplex::Abs(  c7eff );
  Double_t c9c9   = TComplex::Abs(  c9eff ) * TComplex::Abs(  c9eff );
  Double_t c10c10 = TComplex::Abs( c10eff ) * TComplex::Abs( c10eff );
  Double_t c7c9   = ( TComplex(c7eff*(TComplex::Conjugate( c9eff))) ).Re();
  Double_t c7c10  = ( TComplex(c7eff*(TComplex::Conjugate(c10eff))) ).Re();
  Double_t c9c10  = ( TComplex(c9eff*(TComplex::Conjugate(c10eff))) ).Re();
  
  Double_t f710   = Getf710  ( sh,    mbeff );
  Double_t f910   = Getf910  ( sh,    mbeff );
  Double_t f79    = Getf79   ( sh, z, mbeff );
  Double_t f77    = Getf77   ( sh, z, mbeff );
  Double_t f99    = Getf99   ( sh, z, mbeff );
  
  Double_t f=0;
  if( p[0]==0 || p[0]==1 ) f += 3.0/4.0*( (1.0-z*z)+sh*(1.0+z*z) )*( c9c9+c10c10 )*( 1.0+2.0*alphas/TMath::Pi()*f99 );
  if( p[0]==0 || p[0]==2 ) f += 3.0/sh*( (1.0+z*z)+sh*(1.0-z*z) )*c7c7*( 1.0+2.0*alphas/TMath::Pi()*f77 );
  if( p[0]==0 || p[0]==3 ) f -= 3.0*sh*z*c9c10*( 1+2*alphas/TMath::Pi()*f910 );
  if( p[0]==0 || p[0]==4 ) f += 6.0*c7c9*(1.0+2.0*alphas/TMath::Pi()*f79);
  if( p[0]==0 || p[0]==5 ) f -= 6.0*z*c7c10*(1.0+2.0*alphas/TMath::Pi()*f710);
  f *= (1-sh)*(1-sh);

  return f;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const Double_t xmin    =  0.0;
const Double_t xmax[2] = {23.0, 21.16}; // {4.8*4.8, 4.6*4.6} = {23.04, 21.16}
//const Double_t xmax[2] = {23.0, 19.8025}; // {4.8*4.8, 4.45*4.45} = {23.04, 19.8025} // for syst. study
const Double_t ymin    = -1.0;
const Double_t ymax    =  1.0;

TGraph* q2_AFB_PRD66094013( Double_t npoint=115 ){
  TGraph* g = new TGraph();
  TF2* fd2gdsdz = new TF2( "fd2gdsdz", d2gdsdz, xmin, xmax[0], ymin, ymax, 1 );
  fd2gdsdz->SetParameter( 0,0 );
  Double_t step = (xmax[0]-xmin)/npoint;
  for( Int_t i=0; i<npoint; i++ ){
    Double_t nf = fd2gdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step,  0.0, 1.0 ); //std::cout << "nf = " << nf << std::endl;
    Double_t nb = fd2gdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step, -1.0, 0.0 ); //std::cout << "nb = " << nb << std::endl;
    Double_t n  = fd2gdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step, -1.0, 1.0 ); //std::cout << "n  = " << n  << std::endl;

    g->SetPoint(i, xmin+(Double_t)(i+0.5)*step, (nf-nb)/n );
  };
  delete fd2gdsdz;
  return g;
}

TGraphErrors* q2_AFB_bin_PRD66094013( Int_t fl_lep, Double_t npoint=115 ){
  TGraphErrors* g = new TGraphErrors();
  TF2* fd2gdsdz = new TF2( "fd2gdsdz", d2gdsdz, xmin, xmax[0], ymin, ymax, 1 );
  fd2gdsdz->SetParameter( 0,0 );
  Double_t step = (xmax[0]-xmin)/npoint;
  Int_t    cnt  = 0;
  Double_t nf   = 0;
  Double_t nb   = 0;
  Double_t n    = 0;
  for( Int_t i=0; i<npoint; i++ ){
    if( i*step<=0.20 ) continue;
    nf += fd2gdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step,  0.0, 1.0 ); //std::cout << "nf = " << nf << std::endl;
    nb += fd2gdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step, -1.0, 0.0 ); //std::cout << "nb = " << nb << std::endl;
    n  += fd2gdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step, -1.0, 1.0 ); //std::cout << "n  = " << n  << std::endl;
    if( (xmin+(Double_t)(i+0.5)*step > q2_bound[fl_lep][cnt+1]) || i==npoint-1 ){
      g->SetPoint     (cnt, (q2_bound[fl_lep][cnt+1]+q2_bound[fl_lep][cnt])/2.0, (nf-nb)/n );
      g->SetPointError(cnt, (q2_bound[fl_lep][cnt+1]-q2_bound[fl_lep][cnt])/2.0, 0.0       );
      cnt++;
      nf = 0;
      nb = 0;
      n  = 0;
    }
  };
  delete fd2gdsdz;
  return g;
}

TGraph* q2_PRD66094013( Double_t npoint=115 ){
  TGraph* g = new TGraph();
  TF2* fd2gdsdz = new TF2( "fd2gdsdz", d2gdsdz, xmin, xmax[0], ymin, ymax, 1 );
  fd2gdsdz->SetParameter( 0,0 );
  Double_t step = (xmax[0]-xmin)/npoint;
  for( Int_t i=0; i<npoint; i++ ){
    Double_t nf = fd2gdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step,  0.0, 1.0 ); //std::cout << "nf = " << nf << std::endl;
    Double_t nb = fd2gdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step, -1.0, 0.0 ); //std::cout << "nb = " << nb << std::endl;
    Double_t n  = fd2gdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step, -1.0, 1.0 ); //std::cout << "n  = " << n  << std::endl;

    if( i!=0 ) g->SetPoint(i, xmin+(Double_t)(i+0.5)*step, n );
  };
  delete fd2gdsdz;
  return g;
}


TGraph* q2_AFB_eff_PRD66094013( Int_t fl_lep, Double_t npoint=70 ){

  const Int_t nybin = 30; // fixed

  TFile file_eff( Form("eff_table/2d_q2_theta_lep%d_xsid0_setA-U.root",fl_lep) );
  TH2D* hist_eff = (TH2D*)file_eff.Get( "h511_h12" );
  TGraph* g = new TGraph();
  TF2* fd2gdsdz = new TF2( "fd2gdsdz", d2gdsdz, xmin, xmax[0], ymin, ymax, 1 );
  fd2gdsdz->SetParameter( 0,0 );
  Double_t stepx = (xmax[0]-xmin)/npoint;
  Double_t stepy = (ymax-ymin)/(Double_t)nybin;
  Int_t    cnt   = 0;
  for( Int_t i=0; i<npoint; i++ ){
    Double_t nf = 0;
    Double_t nb = 0;
    Double_t n  = 0;
    for( Int_t j=0; j<nybin; j++ ){
      Double_t tmp_n   = fd2gdsdz->Integral( xmin+(Double_t)i*stepx, xmin+(Double_t)(i+1)*stepx, -1.0+(Double_t)j*stepy, -1.0+(Double_t)(j+1)*stepy );
      if( i==0 ) tmp_n = fd2gdsdz->Integral( xmin+0.2,               xmin+(Double_t)(i+1)*stepx, -1.0+(Double_t)j*stepy, -1.0+(Double_t)(j+1)*stepy );
      Int_t    tmp_bin = hist_eff->FindBin( xmin+(Double_t)(i+0.5)*stepx, -1.0+(Double_t)(j+0.5)*stepy );
      Double_t tmp_eff = hist_eff->GetBinContent (tmp_bin);
      if( j<nybin/2 ) nb += tmp_n*tmp_eff;
      else            nf += tmp_n*tmp_eff;
      n += tmp_n*tmp_eff;
      //std::cout << "i = " << i << ", j = " << j << " : bin = " << tmp_bin << ", eff = " << tmp_eff << std::endl;
    }
    if( n ) g->SetPoint(cnt++, xmin+(Double_t)(i+0.5)*stepx, (nf-nb)/n );
  }
  
  delete hist_eff;
  file_eff.Close();
  delete fd2gdsdz;
  return g;
}

TGraphErrors* q2_AFB_bin_eff_PRD66094013( Int_t fl_lep, Double_t npoint=115 ){

  const Int_t nybin = 30; // fixed

  TFile file_eff( Form("eff_table/2d_q2_theta_lep%d_xsid0_setA-U.root",fl_lep) );
  TH2D* hist_eff = (TH2D*)file_eff.Get( "h511_h12" );
  TGraphErrors* g = new TGraphErrors();
  TF2* fd2gdsdz = new TF2( "fd2gdsdz", d2gdsdz, xmin, xmax[0], ymin, ymax, 1 );
  fd2gdsdz->SetParameter( 0,0 );
  Double_t stepx = (xmax[0]-xmin)/npoint;
  Double_t stepy = (ymax-ymin)/(Double_t)nybin;
  Int_t    cnt1  = 0; // for boundary counting
  Int_t    cnt2  = 0; // for point    counting
  Double_t nf    = 0;
  Double_t nb    = 0;
  Double_t n     = 0;
  for( Int_t i=0; i<npoint; i++ ){
    if( i*stepx<=0.20 ) continue;
    for( Int_t j=0; j<nybin; j++ ){
      Double_t tmp_n   = fd2gdsdz->Integral( xmin+(Double_t)i*stepx, xmin+(Double_t)(i+1)*stepx, -1.0+(Double_t)j*stepy, -1.0+(Double_t)(j+1)*stepy );
      Int_t    tmp_bin = hist_eff->FindBin( xmin+(Double_t)(i+0.5)*stepx, -1.0+(Double_t)(j+0.5)*stepy );
      Double_t tmp_eff = hist_eff->GetBinContent (tmp_bin);
      if( j<nybin/2 ) nb += tmp_n*tmp_eff;
      else            nf += tmp_n*tmp_eff;
      n += tmp_n*tmp_eff;
      //std::cout << "i = " << i << ", j = " << j << " : bin = " << tmp_bin << ", eff = " << tmp_eff << std::endl;
    }
    if( (xmin+(Double_t)(i+0.5)*stepx > q2_bound[fl_lep][cnt1+1]) || i==npoint-1 ){
      if( cnt1!=2 && cnt1!=4 ){
	g->SetPoint     (cnt2, (q2_bound[fl_lep][cnt1+1]+q2_bound[fl_lep][cnt1])/2.0, (nf-nb)/n );
	g->SetPointError(cnt2, (q2_bound[fl_lep][cnt1+1]-q2_bound[fl_lep][cnt1])/2.0, 0.0       );
	cnt2++;
      }
      cnt1++;
      nf = 0;
      nb = 0;
      n  = 0;
    }
  }
  
  delete hist_eff;
  file_eff.Close();
  delete fd2gdsdz;
  return g;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t GetU( Double_t s, Double_t ml, Double_t mbeff, Double_t mseff ){
  Double_t f = sqrt(
		    (  s-(mbeff+mseff)*(mbeff+mseff) )
		    *( s-(mbeff-mseff)*(mbeff-mseff) )
		    *(1-4.0*ml*ml/s)
		    );
  return f;
}

Double_t GetA1( Double_t s, Double_t z, Double_t ml, Double_t mbeff, Double_t mseff ){
  Double_t f = 0;
  f += -16.0*mbeff*mbeff*mseff*mseff*s*(2*ml*ml*+s);
  f += (mbeff*mbeff+mseff*mseff)
    * (
       ml*ml*(-8*(mbeff*mbeff+mseff*mseff)*s+8.0*(mbeff*mbeff-mseff*mseff)*(mbeff*mbeff-mseff*mseff))
       + 2*s*(-s*s+GetU(s,ml,mbeff,mseff)*GetU(s,ml,mbeff,mseff)*z*z + (mbeff*mbeff-mseff*mseff)*(mbeff*mbeff-mseff*mseff))
       );
  f *= GetU(s,ml,mbeff,mseff)/s/s;
  return f;
}

Double_t GetFuncA1( Double_t* x, Double_t* p ){
  Double_t s     = x[0];
  Double_t z     = x[1];
  Double_t ml    = x[2];
  Double_t mbeff = x[3];
  Double_t mseff = x[4];
  return GetA1(s, z, ml, mbeff, mseff );
}

Double_t GetA2( Double_t s, Double_t z, Double_t ml, Double_t mbeff, Double_t mseff ){
  Double_t f = 0;
  f -= GetU(s,ml,mbeff,mseff)*GetU(s,ml,mbeff,mseff)*z*z;
  f -= s*s;
  f += (mbeff*mbeff-mseff*mseff)*(mbeff*mbeff-mseff*mseff);
  f *= GetU(s,ml,mbeff,mseff);
  return f;
}

Double_t GetFuncA2( Double_t* x, Double_t* p ){
  Double_t s     = x[0];
  Double_t z     = x[1];
  Double_t ml    = x[2];
  Double_t mbeff = x[3];
  Double_t mseff = x[4];
  return GetA2(s, z, ml, mbeff, mseff );
}

Double_t GetA3( Double_t s, Double_t z, Double_t ml, Double_t mbeff, Double_t mseff ){
  Double_t f = 2*GetU(s,ml,mbeff,mseff)*GetU(s,ml,mbeff,mseff)*z*s;
  return f;
}

Double_t GetFuncA3( Double_t* x, Double_t* p ){
  Double_t s     = x[0];
  Double_t z     = x[1];
  Double_t ml    = x[2];
  Double_t mbeff = x[3];
  Double_t mseff = x[4];
  return GetA3(s, z, ml, mbeff, mseff );
}
 
Double_t GetA4( Double_t s, Double_t z, Double_t ml, Double_t mbeff, Double_t mseff ){
  Double_t f = 0;
  f += (mbeff*mbeff+mseff*mseff)*s;
  f -= (mbeff*mbeff-mseff*mseff)*(mbeff*mbeff-mseff*mseff);
  f *= 2*(2*ml*ml+s);
  f *= GetU(s,ml,mbeff,mseff)/s;
  return f;
}

Double_t GetFuncA4( Double_t* x, Double_t* p ){
  Double_t s     = x[0];
  Double_t z     = x[1];
  Double_t ml    = x[2];
  Double_t mbeff = x[3];
  Double_t mseff = x[4];
  return GetA4(s, z, ml, mbeff, mseff );
}

Double_t GetA5( Double_t s, Double_t z, Double_t ml, Double_t mbeff, Double_t mseff ){
  Double_t f = -2*(mbeff*mbeff+mseff*mseff)*GetU(s,ml,mbeff,mseff)*GetU(s,ml,mbeff,mseff)*z;
  return f;
}

Double_t GetFuncA5( Double_t* x, Double_t* p ){
  Double_t s     = x[0];
  Double_t z     = x[1];
  Double_t ml    = x[2];
  Double_t mbeff = x[3];
  Double_t mseff = x[4];
  return GetA5(s, z, ml, mbeff, mseff );
}

Double_t GetL1( Double_t s, Double_t z, Double_t ml, Double_t mbeff, Double_t mseff ){
  Double_t f = -4*GetU(s,ml,mbeff,mseff)*ml*ml*(mbeff*mbeff+mseff*mseff-s);
  return f;
}

Double_t GetFuncL1( Double_t* x, Double_t* p ){
  Double_t s     = x[0];
  Double_t z     = x[1];
  Double_t ml    = x[2];
  Double_t mbeff = x[3];
  Double_t mseff = x[4];
  return GetL1(s, z, ml, mbeff, mseff );
}

Double_t d2bdsdz( const Double_t* x, const Double_t* p){
  // x[0] = s(=q^2), x[1] = z(=cos_theta), x[2] = ml
  // [PRD 59, 074013(1999)]

  Bool_t btod = false;
  Bool_t nnlo = true;
  Double_t mbeff  = 4.8;
  //Double_t mbeff  = 4.95; // for syst. study
  //Double_t mbeff  = 4.65; // for syst. study // * xmax should be adjusted
  Double_t mseff  = 0.2;
  Double_t ml     = 0.0;
  
  Double_t s      = x[0];
  Double_t sh     = x[0]/(mbeff*mbeff);
  Double_t z      = x[1];

  TComplex c7eff0 = GetC7Eff0( sh,        nnlo       );
  TComplex c7eff1 = GetC7Eff1( sh, mbeff, nnlo       );
  TComplex c9eff0 = GetC9Eff0( sh, mbeff, nnlo, btod );
  TComplex c9eff1 = GetC9Eff1( sh, mbeff, nnlo, btod );
  TComplex c10eff = GetC10Eff( sh,        nnlo       );
 
  TComplex c7eff     = c7eff0 + c7eff1;
  TComplex c9eff     = c9eff0 + c9eff1;
  TComplex c9pc10    = c9eff + c10eff; 
  TComplex c9mc10    = c9eff - c10eff;
  Double_t c7_c9pc10 = TComplex(TComplex::Conjugate(c7eff)*c9pc10).Re();
  Double_t c7_c9mc10 = TComplex(TComplex::Conjugate(c7eff)*c9mc10).Re();

  Double_t A1 = GetA1( s, z, ml, mbeff, mseff);
  Double_t A2 = GetA2( s, z, ml, mbeff, mseff);
  Double_t A3 = GetA3( s, z, ml, mbeff, mseff);
  Double_t A4 = GetA4( s, z, ml, mbeff, mseff);
  Double_t A5 = GetA5( s, z, ml, mbeff, mseff);
  Double_t L1 = GetL1( s, z, ml, mbeff, mseff);
  
  Double_t f = 0;
  f += A1*( 4*TComplex::Abs(c7eff)*TComplex::Abs(c7eff) );
  f += A2*( TComplex::Abs(c9mc10)*TComplex::Abs(c9mc10) + TComplex::Abs(c9pc10)*TComplex::Abs(c9pc10) );
  f += A3*( TComplex::Abs(c9mc10)*TComplex::Abs(c9mc10) - TComplex::Abs(c9pc10)*TComplex::Abs(c9pc10) );
  f += A4*( -4*c7_c9mc10 - 4*c7_c9pc10 );
  f += A5*( -4*c7_c9mc10 + 4*c7_c9pc10 );
  f += L1*2*TComplex( c9mc10*TComplex::Conjugate(c9pc10) ).Re();

  return f;
}

TGraph* q2_AFB_PRD59074013( Double_t npoint=210 ){
  TGraph* g = new TGraph();
  TF2* fd2bdsdz = new TF2( "fd2bdsdz", d2bdsdz, xmin, xmax[1], ymin, ymax, 1 );

  Double_t step = (xmax[0]-xmin)/npoint;
  for( Int_t i=0; i<npoint; i++ ){
    if( xmin+(Double_t)(i+1)*step > xmax[1] ) break;
    Double_t nf = fd2bdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step,  0.0, 1.0 ); // std::cout << "nf = " << nf << std::endl;
    Double_t nb = fd2bdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step, -1.0, 0.0 ); // std::cout << "nb = " << nb << std::endl;
    Double_t n  = fd2bdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step, -1.0, 1.0 ); // std::cout << "n  = " << n  << std::endl;

    if( i!=0    ) g   ->SetPoint(i, xmin+(Double_t)(i+0.5)*step, (nf-nb)/n );
  };
  delete fd2bdsdz;
  return g;
}

TGraphErrors* q2_AFB_bin_PRD59074013( Int_t fl_lep, Double_t npoint=210 ){
  TGraphErrors* g = new TGraphErrors();
  TF2* fd2bdsdz = new TF2( "fd2bdsdz", d2bdsdz, xmin, xmax[1], ymin, ymax, 1 );
  Double_t step = (xmax[0]-xmin)/npoint;
  Int_t    cnt  = 0;
  Double_t nf   = 0;
  Double_t nb   = 0;
  Double_t n    = 0;
  for( Int_t i=0; i<npoint; i++ ){
    if( i*step<=0.20 ) continue;
    if( xmin+(Double_t)(i+1)*step > xmax[1] ) break;
    nf += fd2bdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step,  0.0, 1.0 ); //std::cout << "nf = " << nf << std::endl;
    nb += fd2bdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step, -1.0, 0.0 ); //std::cout << "nb = " << nb << std::endl;
    n  += fd2bdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step, -1.0, 1.0 ); //std::cout << "n  = " << n  << std::endl;
    if( (xmin+(Double_t)(i+0.5)*step > q2_bound[fl_lep][cnt+1]) || i==npoint-1 ){
      g->SetPoint     (cnt, (q2_bound[fl_lep][cnt+1]+q2_bound[fl_lep][cnt])/2.0, (nf-nb)/n );
      g->SetPointError(cnt, (q2_bound[fl_lep][cnt+1]-q2_bound[fl_lep][cnt])/2.0, 0.0       );
      std::cout << Form("[Dump(true)] fl_lep=%d] cnt=%d : %f = (%f-%f)/%f", fl_lep, cnt, (nf-nb)/n, nf, nb, n) << std::endl;
      cnt++;
      nf = 0;
      nb = 0;
      n  = 0;
    }
  };
  delete fd2bdsdz;
  return g;
}

TGraph* q2_PRD59074013( Double_t npoint=210 ){
  TGraph* g = new TGraph();
  TF2* fd2bdsdz = new TF2( "fd2bdsdz", d2bdsdz, xmin, xmax[1], ymin, ymax, 1 );
  fd2bdsdz->SetParameter( 0,0 );
  Double_t step = (xmax[0]-xmin)/npoint;
  for( Int_t i=0; i<npoint; i++ ){
    if( xmin+(Double_t)(i+1)*step > xmax[1] ) break;
    Double_t nf = fd2bdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step,  0.0, 1.0 ); // std::cout << "nf = " << nf << std::endl;
    Double_t nb = fd2bdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step, -1.0, 0.0 ); // std::cout << "nb = " << nb << std::endl;
    Double_t n  = fd2bdsdz->Integral( xmin+(Double_t)i*step, xmin+(Double_t)(i+1)*step, -1.0, 1.0 ); // std::cout << "n  = " << n  << std::endl;
    if( i!=0 ) g->SetPoint(i, xmin+(Double_t)(i+0.5)*step, n         );
  };
  delete fd2bdsdz;
  return g;
}

TGraph* q2_AFB_eff_PRD59074013( Int_t fl_lep, Double_t npoint=70 ){
  const Int_t nybin = 30; // fixed

  TFile file_eff_lep1( "eff_table/2d_q2_theta_lep1_xsid0_setA-U.root" );
  TFile file_eff_lep0( "eff_table/2d_q2_theta_lep0_xsid0_setA-U.root" );
  TH2D* hist_eff_lep1 = (TH2D*)file_eff_lep1.Get( "h511_h12" );
  TH2D* hist_eff_lep0 = (TH2D*)file_eff_lep0.Get( "h511_h12" );
  TGraph* g = new TGraph();
  TF2* fd2bdsdz = new TF2( "fd2bdsdz", d2bdsdz, xmin, xmax[1], ymin, ymax, 1 );
  fd2bdsdz->SetParameter( 0,0 );
  Double_t stepx = (xmax[0]-xmin)/npoint;
  Double_t stepy = (ymax-ymin)/(Double_t)nybin;
  Int_t    cnt   = 0;
  for( Int_t i=0; i<npoint; i++ ){
    Double_t nf = 0;
    Double_t nb = 0;
    Double_t n  = 0;
    if( xmin+(Double_t)(i+1)*stepx > xmax[1] ) break;
    for( Int_t j=0; j<nybin; j++ ){
      Double_t tmp_n   = fd2bdsdz->Integral( xmin+(Double_t)i*stepx, xmin+(Double_t)(i+1)*stepx, -1.0+(Double_t)j*stepy, -1.0+(Double_t)(j+1)*stepy );
      if( i==0 ) tmp_n = fd2bdsdz->Integral( xmin+0.2,               xmin+(Double_t)(i+1)*stepx, -1.0+(Double_t)j*stepy, -1.0+(Double_t)(j+1)*stepy );
      Int_t    tmp_bin_lep1 = hist_eff_lep1->FindBin( xmin+(Double_t)(i+0.5)*stepx, -1.0+(Double_t)(j+0.5)*stepy );
      Int_t    tmp_bin_lep0 = hist_eff_lep0->FindBin( xmin+(Double_t)(i+0.5)*stepx, -1.0+(Double_t)(j+0.5)*stepy );
      Double_t tmp_eff_lep1 = hist_eff_lep1->GetBinContent (tmp_bin_lep1);
      Double_t tmp_eff_lep0 = hist_eff_lep0->GetBinContent (tmp_bin_lep0);
      if( j<nybin/2 ){
	if     ( fl_lep==1 ) nb += tmp_n*tmp_eff_lep1;
	else if( fl_lep==0 ) nb += tmp_n*tmp_eff_lep0;
	else                 nb += tmp_n*(tmp_eff_lep1+tmp_eff_lep0);
      }else{
	if     ( fl_lep==1 ) nf += tmp_n*tmp_eff_lep1;
	else if( fl_lep==0 ) nf += tmp_n*tmp_eff_lep0;
	else                 nf += tmp_n*(tmp_eff_lep1+tmp_eff_lep0);
      }
      if     ( fl_lep==1 ) n += tmp_n*tmp_eff_lep1;
      else if( fl_lep==0 ) n += tmp_n*tmp_eff_lep0;
      else                 n += tmp_n*(tmp_eff_lep1+tmp_eff_lep0);
    }
    if( n ) g->SetPoint(cnt++, xmin+(Double_t)(i+0.5)*stepx, (nf-nb)/n );
  }

  delete hist_eff_lep1;
  delete hist_eff_lep0;
  file_eff_lep1.Close();
  file_eff_lep0.Close();
  delete fd2bdsdz;
  return g;
}

TGraphErrors* q2_AFB_bin_eff_PRD59074013( Int_t fl_lep, Double_t npoint=210 ){

  const Int_t nybin = 30; // fixed

  TFile file_eff( Form("eff_table/2d_q2_theta_lep%d_xsid0_setA-U.root",fl_lep) );
  TH2D* hist_eff = (TH2D*)file_eff.Get( "h511_h12" );
  TGraphErrors* g = new TGraphErrors();
  TF2* fd2bdsdz = new TF2( "fd2bdsdz", d2bdsdz, xmin, xmax[1], ymin, ymax, 1 );
  fd2bdsdz->SetParameter( 0,0 );
  Double_t stepx = (xmax[0]-xmin)/npoint;
  Double_t stepy = (ymax-ymin)/(Double_t)nybin;
  Int_t    cnt1  = 0; // for boundary counting
  Int_t    cnt2  = 0; // for point    counting
  Double_t nf    = 0;
  Double_t nb    = 0;
  Double_t n     = 0;
  for( Int_t i=0; i<npoint; i++ ){
    if( xmin+(Double_t)(i+1)*stepx > xmax[1] ) break;
    if( i*stepx<=0.20 ) continue;
    for( Int_t j=0; j<nybin; j++ ){
      Double_t tmp_n   = fd2bdsdz->Integral( xmin+(Double_t)i*stepx, xmin+(Double_t)(i+1)*stepx, -1.0+(Double_t)j*stepy, -1.0+(Double_t)(j+1)*stepy );
      Int_t    tmp_bin = hist_eff->FindBin( xmin+(Double_t)(i+0.5)*stepx, -1.0+(Double_t)(j+0.5)*stepy );
      Double_t tmp_eff = hist_eff->GetBinContent (tmp_bin);
      if( j<nybin/2 ) nb += tmp_n*tmp_eff;
      else            nf += tmp_n*tmp_eff;
      n += tmp_n*tmp_eff;
      //std::cout << "i = " << i << ", j = " << j << " : bin = " << tmp_bin << ", eff = " << tmp_eff << std::endl;
    }
    if( (xmin+(Double_t)(i+0.5)*stepx > q2_bound[fl_lep][cnt1+1]) || i==npoint-1 ){
      if( cnt1!=2 && cnt1!=4 ){
	g->SetPoint     (cnt2, (q2_bound[fl_lep][cnt1+1]+q2_bound[fl_lep][cnt1])/2.0, (nf-nb)/n );
	g->SetPointError(cnt2, (q2_bound[fl_lep][cnt1+1]-q2_bound[fl_lep][cnt1])/2.0, 0.0       );
	std::cout << Form("[Dump(obs)] fl_lep=%d] cnt=%d : %f = (%f-%f)/%f", fl_lep, cnt2, (nf-nb)/n, nf, nb, n) << std::endl;
	cnt2++;
      }
      cnt1++;
      nf = 0;
      nb = 0;
      n  = 0;
    }
  }
  
  delete hist_eff;
  file_eff.Close();
  delete fd2bdsdz;
  return g;
}


 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const Bool_t flag_save = true; // outfile.eps and outfile.root
Int_t main( Int_t argc, Char_t** argv ){
  
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::cout << "start!" << std::endl;

  TCanvas* c1 = Canvas( "c1", "c1", 5, 3 );
  
  TF2* f99  = new TF2( "f99",  GetFuncf99,  xmin, xmax[0], ymin, ymax, 0 ); f99 ->SetTitle("f99" ); c1->cd(1); f99 ->Draw("colz");  
  TF2* f77  = new TF2( "f77",  GetFuncf77,  xmin, xmax[0], ymin, ymax, 0 ); f77 ->SetTitle("f77" ); c1->cd(2); f77 ->Draw("colz");
  TF2* f910 = new TF2( "f910", GetFuncf910, xmin, xmax[0], ymin, ymax, 0 ); f910->SetTitle("f910"); c1->cd(3); f910->Draw("colz");
  TF2* f79  = new TF2( "f79",  GetFuncf79,  xmin, xmax[0], ymin, ymax, 0 ); f79 ->SetTitle("f79" ); c1->cd(4); f79 ->Draw("colz");
  TF2* f710 = new TF2( "f710", GetFuncf710, xmin, xmax[0], ymin, ymax, 0 ); f710->SetTitle("f710"); c1->cd(5); f710->Draw("colz");

  TF2* f1_1 = new TF2( "f1_1", d2gdsdz, xmin, xmax[0], ymin, ymax, 1 ); f1_1->SetTitle( "q^{2} - cos#theta (1st term:c9,c10)" ); f1_1->SetParameter( 0,1 ); c1->cd( 6); f1_1->Draw("colz");
  TF2* f1_2 = new TF2( "f1_2", d2gdsdz, xmin, xmax[0], ymin, ymax, 1 ); f1_2->SetTitle( "q^{2} - cos#theta (2nd term:c7)"     ); f1_2->SetParameter( 0,2 ); c1->cd( 7); f1_2->Draw("colz");
  TF2* f1_3 = new TF2( "f1_3", d2gdsdz, xmin, xmax[0], ymin, ymax, 1 ); f1_3->SetTitle( "q^{2} - cos#theta (3rd term:c9c10)"  ); f1_3->SetParameter( 0,3 ); c1->cd( 8); f1_3->Draw("colz");
  TF2* f1_4 = new TF2( "f1_4", d2gdsdz, xmin, xmax[0], ymin, ymax, 1 ); f1_4->SetTitle( "q^{2} - cos#theta (4th term:c7c9)"   ); f1_4->SetParameter( 0,4 ); c1->cd( 9); f1_4->Draw("colz");
  TF2* f1_5 = new TF2( "f1_5", d2gdsdz, xmin, xmax[0], ymin, ymax, 1 ); f1_5->SetTitle( "q^{2} - cos#theta (5th term:c7c10)"  ); f1_5->SetParameter( 0,5 ); c1->cd(10); f1_5->Draw("colz");
  TF2* f1   = new TF2( "f1",   d2gdsdz, xmin, xmax[0], ymin, ymax, 1 ); f1  ->SetTitle( "q^{2} - cos#theta"                   ); f1  ->SetParameter( 0,0 ); c1->cd(11); f1  ->Draw("colz"); // surf3
  c1->cd(12); fLi2->Draw();
  c1->cd(13);

  TGraph*       g1_afb_true          = q2_AFB_PRD66094013        ();  g1_afb_true         ->SetLineColor(1);
  TGraphErrors* g1_afb_bin_true_lep1 = q2_AFB_bin_PRD66094013    (1); g1_afb_bin_true_lep1->SetLineColor(2); g1_afb_bin_true_lep1->SetMarkerColor(2);
  TGraphErrors* g1_afb_bin_true_lep0 = q2_AFB_bin_PRD66094013    (0); g1_afb_bin_true_lep0->SetLineColor(3); g1_afb_bin_true_lep0->SetMarkerColor(3);
  TGraph*       g1_afb_obs_lep1      = q2_AFB_eff_PRD66094013    (1); g1_afb_obs_lep1     ->SetLineColor(2);
  TGraph*       g1_afb_obs_lep0      = q2_AFB_eff_PRD66094013    (0); g1_afb_obs_lep0     ->SetLineColor(3);
  TGraphErrors* g1_afb_bin_obs_lep1  = q2_AFB_bin_eff_PRD66094013(1); g1_afb_bin_obs_lep1 ->SetLineColor(2); g1_afb_bin_obs_lep1->SetMarkerColor(2); g1_afb_bin_obs_lep1->SetMarkerStyle(24);
  TGraphErrors* g1_afb_bin_obs_lep0  = q2_AFB_bin_eff_PRD66094013(0); g1_afb_bin_obs_lep0 ->SetLineColor(3); g1_afb_bin_obs_lep0->SetMarkerColor(3); g1_afb_bin_obs_lep0->SetMarkerStyle(24);
  g1_afb_true         ->Draw("AC"   );
  g1_afb_bin_true_lep1->Draw("Psame");
  g1_afb_bin_true_lep0->Draw("Psame");
  g1_afb_obs_lep1     ->Draw("Csame");
  g1_afb_obs_lep0     ->Draw("Csame");
  g1_afb_bin_obs_lep1 ->Draw("Psame");
  g1_afb_bin_obs_lep0 ->Draw("Psame");
  c1->cd(14);
  TGraph* g1_q2 = q2_PRD66094013( 100 );
  g1_q2->Draw("AC");

  c1->Print( Form("pic/theory_curve_renormalization_%d_c1.eps",fl_renorm) );  


  TCanvas* c2 = Canvas( "c2", "c2", 2, 2 );

  TF2* f2 = new TF2( "f2",  d2bdsdz, xmin, xmax[1], ymin, ymax, 0 );
  f2->SetTitle( "q^{2} - cos#theta" );
  c2->cd(1);
  f2->Draw("colz"); // surf3

  c2->cd(2);
  TGraph* g2_afb_true          = q2_AFB_PRD59074013        ();  g2_afb_true         ->SetName( "afb_true"          ); g2_afb_true         ->SetLineColor(1);
  TGraph* g2_afb_bin_true_lep1 = q2_AFB_bin_PRD59074013    (1); g2_afb_bin_true_lep1->SetName( "afb_bin_true_lep1" ); g2_afb_bin_true_lep1->SetLineColor(2); g2_afb_bin_true_lep1->SetMarkerColor(2);
  TGraph* g2_afb_bin_true_lep0 = q2_AFB_bin_PRD59074013    (0); g2_afb_bin_true_lep0->SetName( "afb_bin_true_lep0" ); g2_afb_bin_true_lep0->SetLineColor(3); g2_afb_bin_true_lep0->SetMarkerColor(3);
  TGraph* g2_afb_obs           = q2_AFB_eff_PRD59074013    (2); g2_afb_obs          ->SetName( "afb_obs"           ); g2_afb_obs          ->SetLineColor(1);
  TGraph* g2_afb_obs_lep1      = q2_AFB_eff_PRD59074013    (1); g2_afb_obs_lep1     ->SetName( "afb_obs_lep1"      ); g2_afb_obs_lep1     ->SetLineColor(2);
  TGraph* g2_afb_obs_lep0      = q2_AFB_eff_PRD59074013    (0); g2_afb_obs_lep0     ->SetName( "afb_obs_lep0"      ); g2_afb_obs_lep0     ->SetLineColor(3);
  TGraph* g2_afb_bin_obs_lep1  = q2_AFB_bin_eff_PRD59074013(1); g2_afb_bin_obs_lep1 ->SetName( "afb_bin_obs_lep1"  ); g2_afb_bin_obs_lep1 ->SetLineColor(2); g2_afb_bin_obs_lep1->SetMarkerColor(2); g2_afb_bin_obs_lep1->SetMarkerStyle(24);
  TGraph* g2_afb_bin_obs_lep0  = q2_AFB_bin_eff_PRD59074013(0); g2_afb_bin_obs_lep0 ->SetName( "afb_bin_obs_lep0"  ); g2_afb_bin_obs_lep0 ->SetLineColor(3); g2_afb_bin_obs_lep0->SetMarkerColor(3); g2_afb_bin_obs_lep0->SetMarkerStyle(24);

  g2_afb_true         ->Draw("AC");
  g2_afb_bin_true_lep1->Draw("Psame");
  g2_afb_bin_true_lep0->Draw("Psame");
  g2_afb_obs          ->Draw("Csame");
  g2_afb_obs_lep1     ->Draw("Csame");
  g2_afb_obs_lep0     ->Draw("Csame");
  g2_afb_bin_obs_lep1 ->Draw("Psame");
  g2_afb_bin_obs_lep0 ->Draw("Psame");
  
  c2->cd(3);
  TGraph* g2_q2 = q2_PRD59074013( 100 );
  g2_q2->Draw("AC");

  c2->Print( Form("pic/theory_curve_renormalization_%d_c2.eps",fl_renorm) );  

  TCanvas* c3 = Canvas( "c3", "c3", 2, 1 );
  c3->cd(1);
  g1_afb_true         ->Draw("AC"   );
  g1_afb_bin_true_lep1->Draw("Psame");
  g1_afb_bin_true_lep0->Draw("Psame");
  g1_afb_obs_lep1     ->Draw("Csame");
  g1_afb_obs_lep0     ->Draw("Csame");
  g1_afb_bin_obs_lep1 ->Draw("Psame");
  g1_afb_bin_obs_lep0 ->Draw("Psame");
  c3->cd(2);
  g2_afb_true         ->Draw("AC");
  g2_afb_bin_true_lep1->Draw("Psame");
  g2_afb_bin_true_lep0->Draw("Psame");
  g2_afb_obs          ->Draw("Csame");
  g2_afb_obs_lep1     ->Draw("Csame");
  g2_afb_obs_lep0     ->Draw("Csame");
  g2_afb_bin_obs_lep1 ->Draw("Psame");
  g2_afb_bin_obs_lep0 ->Draw("Psame");
  c3->Print( Form("pic/theory_curve_renormalization_%d_c3.eps",fl_renorm) );  

  TFile outfile( Form("pic/theory_curve_renormalization_%d.root",fl_renorm), "RECREATE" );
  c1->Write();
  c2->Write();
  g2_afb_true         ->Write();
  g2_afb_obs          ->Write();
  g2_afb_obs_lep1     ->Write();
  g2_afb_obs_lep0     ->Write();
  g2_afb_bin_true_lep1->Write();
  g2_afb_bin_true_lep0->Write();
  g2_afb_bin_obs_lep1 ->Write();
  g2_afb_bin_obs_lep0 ->Write();
  std::cout << "AFB(true,ee)" << std::endl; g2_afb_bin_true_lep1->Print();
  std::cout << "AFB(true,mm)" << std::endl; g2_afb_bin_true_lep0->Print();
  std::cout << "AFB(obs, ee)" << std::endl; g2_afb_bin_obs_lep1 ->Print();
  std::cout << "AFB(obs, mm)" << std::endl; g2_afb_bin_obs_lep0 ->Print();
  outfile.Close();

  std::cout << "finish!" << std::endl;
  app.Run();
  
  return 0;
}

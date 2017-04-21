#ifndef FORM_FACTOR_H
#define FORM_FACTOR_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <map>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Set/Nominal_cut_selection.h"
#include "../Set/Branch.h"


#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TCut.h>
#include <TFile.h>


// x = q2
// p[0] = r1, p[1] = m_R^2, p[2] = r2, p[3] = m_fit^2
Double_t func_ball_59(const Double_t* x, const Double_t* p){
  Double_t r1    = p[0];
  Double_t mR2   = p[1];
  Double_t r2    = p[2];
  Double_t mfit2 = p[3];
  Double_t f = r1/(1-x[0]/mR2) + r2/(1-x[0]/mfit2);
  return f;
}

// x = q2
// p[0] = r1, p[1] = r2, p[2] = m_fit^2
Double_t func_ball_60(const Double_t* x, const Double_t* p){
  Double_t r1    = p[0];
  Double_t r2    = p[1];
  Double_t mfit2 = p[2];
  Double_t f = r1/(1-x[0]/mfit2) + r2/(1-x[0]/mfit2)/(1-x[0]/mfit2);
  return f;
}

// x = q2
// p[0] = r1, p[1] = r2, p[2] = m_fit^2, p[3] = m_fit^2'
Double_t func_ball_60_mod(const Double_t* x, const Double_t* p){
  Double_t r1      = p[0];
  Double_t r2      = p[1];
  Double_t mfit2_1 = p[2];
  Double_t mfit2_2 = p[3];
  Double_t f = r1/(1-x[0]/mfit2_1) + r2/(1-x[0]/mfit2_2)/(1-x[0]/mfit2_2);
  return f;
}

// x = q2
// p[0] = r2, p[1] = m_fit^2
Double_t func_ball_61(const Double_t* x, const Double_t* p){
  Double_t r2    = p[0];
  Double_t mfit2 = p[1];
  Double_t f = r2/(1-x[0]/mfit2);
  return f;
}

// x = q2
// p[0] = r1, p[1] = r2, p[2] = m_fit^2 <- T3~
// p[3] = r2, p[4] = m_fit^2            <- T2
Double_t func_ball_T3(const Double_t* x, const Double_t* p){
  Double_t mB = 5.280;
  Double_t mV = 0.892;
  Double_t f  = ( func_ball_60(x,&p[0]) - func_ball_61(x,&p[3]) );
  f *= (mB*mB - mV*mV)/x[0];
  return f;
}

// x = q2
// p[0] = r1, p[1] = r2, p[2] = m_fit^2_1, p[3] = m_fit^2_2 <- T3~
// p[4] = r2, p[5] = m_fit^2                                <- T2
Double_t func_ball_T3_mod(const Double_t* x, const Double_t* p){
  Double_t mB = 5.280;
  Double_t mV = 0.892;
  Double_t f  = ( func_ball_60_mod(x,&p[0]) - func_ball_61(x,&p[4]) );
  f *= (mB*mB - mV*mV)/x[0];
  return f;
}

// x = q2
// p[0] = r2, p[1] = m_fit^2
// p[2] = r1, p[3] = r2, p[4] = m_fit^2
Double_t func_ball_A3(const Double_t* x, const Double_t* p){
  Double_t mB = 5.280;
  Double_t mV = 0.892;
  Double_t f  = (mB + mV)/2/mV * func_ball_61(x,&p[0]) - 
                (mB - mV)/2/mV * func_ball_60(x,&p[2]);

  return f;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// x = q2
// p[0] = F(0), p[1] = c1, p[2] = c2, p[3] = c3
Double_t func_ali(const Double_t* x, const Double_t* p){
  Double_t F0   = p[0];
  Double_t c1   = p[1];
  Double_t c2   = p[2];
  Double_t c3   = p[3];
  Double_t mB   = 5.280;
  Double_t shat = x[0]/mB/mB;
  Double_t f = F0*exp( c1*shat + c2*shat*shat + c3*shat*shat*shat );
  return f;
}

// x = q2
// p[0] = F(0), p[1] = c1, p[2] = c2, p[3] = c3
// p[4] = F(0), p[5] = c1, p[6] = c2, p[7] = c3
Double_t func_ali_A3(const Double_t* x, const Double_t* p){
  Double_t F0   = p[0];
  Double_t c1   = p[1];
  Double_t c2   = p[2];
  Double_t c3   = p[3];
  Double_t mB   = 5.280;
  Double_t mV   = 0.892;
  Double_t f  = (mB+mV)/2.0/mV * func_ali(x,&p[0]) -
                (mB-mV)/2.0/mV * func_ali(x,&p[4]);
  return f;
}

#endif

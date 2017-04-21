#ifndef Nominal_H
#define Nominal_H

#include "Branch.h"
#include "../Util/MChain.h"
#include <TROOT.h>
#include <TSystem.h>
#include <string.h>

Nominal_Cut nominal_cut_selection( MChain* chain, int fl_mode );
// 1(e), 0(mu) <- only de-Mbc cut
//  110 : double lepton veto(e+e-) 
//  111 : single lepton veto(e+)
// -111 : single lepton veto(e-)
//  100 : double lepton veto(m+m-) 
//  101 : single lepton veto(m+)
// -101 : single lepton veto(m-)
#endif

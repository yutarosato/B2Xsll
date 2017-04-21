#ifndef BCS_LR_BKG_H
#define BCS_LR_BKG_H

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

#include "draws_.h"

#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TCut.h>
#include <TFile.h>

typedef std::pair    < double, int > type_di;
typedef std::multimap< double, int >::iterator type_MMapIter;

Int_t main  ( Int_t argc, Char_t** argv );
Int_t select( std::multimap<double, int>& event );

#endif

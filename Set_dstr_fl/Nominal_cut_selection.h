#ifndef Nominal_H
#define Nominal_H

#include "Branch.h"
#include "../Util/MChain.h"
#include <TROOT.h>
#include <TSystem.h>
#include <string.h>

Nominal_Cut nominal_cut_selection( MChain* chain, int fl_mode=0 );

#endif

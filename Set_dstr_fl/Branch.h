#ifndef Branch_H
#define Branch_H

#include "../Util/const.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include <TROOT.h>
#include <vector>
#include <iostream>

typedef void(*Nominal_Cut)(MCut_array*, const Char_t*);

MCut_array branch_table();

void make_cut_common ( MCut_array* cut, const Char_t* h );

#endif

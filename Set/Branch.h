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

void make_cut_common        ( MCut_array* cut, const Char_t* h );
void make_cut_ee            ( MCut_array* cut, const Char_t* h );
void make_cut_mm            ( MCut_array* cut, const Char_t* h );
void make_cut_emu           ( MCut_array* cut, const Char_t* h );
void make_cut_ee_dembc      ( MCut_array* cut, const Char_t* h );
void make_cut_mm_dembc      ( MCut_array* cut, const Char_t* h );
void make_cut_emu_dembc     ( MCut_array* cut, const Char_t* h );
void make_cut_jpsi_ee       ( MCut_array* cut, const Char_t* h );
void make_cut_jpsi_mm       ( MCut_array* cut, const Char_t* h );
void make_cut_jpsi_ee_dembc ( MCut_array* cut, const Char_t* h );
void make_cut_jpsi_mm_dembc ( MCut_array* cut, const Char_t* h );

void make_cut_cc_ee  ( MCut_array* cut );
void make_cut_cc_mm  ( MCut_array* cut );
void make_cut_QQ     ( MCut_array* cut );

#endif

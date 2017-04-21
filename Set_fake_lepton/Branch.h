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

void make_cut_common   ( MCut_array* cut, const Char_t* h );
void make_cut_ee_dembc ( MCut_array* cut, const Char_t* h ); // de&mbc cut for ee
void make_cut_mm_dembc ( MCut_array* cut, const Char_t* h ); // de&mbc cut for mm
void make_cut_double_ee( MCut_array* cut, const Char_t* h ); // double lepton veto
void make_cut_double_mm( MCut_array* cut, const Char_t* h ); // double lepton veto
void make_cut_single_ep( MCut_array* cut, const Char_t* h ); // signle lepton veto( e+ veto)
void make_cut_single_mp( MCut_array* cut, const Char_t* h ); // single lepton veto(mu+ veto)
void make_cut_single_em( MCut_array* cut, const Char_t* h ); // signle lepton veto( e- veto)
void make_cut_single_mm( MCut_array* cut, const Char_t* h ); // single lepton veto(mu- veto)

void make_cut_cc_ee  ( MCut_array* cut );
void make_cut_cc_mm  ( MCut_array* cut );

#endif

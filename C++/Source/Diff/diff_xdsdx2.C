/*
 *  Methods for the class Diff_xdsdx2
 *
 *    (see file diff.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

char diff_xdsdx2_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2005/01/10 16:34:52  j_novak
 * New class for 1D mono-domain differential operators.
 *
 *
 * $Header$
 *
 */

// C headers
#include <assert.h>
#include <stdlib.h>

// Lorene headers
#include "diff.h"
#include "proto.h"

namespace {
    int nap = 0 ;
    Matrice* tab[MAX_BASE*Diff::max_points] ;
    int nr_done[Diff::max_points] ;
}

Diff_xdsdx2::Diff_xdsdx2(int base_r, int nr) : Diff(base_r, nr) {
    initialize() ;
    assert ((base != R_CHEBP)&&(base != R_CHEBI)) ;
}

Diff_xdsdx2::Diff_xdsdx2(const Diff_xdsdx2& diff_in) : Diff(diff_in) {
    assert (nap != 0) ;
    assert ((base != R_CHEBP)&&(base != R_CHEBI)) ;
} 

Diff_xdsdx2::~Diff_xdsdx2() {}

void Diff_xdsdx2::initialize() {
    if (nap == 0) {
	for (int i=0; i<max_points; i++) {
	    nr_done[i] = -1 ;
	    for (int j=0; j<MAX_BASE; j++) 
		tab[j*max_points+i] = 0x0 ;
	}
	nap = 1 ;
    }
    return ;
}

void Diff_xdsdx2::operator=(const Diff_xdsdx2& diff_in) {
    assert (nap != 0) ;
    Diff::operator=(diff_in) ;

}

const Matrice& Diff_xdsdx2::get_matrice() const {
    
    bool done = false ;
    int indice ;
    for (indice =0; indice<max_points; indice++) {
	if (nr_done[indice] == npoints) {
	    if (tab[base*max_points + indice] != 0x0) done = true ;
	    break ;
	}
	if (nr_done[indice] == -1)
	    break ;
    }
    if (!done) { //The computation must be done ...
	if (indice == max_points) {
	    cerr << "Diff_xdsdx2::get_matrice() : no space left!!" << '\n'
		 << "The value of Diff.max_points must be increased..." << endl ;
	    abort() ;
	}
	nr_done[indice] = npoints ;
	tab[base*max_points + indice] = new Matrice(npoints, npoints) ;
	Matrice& resu = *tab[base*max_points + indice] ;
	resu.set_etat_qcq() ;

	double* vect = new double[npoints] ;
	double* cres = new double[npoints] ;
	for (int i=0; i<npoints; i++) {
	    for (int j=0; j<npoints; j++)
		vect[j] = 0. ;
	    vect[i] = 1. ;
	    d2sdx2_1d(npoints, &vect, base << TRA_R) ;
	    if (base == R_CHEBU) {
		mult_xm1_1d_cheb(npoints, vect, cres) ;
		for (int j=0; j<npoints; j++)
		    resu.set(j,i) = cres[j] ;
	    }
	    else {
		multx_1d(npoints, &vect, base << TRA_R) ;
		for (int j=0; j<npoints; j++)
		    resu.set(j,i) = vect[j] ;
	    }
	}
	delete [] vect ;
	delete [] cres ;
    }
    return *tab[base*max_points + indice] ;
}

ostream& Diff_xdsdx2::operator>>(ostream& ost) const {

    ost << " x d2 / dx2 " << endl ;

    return ost ;

}

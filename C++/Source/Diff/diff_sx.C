/*
 *  Methods for the class Diff_sx
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.6  2016/12/05 16:17:50  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:52:51  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:05  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2007/12/11 15:28:11  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.2  2006/04/10 15:20:52  j_novak
 * Operators dsdx and sx can now be used in the nucleus.
 *
 * Revision 1.1  2005/01/10 16:34:52  j_novak
 * New class for 1D mono-domain differential operators.
 *
 *
 * $Header$
 *
 */

// C headers
#include <cassert>
#include <cstdlib>

// Lorene headers
#include "diff.h"
#include "proto.h"

namespace Lorene {
void sxpun_1d(int, double**, int) ;

namespace {
    int nap = 0 ;
    Matrice* tab[MAX_BASE*Diff::max_points] ;
    int nr_done[Diff::max_points] ;
}

Diff_sx::Diff_sx(int base_r, int nr) : Diff(base_r, nr) {
    initialize() ;
    assert(base != R_CHEB) ;
}

Diff_sx::Diff_sx(const Diff_sx& diff_in) : Diff(diff_in) {
    assert (nap != 0) ;
    assert(base != R_CHEB) ;
} 

Diff_sx::~Diff_sx() {}

void Diff_sx::initialize() {
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

void Diff_sx::operator=(const Diff_sx& diff_in) {
    assert (nap != 0) ;
    Diff::operator=(diff_in) ;

}

const Matrice& Diff_sx::get_matrice() const {
    
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
	    cerr << "Diff_sx::get_matrice() : no space left!!" << '\n'
		 << "The value of Diff.max_points must be increased..." << endl ;
	    abort() ;
	}
	nr_done[indice] = npoints ;
	tab[base*max_points + indice] = new Matrice(npoints, npoints) ;
	Matrice& resu = *tab[base*max_points + indice] ;
	resu.set_etat_qcq() ;

	double* vect = new double[npoints] ;
	for (int i=0; i<npoints; i++) {
	    for (int j=0; j<npoints; j++)
		vect[j] = 0. ;
	    vect[i] = 1. ;
	    if (base == R_JACO02) {
	    sxpun_1d(npoints, &vect, base << TRA_R) ;
	    for (int j=0; j<npoints; j++)  
		resu.set(j,i)=vect[j];
	    }
	    else {
	    sx_1d(npoints, &vect, base << TRA_R) ;
	    for (int j=0; j<npoints; j++)
		resu.set(j,i) = vect[j] ;
	    }
	}
	delete [] vect ;
    }
    return *tab[base*max_points + indice] ;
}

ostream& Diff_sx::operator>>(ostream& ost) const {

    ost << " Division by x " << endl ;

    return ost ;

}
}

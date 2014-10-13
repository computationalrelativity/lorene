/*
 *  Methods for the class Diff_xdsdx
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

char diff_xdsdx_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:52:51  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:05  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2007/12/11 15:28:11  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.2  2005/01/11 15:16:10  j_novak
 * More Diff operators.
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
void xpundsdx_1d(int, double**, int) ;

namespace {
    int nap = 0 ;
    Matrice* tab[MAX_BASE*Diff::max_points] ;
    int nr_done[Diff::max_points] ;
}

Diff_xdsdx::Diff_xdsdx(int base_r, int nr) : Diff(base_r, nr) {
    initialize() ;
}

Diff_xdsdx::Diff_xdsdx(const Diff_xdsdx& diff_in) : Diff(diff_in) {
    assert (nap != 0) ;
} 

Diff_xdsdx::~Diff_xdsdx() {}

void Diff_xdsdx::initialize() {
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

void Diff_xdsdx::operator=(const Diff_xdsdx& diff_in) {
    assert (nap != 0) ;
    Diff::operator=(diff_in) ;

}

const Matrice& Diff_xdsdx::get_matrice() const {
    
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
	    cerr << "Diff_xdsdx::get_matrice() : no space left!!" << '\n'
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

	switch (base) {

	    case R_JACO02 : {
		xpundsdx_1d(npoints, &vect, base << TRA_R) ;
		for (int j=0 ; j<npoints; j++)
		    resu.set(j,i) = vect[j] ;
		}
		break ;
	
	    case R_CHEBU : {
		dsdx_1d(npoints, &vect, R_CHEBU) ;
		mult_xm1_1d_cheb(npoints, vect, cres) ;
		for (int j=0; j<npoints; j++)
		    resu.set(j,i) = cres[j] ;
		}
		break ;

	    default : {
		xdsdx_1d(npoints, &vect, base << TRA_R) ;
		for (int j=0; j<npoints; j++)
		    resu.set(j,i) = vect[j] ;
		}
		break ;
	    }	
	}
	delete [] vect ;
	delete [] cres ;
    }
	
    return *tab[base*max_points + indice] ;
}

ostream& Diff_xdsdx::operator>>(ostream& ost) const {

    ost << " xi * d / dx " << endl ;
    return ost ;

}
}

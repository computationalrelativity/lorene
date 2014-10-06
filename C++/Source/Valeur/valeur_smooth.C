/*
 * Method of the class Valeur to make a function be smooth
 *  between the nucleus and the first shell.
 *
 * (see file valeur.h for the documentation).
 */

/*
 *   Copyright (c) 2000-2001 Keisuke Taniguchi
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


char valeur_smooth_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:13:24  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2002/10/16 14:37:16  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.4  2001/10/10  13:56:07  eric
 * Modif Joachim: pow(-1,i) --> pow(-1.,i).
 *
 * Revision 1.3  2001/01/16  16:11:41  keisuke
 * Correct the initialization of the summation
 *  and insert some explanations.
 *
 * Revision 1.2  2001/01/16  15:15:12  keisuke
 * change the argument and correct some errors.
 *
 * Revision 1.1  2001/01/16  14:54:54  keisuke
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cmath>

// Headers Lorene
#include "valeur.h"

//********************************************************************

void Valeur::smooth(int nzet, Valeur& uuva) const {

    int nucl = nzet - 1 ;
    int nr = mg->get_nr(nucl) ;
    int nt = mg->get_nt(nucl) ;
    int np = mg->get_np(nucl) ;

    // Protections
    // -----------
    assert(etat == ETATQCQ) ;
    assert(nzet > 0) ;
    assert(nr == mg->get_nr(nzet)) ;
    assert(nt == mg->get_nt(nzet)) ;
    assert(np == mg->get_np(nzet)) ;

    Valeur pot(mg) ;
    (*this).coef() ;	// the spectral coefficients are required
    pot = *((*this).c_cf) ;

    Tbl& ccf_nucl = *((pot.c_cf)->t[nucl]) ;
    Tbl& ccf_shell = *((pot.c_cf)->t[nzet]) ;


    // Get the values at the outer boundary of the nucleus
    //-----------------------------------------------------

    Tbl nucl_kj(np, nt) ;
    nucl_kj.set_etat_qcq() ;

    for (int k=0 ; k<np ; k++) {
      for (int j=0 ; j<nt ; j++) {

	double tmp = 0. ;
	for (int i=0 ; i<nr ; i++) {

	  tmp += ccf_nucl(k, j, i) ;

	}
	nucl_kj.set(k, j) = tmp ;
      }
    }


    // Get the values at the inner boundary of the first shell
    //  without the last coefficient
    //---------------------------------------------------------

    Tbl shell_kj(np, nt) ;
    shell_kj.set_etat_qcq() ;

    for (int k=0 ; k<np ; k++) {
      for (int j=0 ; j<nt ; j++) {

	double tmp2 = 0. ;
	for (int i=0 ; i<nr-1 ; i++) {

	  tmp2 += pow(-1., i) * ccf_shell(k, j, i) ;

	}
	shell_kj.set(k, j) = tmp2 ;
      }
    }


    // Set the last coefficient of the first shell
    //---------------------------------------------

    uuva.set_etat_cf_qcq() ;
    uuva.c_cf->set_etat_qcq() ;
    uuva.c_cf->t[nzet]->set_etat_qcq() ;

    Mtbl_cf& uuva_cf = *(uuva.c_cf) ;

    for (int k=0 ; k<np ; k++) {
      for (int j=0 ; j<nt ; j++) {

	uuva_cf.set(nzet, k, j, nr-1) = nucl_kj(k, j) - shell_kj(k, j) ;

      }
    }

    uuva.coef_i() ;

}

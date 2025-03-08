/*
 *  Method Mtbl_cf::display
 *
 *    (see file mtbl_cf.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003  Eric Gourgoulhon
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
 * Revision 1.5  2023/06/28 10:04:32  j_novak
 * Use of C++ strings and flows instead of C types.
 *
 * Revision 1.4  2016/12/05 16:17:59  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:08  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:15  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2003/10/19 19:51:58  e_gourgoulhon
 * First version
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>
#include <cmath>

// Lorene headers
#include "mtbl_cf.h"

namespace Lorene {
  void Mtbl_cf::display(double thres, int precis, ostream& ost) const {

    ost << "Spectral expansion (Mtbl_cf, threshold for display = " 
	<< thres << ")" << endl ; 
    ost << base << endl ; 

    if (etat == ETATNONDEF) {
      ost << "    state: UNDEFINED" << endl ;
      return ;
    }

    if (etat == ETATZERO) {
      ost << "    state: ZERO" << endl ;
      return ;
    }
	
    ost.precision(precis);
    ost.setf(ios::showpoint);
    assert(etat == ETATQCQ) ; 
    string namep ; 
    string namet ; 
    string namer ; 
	
    for (int l=0; l<nzone; l++) {

      int nr = mg->get_nr(l) ; 
      int nt = mg->get_nt(l) ; 
      int np = mg->get_np(l) ;

      ost << " --------- Domain no. " << l << " ------- nr x nt x np = "
	  << nr << " x " << nt << " x " << np << " ------" << endl ; 
      const Tbl& tcf = *(t[l]) ; 
      if (tcf.get_etat() == ETATZERO) {
	ost << "*** identically ZERO ***" << endl << endl ; 
	continue ; 
      }
      if (tcf.get_etat() == ETATNONDEF) {
	ost << "*** UNDEFINED ***" << endl << endl ; 
	continue ; 
      }
      assert( tcf.get_etat() == ETATQCQ ) ; 

      for (int k=0; k<=np; k++) {
	base.name_phi(l, k, namep) ; 
	if (namep[0] == 'u') continue ; // unused phi coefficient

	for (int j=0; j<nt; j++) {
				
	  bool test_display = false ; 
	  for (int i=0; i<nr; i++) {
	    if (fabs( tcf(k, j, i) ) >= thres) test_display = true ; 
	  }
				
	  base.name_theta(l, k, j, namet) ;
				
	  test_display = test_display && ( namet[0] != 'u' ) ;
				
	  if (test_display) {
	    ost << "# " << namep << " " << namet << " :" ;
	    for (int i=0; i<nr; i++) {
	      double cx = tcf(k, j, i) ;
	      if (fabs( cx ) >= thres) {
		base.name_r(l, k, j, i, namer) ;
		if (namer[0] == 'u') continue ; // unused r coefficient
		if ( (i>0) && (cx >= 0.) ) {
		  ost <<  " +" << setw(precis) << cx 
		      << " " << namer ; 
		}
		else {
		  ost <<  " " << setw(precis) << cx 
		      << " " << namer ; 
		}
	      }
	    }
	    ost << endl ; 	
	  }

	} // end of theta loop (index j)
			
      } // end of phi loop (index k)
		
      ost << endl ; 
		
    } // end of loop on the domains (index l)

  }
}

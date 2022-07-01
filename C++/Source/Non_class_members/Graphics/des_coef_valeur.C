/*
 * Plots the spectral coefficients of a Valeur.
 *
 * (see file graphique.h for the documentation).
 *
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


 



/*
 * $Id$
 * $Log$
 * Revision 1.7  2022/07/01 08:10:19  j_novak
 * Added a missing 'include'
 *
 * Revision 1.6  2022/02/10 16:56:57  j_novak
 * Using C++ strings to avoid warnings
 *
 * Revision 1.5  2016/12/05 16:18:06  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:21  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:16:04  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.5  2000/02/25  10:28:02  eric
 * Suppression des appels a Mtbl_cf::nettoie().
 *
 * Revision 1.4  1999/12/20  14:27:21  eric
 * Amelioration des legendes.
 *
 * Revision 1.3  1999/12/20  10:57:33  eric
 * Ajout des arguments device, newgraph, nxpage et nypage.
 *
 * Revision 1.2  1999/12/10  12:30:44  eric
 * *** empty log message ***
 *
 * Revision 1.1  1999/12/10  12:14:28  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Header C
#include <cstdlib>
#include <cstring>

// Header C++
#include <sstream>

// Header Lorene
#include "valeur.h"
#include "graphique.h"

			//-------------------------//
			//	xi coefficients	   //
			//-------------------------//

namespace Lorene {
void des_coef_xi(const Valeur& uu, int l, int k, int j, double pzero, 
		 const char* nomy, const char* title, const char* device, 
	         int newgraph, int nxpage, int nypage) {

    assert(uu.get_etat() != ETATNONDEF) ; 	
    uu.coef() ; 
    
    int nr = uu.get_mg()->get_nr(l) ; 
    
    double* cf = new double[nr] ; 
    
    // Are all the coefficients zero ?
    // -------------------------------
    if (uu.get_etat() == ETATZERO) {
	for (int i=0; i<nr; i++) {
	    cf[i] = 0 ; 
	}
    }
    else{
	assert(uu.get_etat() == ETATQCQ) ;
	for (int i=0; i<nr; i++) {
	    cf[i] = (*(uu.c_cf))(l, k, j, i) ; 
	}
    }

    const char* nomx = "i" ;

    string title1 ;
    if (title == 0x0) {
      ostringstream str_tit1 ;
      str_tit1 << "\\gc coef. for k=" << k << ", j=" << j << " (domain "
	       << l << ")" ;// << endl ;
      title1 = str_tit1.str() ;
    }
    else
      title1 = title ;
      
    string nomy1 ;
    if (nomy == 0x0) {
      ostringstream str_nomy ;
      str_nomy << "log| c\\d" << k << ',' << j << ",i\\u |" ;
      nomy1 = str_nomy.str() ;
    }
    else{
      nomy1 = nomy ;
    }

    des_coef(cf, nr, pzero, nomx, nomy1.c_str(), title1.c_str(), device, newgraph, 
	     nxpage, nypage) ;    
    
    delete [] cf ; 
    
} 

			//------------------------------//
			//	theta coefficients	//
			//------------------------------//

void des_coef_theta(const Valeur& uu, int l, int k, int i, double pzero, 
		 const char* nomy, const char* title, const char* device, 
	         int newgraph, int nxpage, int nypage) {

    assert(uu.get_etat() != ETATNONDEF) ; 	
    uu.coef() ; 
    
    int nt = uu.get_mg()->get_nt(l) ; 
    
    double* cf = new double[nt] ; 
    
    // Are all the coefficients zero ?
    // -------------------------------
    if (uu.get_etat() == ETATZERO) {
	for (int j=0; j<nt; j++) {
	    cf[j] = 0 ; 
	}
    }
    else{
	assert(uu.get_etat() == ETATQCQ) ;
	for (int j=0; j<nt; j++) {
	    cf[j] = (*(uu.c_cf))(l, k, j, i) ; 
	}
    }

    const char* nomx = "j" ; 

    string title1 ;
    if (title == 0x0) {
      ostringstream str_tit1 ;
      str_tit1 << "\\gh coef. for k=" << k << ", i=" << i << " (domain "
	       << l << ")" ;// << endl ;
      title1 = str_tit1.str() ;
    }
    else
      title1 = title ;
      
    string nomy1 ;
    if (nomy == 0x0) {
      ostringstream str_nomy ;
      str_nomy << "log| c\\d" << k << ",j," << i << "\\u |" ;
      nomy1 = str_nomy.str() ;
    }
    else{
      nomy1 = nomy ;
    }

    des_coef(cf, nt, pzero, nomx, nomy1.c_str(), title1.c_str(), device, newgraph, 
	     nxpage, nypage) ;    
    
    delete [] cf ; 
    
} 


			//------------------------------//
			//	phi coefficients	//
			//------------------------------//

void des_coef_phi(const Valeur& uu, int l, int j, int i, double pzero, 
		 const char* nomy, const char* title, const char* device, 
	         int newgraph, int nxpage, int nypage) {

    assert(uu.get_etat() != ETATNONDEF) ; 	
    uu.coef() ; 
    
    int np = uu.get_mg()->get_np(l) + 2 ; 
    
    double* cf = new double[np] ; 
    
    // Are all the coefficients zero ?
    // -------------------------------
    if (uu.get_etat() == ETATZERO) {
	for (int k=0; k<np; k++) {
	    cf[k] = 0 ; 
	}
    }
    else{
	assert(uu.get_etat() == ETATQCQ) ;
	for (int k=0; k<np; k++) {
	    cf[k] = (*(uu.c_cf))(l, k, j, i) ; 
	}
    }

    const char* nomx = "k" ; 

    string title1 ;
    if (title == 0x0) {
      ostringstream str_tit1 ;
      str_tit1 << "\\gf coef. for j=" << j << ", i=" << i << " (domain "
	       << l << ")" ;// << endl ;
      title1 = str_tit1.str() ;
    }
    else
      title1 = title ;
      
    string nomy1 ;
    if (nomy == 0x0) {
      ostringstream str_nomy ;
      str_nomy << "log| c\\dk," << j << ',' << i << "\\u |" ;
      nomy1 = str_nomy.str() ;
    }
    else{
      nomy1 = nomy ;
    }

    des_coef(cf, np, pzero, nomx, nomy1.c_str(), title1.c_str(), device, newgraph, 
	     nxpage, nypage) ;    
    
    delete [] cf ; 
    
} 
}

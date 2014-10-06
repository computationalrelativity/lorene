/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


char division_xpun_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:16:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/10/03 15:58:49  j_novak
 * Cleaning of some headers
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2000/09/07  13:19:20  phil
 * *** empty log message ***
 *
 * Revision 1.2  2000/06/06  12:23:12  phil
 * suppression des fichiers include locaux
 *
 * Revision 1.1  2000/06/06  12:18:59  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */

//standard
#include <cstdlib>

// Lorene
#include "cmp.h"
#include "proto.h"

Cmp division_xpun (const Cmp& source, int num_front) {
    
    assert (source.get_etat() != ETATNONDEF) ;
    Cmp resultat (source) ;
   
    if (resultat.get_etat() == ETATZERO)
	return resultat ;
    else  {
	resultat.va.coef() ;
	resultat.va.set_etat_cf_qcq() ;
	int base_r = source.va.base.b[num_front+1] & MSQ_R ;
    
	int nr = source.get_mp()->get_mg()->get_nr(num_front+1)  ;
	double* coef = new double[nr] ;
    
	for (int k=0 ; k<source.get_mp()->get_mg()->get_np(num_front+1)+1 ; k++)
	    if (k != 1)
		for (int j=0 ; j<source.get_mp()->get_mg()->get_nt(num_front) ; j++) {
		    for (int i=0 ; i<nr ; i++)
			coef[i] = (*resultat.va.c_cf)(num_front+1, k, j, i) ;
		    sxpun_1d (nr, &coef, base_r) ;
		    for (int i=0 ; i<nr ; i++)
			resultat.va.c_cf->set(num_front+1, k, j, i) = coef[i] ;
		}
    
	delete [] coef ;
	return resultat ;
    }
}

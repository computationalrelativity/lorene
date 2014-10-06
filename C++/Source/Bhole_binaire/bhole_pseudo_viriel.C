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


char bhole_pseudo_viriel_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/06 15:12:58  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2005/08/29 15:10:14  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.2  2003/10/03 15:58:44  j_novak
 * Cleaning of some headers
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2000/12/15  13:52:38  phil
 * modification critere sans les derivees
 *
 * Revision 2.5  2000/12/14  14:09:11  phil
 * simplification
 *
 * Revision 2.4  2000/12/14  10:45:30  phil
 * ATTENTION : PASSAGE DE PHI A PSI
 *
 * Revision 2.3  2000/11/15  18:27:08  phil
 * retour ancienne version (avec signe)
 *
 * Revision 2.2  2000/11/15  15:13:44  phil
 * *** empty log message ***
 *
 * Revision 2.1  2000/11/15  09:45:04  phil
 * *** empty log message ***
 *
 * Revision 2.0  2000/11/15  09:43:49  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

#include <cmath>

// Lorene
#include "tenseur.h"
#include "bhole.h"

double Bhole::viriel_seul () const{
    
    int nz = mp.get_mg()->get_nzone() ;
	    
    Valeur** devel_psi (psi_auto().asymptot(1)) ;
    Valeur** devel_n (n_auto().asymptot(1)) ;
    
    double erreur = (2*(*devel_psi[1])(nz-1, 0, 0, 0)
	+ (*devel_n[1])(nz-1, 0, 0, 0))/fabs ((*devel_n[1])(nz-1, 0, 0, 0)) ;
    
   return erreur ;
}


double Bhole_binaire::viriel () const{
    
    int nz_un = hole1.mp.get_mg()->get_nzone() ;
    int nz_deux = hole2.mp.get_mg()->get_nzone() ;
    
    Valeur** devel_psi_un (hole1.psi_auto().asymptot(1)) ;
    Valeur** devel_psi_deux (hole2.psi_auto().asymptot(1)) ;
    Valeur** devel_n_un (hole1.n_auto().asymptot(1)) ;
    Valeur** devel_n_deux (hole2.n_auto().asymptot(1)) ;
    
    double res = 
	(2*(*devel_psi_un[1])(nz_un-1, 0, 0, 0)+
	    2*(*devel_psi_deux[1])(nz_deux-1, 0, 0, 0)+
	    (*devel_n_deux[1])(nz_deux-1, 0, 0, 0) +
	    (*devel_n_un[1])(nz_un-1, 0, 0, 0))
	/ fabs ((*devel_n_deux[1])(nz_deux-1, 0, 0, 0) +
	    (*devel_n_un[1])(nz_un-1, 0, 0, 0)) ;
    
    return res ;
}


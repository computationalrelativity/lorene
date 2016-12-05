/*
 *   Copyright (c) 2005 Francois Limousin
 *                      Jose Luis Jaramillo
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
 * Revision 1.6  2016/12/05 16:17:46  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:52:42  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:01  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2007/04/13 15:28:55  f_limousin
 * Lots of improvements, generalisation to an arbitrary state of
 * rotation, implementation of the spatial metric given by Samaya.
 *
 * Revision 1.2  2005/04/08 12:35:07  f_limousin
 * Just to avoid warnings...
 *
 * Revision 1.1  2005/02/11 18:22:06  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

#include <cmath>

// Lorene
#include "tensor.h"
#include "isol_hor.h"

namespace Lorene {
double Single_hor::viriel_seul () const{
    
    int nz1 = mp.get_mg()->get_nzone() ;
	    
    Valeur** devel_psi (psi_auto.asymptot(1)) ;
    Valeur** devel_n (n_auto.asymptot(1)) ;
    
    double erreur = (2*(*devel_psi[1])(nz1-1, 0, 0, 0)
	+ (*devel_n[1])(nz1-1, 0, 0, 0))/fabs ((*devel_n[1])(nz1-1, 0, 0, 0)) ;
    
   return erreur ;
}


double Bin_hor::viriel () const{
    
    int nz_un = hole1.mp.get_mg()->get_nzone() ;
    int nz_deux = hole2.mp.get_mg()->get_nzone() ;
    
    Valeur** devel_psi_un (hole1.psi_auto.asymptot(1)) ;
    Valeur** devel_psi_deux (hole2.psi_auto.asymptot(1)) ;
    Valeur** devel_n_un (hole1.n_auto.asymptot(1)) ;
    Valeur** devel_n_deux (hole2.n_auto.asymptot(1)) ;
    
    double res = 
	(2*(*devel_psi_un[1])(nz_un-1, 0, 0, 0)+
	    2*(*devel_psi_deux[1])(nz_deux-1, 0, 0, 0)+
	    (*devel_n_deux[1])(nz_deux-1, 0, 0, 0) +
	    (*devel_n_un[1])(nz_un-1, 0, 0, 0))
	/ fabs ((*devel_n_deux[1])(nz_deux-1, 0, 0, 0) +
	    (*devel_n_un[1])(nz_un-1, 0, 0, 0)) ;
    
    return res ;
}

}

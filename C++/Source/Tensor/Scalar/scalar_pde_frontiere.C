/*
 * Methods Scalar::poisson_*
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 2000-2001 Philippe Grandclement (for preceding Cmp version)
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


char scalar_pde_frontiere_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2007/06/21 20:01:32  k_taniguchi
 * Modification of the method to convert Scalar to Cmp.
 *
 * Revision 1.4  2004/11/23 12:46:57  f_limousin
 * Intoduce function poisson_dir_neu(...) to solve a scalar poisson
 * equation with a mixed boundary condition (Dirichlet + Neumann).
 *
 * Revision 1.3  2003/10/03 15:58:52  j_novak
 * Cleaning of some headers
 *
 * Revision 1.2  2003/09/25 08:57:27  e_gourgoulhon
 * modif comments
 *
 * Revision 1.1  2003/09/25 08:06:56  e_gourgoulhon
 * First versions (use Cmp as intermediate quantities).
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "tensor.h"
#include "cmp.h"


Scalar Scalar::poisson_dirichlet(const Valeur& limite, int num_front) const {
    
    //    Cmp csource(*this) ; 
    Cmp csource(mp) ;
    csource = (*this).va ;
    csource.set_dzpuis((*this).get_dzpuis()) ;
    (csource.va).set_base( ((*this).va).get_base() ) ;

    Cmp cresu(mp) ;
	
    mp->poisson_frontiere(csource, limite, 1, num_front, cresu) ; 
	
    Scalar resu(cresu) ; 
    return resu ;          
}

Scalar Scalar::poisson_dir_neu(const Valeur& limite, int num_front, 
	        double fact_dir, double fact_neu) const {
    
    //    Cmp csource(*this) ; 
    Cmp csource(mp) ;
    csource = (*this).va ;
    csource.set_dzpuis((*this).get_dzpuis()) ;
    (csource.va).set_base( ((*this).va).get_base() ) ;

    Cmp cresu(mp) ;
 	
    mp->poisson_frontiere(csource, limite, 3, num_front, cresu, fact_dir,
			  fact_neu) ; 
	
    Scalar resu(cresu) ; 
    return resu ;          
}

Scalar Scalar::poisson_neumann(const Valeur& limite, int num_front) const {
    
    //    Cmp csource(*this) ; 
    Cmp csource(mp) ;
    csource = (*this).va ;
    csource.set_dzpuis((*this).get_dzpuis()) ;
    (csource.va).set_base( ((*this).va).get_base() ) ;

    Cmp cresu(mp) ;

    mp->poisson_frontiere (csource, limite, 2, num_front, cresu) ; 

    Scalar resu(cresu) ; 
    return resu ;    
}


Scalar Scalar::poisson_frontiere_double (const Valeur& lim_func, 
			 const Valeur& lim_der, int num_zone) const {
	
    Cmp csource(*this) ; 
    Cmp cresu(mp) ;

    mp->poisson_frontiere_double(csource, lim_func, lim_der, num_zone, cresu) ; 

    Scalar resu(cresu) ; 
    return resu ;    
}		




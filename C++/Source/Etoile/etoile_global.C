/*
 * Methods of class Etoile to compute global quantities
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char etoile_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2005/01/18 22:37:30  k_taniguchi
 * Modify the method of ray_eq(int kk).
 *
 * Revision 1.4  2005/01/18 20:35:46  k_taniguchi
 * Addition of ray_eq(int kk).
 *
 * Revision 1.3  2005/01/17 20:40:56  k_taniguchi
 * Addition of ray_eq_3pis2().
 *
 * Revision 1.2  2003/12/05 14:50:26  j_novak
 * To suppress some warnings...
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  2000/07/21  12:02:14  eric
 * Etoile::ray_eq_pi() : traitement du cas type_p = SYM.
 *
 * Revision 1.1  2000/01/28  17:18:45  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "etoile.h"

			//--------------------------//
			//	Stellar surface	    //
			//--------------------------//

const Itbl& Etoile::l_surf() const {

    if (p_l_surf == 0x0) {    // a new computation is required
    
	assert(p_xi_surf == 0x0) ;  // consistency check
	
	int np = mp.get_mg()->get_np(0) ;   
	int nt = mp.get_mg()->get_nt(0) ;   
	
	p_l_surf = new Itbl(np, nt) ;
	p_xi_surf = new Tbl(np, nt) ;
	
	double ent0 = 0 ;	// definition of the surface
	double precis = 1.e-15 ; 
	int nitermax = 100 ; 
	int niter ; 
	
	(ent().va).equipot(ent0, nzet, precis, nitermax, niter, *p_l_surf, 
		    *p_xi_surf) ; 
    
    }
   
    return *p_l_surf ; 
    
}

const Tbl& Etoile::xi_surf() const {

    if (p_xi_surf == 0x0) {    // a new computation is required
    
	assert(p_l_surf == 0x0) ;  // consistency check
	
	l_surf() ;  // the computation is done by l_surf()
    
    }
   
    return *p_xi_surf ; 
    
}


			//--------------------------//
			//	Coordinate radii    //
			//--------------------------//

double Etoile::ray_eq() const {

    if (p_ray_eq == 0x0) {    // a new computation is required
	
	const Mg3d& mg = *(mp.get_mg()) ;
	
	int type_t = mg.get_type_t() ; 
#ifndef NDEBUG
	int type_p = mg.get_type_p() ; 
#endif
	int nt = mg.get_nt(0) ; 	
	
	if ( type_t == SYM ) {
	    assert( (type_p == SYM) || (type_p == NONSYM) ) ; 
	    int k = 0 ; 
	    int j = nt-1 ; 
	    int l = l_surf()(k, j) ; 
	    double xi = xi_surf()(k, j) ; 
	    double theta = M_PI / 2 ; 
	    double phi = 0 ; 
	    
	    p_ray_eq = new double( mp.val_r(l, xi, theta, phi) ) ;

	}
	else {
	    cout << "Etoile::ray_eq : the case type_t = " << type_t
		 << " is not contemplated yet !" << endl ;
	    abort() ; 
	}

    }
    
    return *p_ray_eq ; 

} 


double Etoile::ray_eq_pis2() const {

    if (p_ray_eq_pis2 == 0x0) {    // a new computation is required
	
	const Mg3d& mg = *(mp.get_mg()) ;
	
	int type_t = mg.get_type_t() ; 
	int type_p = mg.get_type_p() ; 
	int nt = mg.get_nt(0) ; 	
	int np = mg.get_np(0) ; 	
	
	if ( type_t == SYM ) {
	
	    int j = nt-1 ; 
	    double theta = M_PI / 2 ; 
	    double phi = M_PI / 2 ;
	    
	    switch (type_p) {
	    
		case SYM : {
		    int k = np / 2  ; 
		    int l = l_surf()(k, j) ; 
		    double xi = xi_surf()(k, j) ; 
		    p_ray_eq_pis2 = new double( mp.val_r(l, xi, theta, phi) ) ;
		    break ; 
		}
	    
		case NONSYM : {
		    assert( np % 4 == 0 ) ; 
		    int k = np / 4  ; 
		    int l = l_surf()(k, j) ; 
		    double xi = xi_surf()(k, j) ; 
		    p_ray_eq_pis2 = new double( mp.val_r(l, xi, theta, phi) ) ;
		    break ; 
		}
	    
		default : {
		    cout << "Etoile::ray_eq_pis2 : the case type_p = " 
			<< type_p << " is not contemplated yet !" << endl ;
		    abort() ; 
		}
	    } 

	}
	else {
	    cout << "Etoile::ray_eq_pis2 : the case type_t = " << type_t
	         << " is not contemplated yet !" << endl ;
	    abort() ; 
	}

    }
    
    return *p_ray_eq_pis2 ; 

} 


double Etoile::ray_eq_pi() const {

    if (p_ray_eq_pi == 0x0) {    // a new computation is required
	
	const Mg3d& mg = *(mp.get_mg()) ;
	
	int type_t = mg.get_type_t() ; 
	int type_p = mg.get_type_p() ; 
	int nt = mg.get_nt(0) ; 	
	int np = mg.get_np(0) ; 	
	
	if ( type_t == SYM ) {

	    switch (type_p) {
		
		case SYM : {
		    p_ray_eq_pi = new double( ray_eq() ) ;
		    break ; 
		}		
		
		case NONSYM : {
		    int k = np / 2  ; 
		    int j = nt-1 ; 
		    int l = l_surf()(k, j) ; 
		    double xi = xi_surf()(k, j) ; 
		    double theta = M_PI / 2 ; 
		    double phi = M_PI ; 
	    
		    p_ray_eq_pi = new double( mp.val_r(l, xi, theta, phi) ) ;
		    break ;
		}
		
		default : {

	    cout << "Etoile::ray_eq_pi : the case type_t = " << type_t
		 << " and type_p = " << type_p << endl ; 
	    cout << " is not contemplated yet !" << endl ;
	    abort() ; 
	    break ; 
		}
	    }
	}

    }
    
    return *p_ray_eq_pi ; 

} 

double Etoile::ray_eq_3pis2() const {

    if (p_ray_eq_3pis2 == 0x0) {    // a new computation is required
	
	const Mg3d& mg = *(mp.get_mg()) ;
	
	int type_t = mg.get_type_t() ; 
	int type_p = mg.get_type_p() ; 
	int nt = mg.get_nt(0) ; 	
	int np = mg.get_np(0) ; 	
	
	if ( type_t == SYM ) {
	
	    int j = nt-1 ; 
	    double theta = M_PI / 2 ; 
	    double phi = 3. * M_PI / 2 ;
	    
	    switch (type_p) {
	    
		case SYM : {
		    p_ray_eq_3pis2 = new double( ray_eq_pis2() ) ;
		    break ; 
		}
	    
		case NONSYM : {
		    assert( np % 4 == 0 ) ; 
		    int k = 3 * np / 4  ; 
		    int l = l_surf()(k, j) ; 
		    double xi = xi_surf()(k, j) ; 
		    p_ray_eq_3pis2 = new double( mp.val_r(l,xi,theta,phi) ) ;
		    break ; 
		}
	    
		default : {
		    cout << "Etoile::ray_eq_3pis2 : the case type_p = " 
			<< type_p << " is not contemplated yet !" << endl ;
		    abort() ; 
		}
	    } 

	}
	else {
	    cout << "Etoile::ray_eq_3pis2 : the case type_t = " << type_t
	         << " is not contemplated yet !" << endl ;
	    abort() ; 
	}

    }
    
    return *p_ray_eq_3pis2 ; 

} 

double Etoile::ray_pole() const {

    if (p_ray_pole == 0x0) {    // a new computation is required
	
#ifndef NDEBUG
	const Mg3d& mg = *(mp.get_mg()) ;
	int type_t = mg.get_type_t() ; 
#endif	
	assert( (type_t == SYM) || (type_t == NONSYM) ) ;  
	
	int k = 0 ; 
	int j = 0 ; 
	int l = l_surf()(k, j) ; 
	double xi = xi_surf()(k, j) ; 
	double theta = 0 ; 
	double phi = 0 ; 
	    
	p_ray_pole = new double( mp.val_r(l, xi, theta, phi) ) ;

    }
    
    return *p_ray_pole ; 

} 

double Etoile::ray_eq(int kk) const {

    const Mg3d& mg = *(mp.get_mg()) ;
	
    int type_t = mg.get_type_t() ; 
    int type_p = mg.get_type_p() ; 
    int nt = mg.get_nt(0) ; 	
    int np = mg.get_np(0) ; 

    assert( kk >= 0 ) ;
    assert( kk < np ) ;
	
    double ray_eq_kk ;
    if ( type_t == SYM ) {
	
        int j = nt-1 ; 
	double theta = M_PI / 2 ; 
	    
	switch (type_p) {
	    
	    case SYM : {
	        cout << "Etoile::ray_eq(kk) : the case type_p = " 
		     << type_p << " is not contemplated yet !" << endl ;
		abort() ; 
	    }
	    
	    case NONSYM : {
	        double phi = 2. * kk * M_PI / np ;
		int l = l_surf()(kk, j) ; 
		double xi = xi_surf()(kk, j) ; 
		ray_eq_kk = mp.val_r(l,xi,theta,phi) ;
		break ; 
	    }
	    
	    default : {
	        cout << "Etoile::ray_eq(kk) : the case type_p = " 
		     << type_p << " is not contemplated yet !" << endl ;
		abort() ; 
	    }
	} 

    }
    else {
        cout << "Etoile::ray_eq(kk) : the case type_t = " << type_t
	     << " is not contemplated yet !" << endl ;
	abort() ; 
    }

    return ray_eq_kk ;
}


			//--------------------------//
			//	Baryon mass	    //
			//--------------------------//

double Etoile::mass_b() const {

    if (p_mass_b == 0x0) {    // a new computation is required
	
	if (relativistic) {
	    cout << 
	    "Etoile::mass_b : in the relativistic case, the baryon mass"
	    << endl << 
	    "computation cannot be performed by the base class Etoile !" 
	    << endl ; 
	    abort() ; 
	}
	
	assert(nbar.get_etat() == ETATQCQ) ; 
	p_mass_b = new double( nbar().integrale() ) ;

    }
    
    return *p_mass_b ; 

} 
		
			//----------------------------//
			//	Gravitational mass    //
			//----------------------------//

double Etoile::mass_g() const {

    if (p_mass_g == 0x0) {    // a new computation is required
	
	if (relativistic) {
	    cout << 
	    "Etoile::mass_g : in the relativistic case, the gravitational mass"
	    << endl << 
	    "computation cannot be performed by the base class Etoile !" 
	    << endl ; 
	    abort() ; 
	}
	
	p_mass_g = new double( mass_b() ) ;   // in the Newtonian case
					      //  M_g = M_b

    }
    
    return *p_mass_g ; 

} 
		
